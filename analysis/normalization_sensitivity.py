"""Reviewer-requested propagation-normalization and degree-correction audit.

The locked primary score is not replaced here. This module asks how the
ACC-context component changes under three alternatives to the original
column-stochastic/min-max implementation:

1. gene-level rank normalization of the column-stochastic RWR;
2. division by the uniform-restart steady-state distribution, followed by
   gene-level rank normalization;
3. symmetric-normalized network propagation, followed by gene-level rank
   normalization.

Every branch uses the same ACC seeds, STRING graph, restart probability,
degree-matched seed sets, 108-drug universe, pseudo-count and BH family.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import scipy
from scipy import sparse, stats

from analysis.acc_primary_pipeline import compute_primary_analysis, load_inputs
from analysis.mechanism_enrichment import benjamini_hochberg
from analysis.method_strengthening import (
    C_ACC_PSEUDO_COUNT,
    CDK46_DRUGS,
    NULL_DRAWS,
    REQUESTED_DEGREE_BINS,
    RNG_SEED,
    RWR_RESTART,
    STRING_THRESHOLD,
    build_association_matrix,
    build_restart_matrix,
    compute_c_acc_matrix,
    empirical_upper_p,
    generate_degree_matched_seed_sets,
    load_disease_seed_weights,
    load_string_graph,
    minmax_columns,
    percentile_columns,
    random_walk_with_restart,
    weighted_pagerank,
)


ANALYSIS_VERSION = "normalization-sensitivity-v1-reviewer-m4"
VARIANT_DESCRIPTIONS = {
    "column_minmax": (
        "Locked column-stochastic RWR with gene-level min-max scaling"
    ),
    "column_gene_rank": (
        "Column-stochastic RWR with gene-level average-rank percentiles"
    ),
    "uniform_ratio_gene_rank": (
        "Column-stochastic RWR divided by the uniform-restart steady state, "
        "then gene-level average-rank percentiles"
    ),
    "symmetric_gene_rank": (
        "Symmetric D^-1/2 W D^-1/2 propagation with gene-level "
        "average-rank percentiles"
    ),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: f"{value:.12g}" if isinstance(value, float) else value
                    for key, value in row.items()
                }
            )


def symmetric_normalized_operator(
    adjacency: sparse.csr_matrix,
) -> sparse.csr_matrix:
    """Return D^-1/2 W D^-1/2 using weighted node strength."""

    strength = np.asarray(adjacency.sum(axis=1)).ravel()
    inverse_sqrt = np.zeros_like(strength, dtype=float)
    nonzero = strength > 0
    inverse_sqrt[nonzero] = 1.0 / np.sqrt(strength[nonzero])
    diagonal = sparse.diags(inverse_sqrt)
    operator = (diagonal @ adjacency @ diagonal).tocsr()
    if operator.shape != adjacency.shape:
        raise AssertionError("Symmetric operator shape changed unexpectedly")
    return operator


def network_smooth_with_restart(
    operator: sparse.csr_matrix,
    restart: np.ndarray,
    alpha: float = RWR_RESTART,
    tolerance: float = 1e-10,
    max_iterations: int = 500,
) -> tuple[np.ndarray, int, float]:
    """Solve f=(1-alpha)Sf+alpha*r for a nonnegative symmetric operator."""

    restart_array = np.asarray(restart, dtype=float)
    if restart_array.ndim == 1:
        restart_array = restart_array[:, None]
    if restart_array.shape[0] != operator.shape[0]:
        raise ValueError("Restart matrix and propagation operator do not align")
    if not np.isfinite(restart_array).all() or np.any(restart_array < 0):
        raise ValueError("Restart matrix must contain finite nonnegative values")
    if not 0 < alpha <= 1:
        raise ValueError("Restart probability must lie in (0, 1]")

    propagated = restart_array.copy()
    final_delta = math.inf
    for iteration in range(1, max_iterations + 1):
        updated = (1.0 - alpha) * (operator @ propagated) + alpha * restart_array
        final_delta = float(
            np.max(np.sum(np.abs(updated - propagated), axis=0))
        )
        propagated = updated
        if final_delta < tolerance:
            break
    else:
        raise RuntimeError(
            "Symmetric propagation failed to converge; "
            f"last max L1 delta={final_delta}"
        )
    if not np.isfinite(propagated).all() or np.any(propagated < -1e-14):
        raise ValueError("Symmetric propagation produced invalid values")
    return propagated, iteration, final_delta


def top_k_jaccard(
    names: Sequence[str],
    first: np.ndarray,
    second: np.ndarray,
    k: int,
) -> float:
    if not 1 <= k <= len(names):
        raise ValueError("k must lie within the score vector")

    def top(values: np.ndarray) -> set[str]:
        return {
            name
            for name, _ in sorted(
                zip(names, np.asarray(values, dtype=float), strict=True),
                key=lambda item: (-item[1], item[0]),
            )[:k]
        }

    first_top = top(first)
    second_top = top(second)
    return len(first_top & second_top) / len(first_top | second_top)


def _gene_variants(
    column_values: np.ndarray,
    symmetric_values: np.ndarray,
    uniform_column_values: np.ndarray,
) -> dict[str, np.ndarray]:
    column = np.asarray(column_values, dtype=float)
    symmetric = np.asarray(symmetric_values, dtype=float)
    uniform = np.asarray(uniform_column_values, dtype=float)
    if column.ndim == 1:
        column = column[:, None]
    if symmetric.ndim == 1:
        symmetric = symmetric[:, None]
    if uniform.ndim == 1:
        uniform = uniform[:, None]
    if not (column.shape == symmetric.shape):
        raise ValueError("Column and symmetric propagated matrices do not align")
    if uniform.shape[0] != column.shape[0] or uniform.shape[1] != 1:
        raise ValueError("Uniform steady state must be a single aligned column")
    if np.any(uniform <= 0):
        raise ValueError("Uniform-restart steady state contains a zero")

    ratio = column / uniform
    return {
        "column_minmax": minmax_columns(column),
        "column_gene_rank": percentile_columns(column),
        "uniform_ratio_gene_rank": percentile_columns(ratio),
        "symmetric_gene_rank": percentile_columns(symmetric),
    }


def _drug_percentiles(
    association_matrix: sparse.csr_matrix,
    primary_indices: np.ndarray,
    gene_values: np.ndarray,
) -> np.ndarray:
    c_acc_all = compute_c_acc_matrix(
        association_matrix,
        gene_values,
        pseudo_count=C_ACC_PSEUDO_COUNT,
    )
    return percentile_columns(c_acc_all[primary_indices, :])


def run_analysis(
    project_root: Path,
    n_null: int = NULL_DRAWS,
    null_batch_size: int = 64,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    root = project_root.resolve()
    inputs = load_inputs(root)
    primary = compute_primary_analysis(inputs)
    primary_rows = list(primary.primary_rows)
    primary_drugs = tuple(str(row["drug"]) for row in primary_rows)
    all_drugs = tuple(sorted(inputs.associations))
    associated_genes = tuple(sorted(set().union(*inputs.associations.values())))
    seed_weights = load_disease_seed_weights(
        root / "data" / "ACC_P0.5C_gene_weights_v1.csv"
    )
    graph = load_string_graph(
        root / "9606.protein.info.v12.0.txt.gz",
        root / "9606.protein.links.v12.0.txt.gz",
        required_nodes=set(associated_genes) | set(seed_weights),
        threshold=STRING_THRESHOLD,
    )
    symmetric_operator = symmetric_normalized_operator(graph.adjacency)
    gene_indices = np.asarray(
        [graph.node_index[gene] for gene in associated_genes],
        dtype=int,
    )
    association_matrix = build_association_matrix(
        all_drugs,
        associated_genes,
        inputs.associations,
    )
    all_drug_index = {drug: index for index, drug in enumerate(all_drugs)}
    primary_indices = np.asarray(
        [all_drug_index[drug] for drug in primary_drugs],
        dtype=int,
    )

    observed_restart = build_restart_matrix(graph, (seed_weights,))
    observed_column, column_iterations, column_delta = (
        random_walk_with_restart(graph, observed_restart)
    )
    observed_symmetric, symmetric_iterations, symmetric_delta = (
        network_smooth_with_restart(symmetric_operator, observed_restart)
    )
    uniform_restart = np.full(
        (len(graph.node_names), 1),
        1.0 / len(graph.node_names),
        dtype=float,
    )
    uniform_column, uniform_iterations, uniform_delta = (
        random_walk_with_restart(graph, uniform_restart)
    )
    observed_gene_variants = _gene_variants(
        observed_column[gene_indices, :],
        observed_symmetric[gene_indices, :],
        uniform_column[gene_indices, :],
    )
    observed_drug_variants = {
        name: _drug_percentiles(
            association_matrix,
            primary_indices,
            values,
        )[:, 0]
        for name, values in observed_gene_variants.items()
    }
    locked_c_acc = np.asarray(
        [float(row["C_ACC_pct"]) for row in primary_rows],
        dtype=float,
    )
    if stats.spearmanr(
        observed_drug_variants["column_minmax"],
        locked_c_acc,
    ).statistic < 0.999999:
        raise ValueError("Locked C_ACC component was not reproduced")

    matched_rows, degree_bin_edges = generate_degree_matched_seed_sets(
        seed_weights,
        graph.node_names,
        graph.degree,
        n_draws=n_null,
        rng_seed=RNG_SEED,
        n_bins=REQUESTED_DEGREE_BINS,
    )
    rows_by_replicate: dict[int, list[dict[str, Any]]] = {}
    for row in matched_rows:
        rows_by_replicate.setdefault(int(row["replicate"]), []).append(row)
    null_drug_variants = {
        name: np.empty((len(primary_drugs), n_null), dtype=float)
        for name in VARIANT_DESCRIPTIONS
    }
    column_null_iterations: list[int] = []
    symmetric_null_iterations: list[int] = []
    column_null_deltas: list[float] = []
    symmetric_null_deltas: list[float] = []
    for batch_start in range(0, n_null, null_batch_size):
        batch_replicates = list(
            range(batch_start, min(batch_start + null_batch_size, n_null))
        )
        batch_weights = [
            {
                row["null_seed"]: float(row["weight"])
                for row in rows_by_replicate[replicate]
            }
            for replicate in batch_replicates
        ]
        restart = build_restart_matrix(graph, batch_weights)
        column, column_iteration, column_final_delta = (
            random_walk_with_restart(graph, restart)
        )
        symmetric, symmetric_iteration, symmetric_final_delta = (
            network_smooth_with_restart(symmetric_operator, restart)
        )
        column_null_iterations.extend(
            [column_iteration] * len(batch_replicates)
        )
        symmetric_null_iterations.extend(
            [symmetric_iteration] * len(batch_replicates)
        )
        column_null_deltas.extend(
            [column_final_delta] * len(batch_replicates)
        )
        symmetric_null_deltas.extend(
            [symmetric_final_delta] * len(batch_replicates)
        )
        gene_variants = _gene_variants(
            column[gene_indices, :],
            symmetric[gene_indices, :],
            uniform_column[gene_indices, :],
        )
        for name, values in gene_variants.items():
            null_drug_variants[name][:, batch_replicates] = (
                _drug_percentiles(
                    association_matrix,
                    primary_indices,
                    values,
                )
            )

    pagerank, _, _ = weighted_pagerank(graph)
    gene_degree = graph.degree[gene_indices]
    gene_strength = graph.strength[gene_indices]
    gene_pagerank = pagerank[gene_indices]
    primary_index = {drug: index for index, drug in enumerate(primary_drugs)}
    cdk_indices = np.asarray(
        [primary_index[drug] for drug in CDK46_DRUGS],
        dtype=int,
    )
    summary_rows: list[dict[str, Any]] = []
    drug_rows: list[dict[str, Any]] = []
    for name, description in VARIANT_DESCRIPTIONS.items():
        observed_gene = observed_gene_variants[name][:, 0]
        observed_drug = observed_drug_variants[name]
        null_matrix = null_drug_variants[name]
        null_mean = null_matrix.mean(axis=1)
        null_sd = null_matrix.std(axis=1, ddof=1)
        empirical_p = np.asarray(
            [
                empirical_upper_p(observed_drug[index], null_matrix[index])
                for index in range(len(primary_drugs))
            ],
            dtype=float,
        )
        q_values = np.asarray(
            benjamini_hochberg(empirical_p.tolist()),
            dtype=float,
        )
        z_values = np.divide(
            observed_drug - null_mean,
            null_sd,
            out=np.full_like(null_sd, np.nan),
            where=null_sd > 0,
        )
        observed_cdk_mean = float(observed_drug[cdk_indices].mean())
        null_cdk_mean = null_matrix[cdk_indices, :].mean(axis=0)
        cdk_p = empirical_upper_p(observed_cdk_mean, null_cdk_mean)
        summary_rows.append(
            {
                "variant": name,
                "description": description,
                "gene_rho_degree": float(
                    stats.spearmanr(observed_gene, gene_degree).statistic
                ),
                "gene_rho_strength": float(
                    stats.spearmanr(observed_gene, gene_strength).statistic
                ),
                "gene_rho_PageRank": float(
                    stats.spearmanr(observed_gene, gene_pagerank).statistic
                ),
                "drug_rho_vs_locked_C_ACC": float(
                    stats.spearmanr(observed_drug, locked_c_acc).statistic
                ),
                "drug_top20_jaccard_vs_locked_C_ACC": top_k_jaccard(
                    primary_drugs,
                    observed_drug,
                    locked_c_acc,
                    k=20,
                ),
                "n_drugs_q_lt_0_05": int(np.count_nonzero(q_values < 0.05)),
                "minimum_q_bh_108": float(q_values.min()),
                "CDK46_observed_mean_C_ACC_pct": observed_cdk_mean,
                "CDK46_degree_matched_empirical_p": cdk_p,
            }
        )
        for index, drug in enumerate(primary_drugs):
            drug_rows.append(
                {
                    "variant": name,
                    "drug": drug,
                    "observed_C_ACC_pct": float(observed_drug[index]),
                    "null_mean_C_ACC_pct": float(null_mean[index]),
                    "null_sd_C_ACC_pct": float(null_sd[index]),
                    "z_degree_matched": float(z_values[index]),
                    "empirical_p_upper": float(empirical_p[index]),
                    "q_bh_108": float(q_values[index]),
                    "null_draws": n_null,
                }
            )

    cdk_q_values = benjamini_hochberg(
        [
            float(row["CDK46_degree_matched_empirical_p"])
            for row in summary_rows
        ]
    )
    for row, q_value in zip(summary_rows, cdk_q_values, strict=True):
        row["CDK46_q_bh_across_variants"] = float(q_value)

    metrics = {
        "analysis_version": ANALYSIS_VERSION,
        "null_draws": n_null,
        "null_rng_seed": RNG_SEED,
        "null_batch_size": null_batch_size,
        "empirical_p_minimum_resolution": 1.0 / (n_null + 1.0),
        "BH_q_minimum_possible": len(primary_drugs) / (n_null + 1.0),
        "BH_resolution_adequate": (
            len(primary_drugs) / (n_null + 1.0) < 0.05
        ),
        "primary_universe_n": len(primary_drugs),
        "associated_gene_n": len(associated_genes),
        "ACC_seed_n": len(seed_weights),
        "degree_bins_effective": len(degree_bin_edges) - 1,
        "column_observed_iterations": column_iterations,
        "column_observed_final_delta": column_delta,
        "uniform_observed_iterations": uniform_iterations,
        "uniform_observed_final_delta": uniform_delta,
        "symmetric_observed_iterations": symmetric_iterations,
        "symmetric_observed_final_delta": symmetric_delta,
        "column_null_iterations_min": min(column_null_iterations),
        "column_null_iterations_max": max(column_null_iterations),
        "column_null_final_delta_max": max(column_null_deltas),
        "symmetric_null_iterations_min": min(symmetric_null_iterations),
        "symmetric_null_iterations_max": max(symmetric_null_iterations),
        "symmetric_null_final_delta_max": max(symmetric_null_deltas),
        "variants": {
            row["variant"]: {
                key: value
                for key, value in row.items()
                if key not in {"variant", "description"}
            }
            for row in summary_rows
        },
    }
    drug_rows.sort(
        key=lambda row: (
            row["variant"],
            row["empirical_p_upper"],
            row["drug"],
        )
    )
    return summary_rows, drug_rows, metrics


def write_outputs(
    project_root: Path,
    output_dir: Path,
    summary_rows: Sequence[Mapping[str, Any]],
    drug_rows: Sequence[Mapping[str, Any]],
    metrics: Mapping[str, Any],
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "summary": output_dir / "normalization_sensitivity_summary.csv",
        "drugs": output_dir / "normalization_sensitivity_null_primary108.csv",
        "metrics": output_dir / "normalization_sensitivity_metrics.json",
        "report": output_dir / "normalization_sensitivity_report.md",
        "manifest": output_dir / "run_manifest.md",
    }
    _write_csv(paths["summary"], summary_rows)
    _write_csv(paths["drugs"], drug_rows)
    paths["metrics"].write_text(
        json.dumps(metrics, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    report_lines = [
        "# Propagation-normalization sensitivity report",
        "",
        f"- Version: `{metrics['analysis_version']}`.",
        f"- Degree-matched draws per variant: {metrics['null_draws']:,}.",
        f"- Best-case BH q floor: {metrics['BH_q_minimum_possible']:.6f}.",
        "",
        "| Variant | rho(gene, degree) | rho(drug, locked C_ACC) | Top-20 Jaccard | Drugs q<0.05 | Minimum q | CDK4/6 P | CDK4/6 q across variants |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        report_lines.append(
            f"| {row['variant']} | {row['gene_rho_degree']:.3f} | "
            f"{row['drug_rho_vs_locked_C_ACC']:.3f} | "
            f"{row['drug_top20_jaccard_vs_locked_C_ACC']:.3f} | "
            f"{row['n_drugs_q_lt_0_05']} | "
            f"{row['minimum_q_bh_108']:.4f} | "
            f"{row['CDK46_degree_matched_empirical_p']:.4f} | "
            f"{row['CDK46_q_bh_across_variants']:.4f} |"
        )
    report_lines.extend(
        [
            "",
            "The alternatives are sensitivity analyses, not outcome-optimized "
            "replacement models. A drug-level FDR signal means that its network "
            "context exceeded the matched-seed null under that transformation; "
            "it does not establish efficacy or clinical validity.",
            "",
        ]
    )
    paths["report"].write_text(
        "\n".join(report_lines),
        encoding="utf-8",
    )
    input_paths = (
        project_root / "data" / "bindex_network" / "bindex_edges_1304.csv",
        project_root / "data" / "bindex_network" / "rACC_399_fullSTRING.csv",
        project_root / "data" / "ACC_P0.5C_gene_weights_v1.csv",
        project_root / "9606.protein.info.v12.0.txt.gz",
        project_root / "9606.protein.links.v12.0.txt.gz",
    )
    paths["manifest"].write_text(
        "\n".join(
            [
                "# Normalization-sensitivity run manifest",
                "",
                f"- Analysis version: `{ANALYSIS_VERSION}`",
                "- Command: `python -m analysis.normalization_sensitivity "
                "--project-root .`",
                f"- Python: `{platform.python_version()}`",
                f"- NumPy: `{np.__version__}`",
                f"- SciPy: `{scipy.__version__}`",
                f"- RNG seed: `{metrics['null_rng_seed']}`",
                f"- Null draws: `{metrics['null_draws']}`",
                "",
                "## Input SHA-256",
                "",
                *[
                    f"- `{path.relative_to(project_root).as_posix()}`: "
                    f"`{_sha256(path)}`"
                    for path in input_paths
                ],
                "",
                "## Output SHA-256",
                "",
                *[
                    f"- `{path.name}`: `{_sha256(path)}`"
                    for key, path in paths.items()
                    if key != "manifest"
                ],
                "",
            ]
        ),
        encoding="utf-8",
    )
    return paths


def run(
    project_root: Path,
    output_dir: Path | None = None,
    n_null: int = NULL_DRAWS,
    null_batch_size: int = 64,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    root = project_root.resolve()
    summary_rows, drug_rows, metrics = run_analysis(
        root,
        n_null=n_null,
        null_batch_size=null_batch_size,
    )
    target = (
        output_dir.resolve()
        if output_dir is not None
        else root / "results" / "normalization_sensitivity"
    )
    write_outputs(root, target, summary_rows, drug_rows, metrics)
    return summary_rows, drug_rows, metrics


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--null-draws", type=int, default=NULL_DRAWS)
    parser.add_argument("--null-batch-size", type=int, default=64)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    summary_rows, _, metrics = run(
        args.project_root,
        output_dir=args.output_dir,
        n_null=args.null_draws,
        null_batch_size=args.null_batch_size,
    )
    print(
        json.dumps(
            {
                "status": "ok",
                "analysis_version": ANALYSIS_VERSION,
                "null_draws": metrics["null_draws"],
                "variants": {
                    row["variant"]: {
                        "rho_degree": row["gene_rho_degree"],
                        "n_q_lt_0_05": row["n_drugs_q_lt_0_05"],
                        "CDK46_p": row[
                            "CDK46_degree_matched_empirical_p"
                        ],
                    }
                    for row in summary_rows
                },
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
