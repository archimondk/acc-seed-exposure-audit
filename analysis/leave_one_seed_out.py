"""Reproduce Protocol Amendment 4 full leave-one-seed-out influence scan."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import scipy

from analysis.acc_primary_pipeline import compute_primary_analysis, load_inputs
from analysis.dirichlet_weight_sensitivity import (
    ordinal_rank_columns,
    rank_spearman_columns,
    top_k_jaccard_columns,
)
from analysis.method_strengthening import (
    build_association_matrix,
    build_restart_matrix,
    compute_c_acc_matrix,
    load_disease_seed_weights,
    load_string_graph,
    minmax_columns,
    percentile_columns,
    random_walk_with_restart,
)
from analysis.normalization_sensitivity import (
    network_smooth_with_restart,
    symmetric_normalized_operator,
)


ANALYSIS_VERSION = "leave-one-seed-out-v1"
PROTOCOL_PATH = Path("experiments/amendment4_leave_one_seed_out_protocol_v1.md")
OUTPUT_DIR = Path("results/leave_one_seed_out")
VARIANTS = (
    "column_minmax",
    "column_gene_rank",
    "uniform_ratio_gene_rank",
    "symmetric_gene_rank",
)
TOP_K = 20
PSEUDO_COUNT = 3.0


@dataclass(frozen=True)
class LeaveOneSeedOutResult:
    rows: tuple[dict[str, Any], ...]
    variant_summary: dict[str, dict[str, Any]]
    seed_summary: tuple[dict[str, Any], ...]
    metrics: dict[str, Any]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _gene_variants(
    column_values: np.ndarray,
    symmetric_values: np.ndarray,
    uniform_column_values: np.ndarray,
) -> dict[str, np.ndarray]:
    column = np.asarray(column_values, dtype=float)
    symmetric = np.asarray(symmetric_values, dtype=float)
    uniform = np.asarray(uniform_column_values, dtype=float)
    if uniform.ndim == 1:
        uniform = uniform[:, None]
    if column.shape != symmetric.shape:
        raise ValueError("Column and symmetric propagated matrices do not align")
    if uniform.shape != (column.shape[0], 1) or np.any(uniform <= 0):
        raise ValueError("Uniform-reference steady state is invalid")
    return {
        "column_minmax": minmax_columns(column),
        "column_gene_rank": percentile_columns(column),
        "uniform_ratio_gene_rank": percentile_columns(column / uniform),
        "symmetric_gene_rank": percentile_columns(symmetric),
    }


def _rank_number(
    records: Sequence[Mapping[str, Any]],
    seed: str,
    metric: str,
) -> int:
    ordered = sorted(
        records,
        key=lambda row: (-float(row[metric]), str(row["seed"])),
    )
    return next(index for index, row in enumerate(ordered, start=1) if row["seed"] == seed)


def run_analysis(project_root: Path) -> LeaveOneSeedOutResult:
    root = project_root.resolve()
    inputs = load_inputs(root)
    primary = compute_primary_analysis(inputs)
    primary_by_drug = {
        str(row["drug"]): row for row in primary.primary_rows
    }
    drug_names = tuple(sorted(primary_by_drug))
    if len(drug_names) != 108:
        raise ValueError(f"Expected 108 primary drugs, got {len(drug_names)}")
    residual_pct = np.asarray(
        [float(primary_by_drug[drug]["residual_pct"]) for drug in drug_names],
        dtype=float,
    )

    raw_seed_weights = load_disease_seed_weights(
        root / "data" / "ACC_P0.5C_gene_weights_v1.csv"
    )
    seed_names = tuple(sorted(raw_seed_weights))
    if len(seed_names) != 45:
        raise ValueError(f"Expected 45 disease-only seeds, got {len(seed_names)}")
    normalized_seed_weights = {
        seed: raw_seed_weights[seed] / sum(raw_seed_weights.values())
        for seed in seed_names
    }
    restart_columns: list[Mapping[str, float]] = [raw_seed_weights]
    restart_columns.extend(
        {
            gene: weight
            for gene, weight in raw_seed_weights.items()
            if gene != omitted
        }
        for omitted in seed_names
    )

    all_drugs = tuple(sorted(inputs.associations))
    associated_genes = tuple(sorted(set().union(*inputs.associations.values())))
    graph = load_string_graph(
        root / "9606.protein.info.v12.0.txt.gz",
        root / "9606.protein.links.v12.0.txt.gz",
        required_nodes=set(associated_genes) | set(seed_names),
    )
    restart_matrix = build_restart_matrix(graph, restart_columns)
    column_values, column_iterations, column_delta = random_walk_with_restart(
        graph, restart_matrix
    )
    symmetric_operator = symmetric_normalized_operator(graph.adjacency)
    symmetric_values, symmetric_iterations, symmetric_delta = (
        network_smooth_with_restart(symmetric_operator, restart_matrix)
    )
    uniform_restart = np.full(
        (len(graph.node_names), 1), 1.0 / len(graph.node_names), dtype=float
    )
    uniform_values, uniform_iterations, uniform_delta = random_walk_with_restart(
        graph, uniform_restart
    )

    gene_indices = np.asarray(
        [graph.node_index[gene] for gene in associated_genes], dtype=int
    )
    association_matrix = build_association_matrix(
        all_drugs, associated_genes, inputs.associations
    )
    all_drug_index = {drug: index for index, drug in enumerate(all_drugs)}
    primary_indices = np.asarray(
        [all_drug_index[drug] for drug in drug_names], dtype=int
    )
    gene_variant_values = _gene_variants(
        column_values[gene_indices, :],
        symmetric_values[gene_indices, :],
        uniform_values[gene_indices, :],
    )

    rows: list[dict[str, Any]] = []
    for variant in VARIANTS:
        c_acc_all = compute_c_acc_matrix(
            association_matrix,
            gene_variant_values[variant],
            pseudo_count=PSEUDO_COUNT,
        )
        c_acc_pct = percentile_columns(c_acc_all[primary_indices, :])
        composite = 0.5 * c_acc_pct + 0.5 * residual_pct[:, None]
        ranks = ordinal_rank_columns(composite, drug_names)
        baseline_ranks = ranks[:, 0]
        loo_ranks = ranks[:, 1:]
        rho = rank_spearman_columns(loo_ranks, baseline_ranks)
        jaccard = top_k_jaccard_columns(loo_ranks, baseline_ranks, TOP_K)
        shifts = np.abs(loo_ranks.astype(int) - baseline_ranks[:, None].astype(int))
        for column, seed in enumerate(seed_names):
            exposed_mask = np.asarray(
                [seed in inputs.associations[drug] for drug in drug_names],
                dtype=bool,
            )
            unexposed_mask = ~exposed_mask
            exposed_shifts = shifts[exposed_mask, column]
            unexposed_shifts = shifts[unexposed_mask, column]
            rows.append(
                {
                    "seed": seed,
                    "variant": variant,
                    "rho": float(rho[column]),
                    "top20_jaccard": float(jaccard[column]),
                    "median_abs_shift": float(np.median(shifts[:, column])),
                    "max_abs_shift": int(np.max(shifts[:, column])),
                    "n_shift_ge_5": int(np.sum(shifts[:, column] >= 5)),
                    "n_shift_ge_10": int(np.sum(shifts[:, column] >= 10)),
                    "exposed_n": int(np.sum(exposed_mask)),
                    "exposed_median_abs_shift": (
                        float(np.median(exposed_shifts))
                        if exposed_shifts.size
                        else None
                    ),
                    "exposed_max_abs_shift": (
                        int(np.max(exposed_shifts)) if exposed_shifts.size else None
                    ),
                    "unexposed_median_abs_shift": float(
                        np.median(unexposed_shifts)
                    ),
                    "unexposed_max_abs_shift": int(np.max(unexposed_shifts)),
                }
            )

    rows.sort(key=lambda row: (row["seed"], row["variant"]))
    variant_summary: dict[str, dict[str, Any]] = {}
    for variant in VARIANTS:
        subset = [row for row in rows if row["variant"] == variant]
        variant_summary[variant] = {
            "minimum_rho": min(float(row["rho"]) for row in subset),
            "median_rho": float(np.median([float(row["rho"]) for row in subset])),
            "minimum_top20_jaccard": min(
                float(row["top20_jaccard"]) for row in subset
            ),
            "largest_abs_shift": max(int(row["max_abs_shift"]) for row in subset),
            "seeds_with_any_shift_ge_10": sum(
                int(row["n_shift_ge_10"]) > 0 for row in subset
            ),
        }

    seed_summary: list[dict[str, Any]] = []
    for seed in seed_names:
        subset = [row for row in rows if row["seed"] == seed]
        exposed_max_values = [
            int(row["exposed_max_abs_shift"])
            for row in subset
            if row["exposed_max_abs_shift"] is not None
        ]
        seed_summary.append(
            {
                "seed": seed,
                "normalized_weight": normalized_seed_weights[seed],
                "exposed_n": int(subset[0]["exposed_n"]),
                "minimum_rho": min(float(row["rho"]) for row in subset),
                "minimum_top20_jaccard": min(
                    float(row["top20_jaccard"]) for row in subset
                ),
                "max_abs_shift": max(int(row["max_abs_shift"]) for row in subset),
                "max_n_shift_ge_10": max(
                    int(row["n_shift_ge_10"]) for row in subset
                ),
                "max_exposed_shift": max(exposed_max_values, default=0),
                "max_unexposed_shift": max(
                    int(row["unexposed_max_abs_shift"]) for row in subset
                ),
            }
        )
    seed_summary.sort(key=lambda row: row["seed"])

    rb1 = next(row for row in seed_summary if row["seed"] == "RB1")
    metrics: dict[str, Any] = {
        "analysis_version": ANALYSIS_VERSION,
        "seed_count": len(seed_names),
        "drug_count": len(drug_names),
        "variant_count": len(VARIANTS),
        "run_count": len(rows),
        "minimum_rho": min(float(row["rho"]) for row in rows),
        "minimum_top20_jaccard": min(
            float(row["top20_jaccard"]) for row in rows
        ),
        "zero_exposure_seed_count": sum(
            int(row["exposed_n"]) == 0 for row in seed_summary
        ),
        "rb1_max_abs_shift": int(rb1["max_abs_shift"]),
        "rb1_worst_shift_rank": _rank_number(
            seed_summary, "RB1", "max_abs_shift"
        ),
        "rb1_exposed_shift_rank": _rank_number(
            seed_summary, "RB1", "max_exposed_shift"
        ),
        "top_influential_seeds": [
            {"seed": row["seed"], "max_abs_shift": row["max_abs_shift"]}
            for row in sorted(
                seed_summary,
                key=lambda row: (-int(row["max_abs_shift"]), row["seed"]),
            )[:8]
        ],
        "convergence": {
            "column_iterations": column_iterations,
            "column_final_max_l1_delta": column_delta,
            "symmetric_iterations": symmetric_iterations,
            "symmetric_final_max_l1_delta": symmetric_delta,
            "uniform_iterations": uniform_iterations,
            "uniform_final_max_l1_delta": uniform_delta,
        },
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "protocol_sha256": _sha256(root / PROTOCOL_PATH),
        "input_sha256": {
            path.as_posix(): _sha256(root / path)
            for path in (
                Path("data/ACC_P0.5C_gene_weights_v1.csv"),
                Path("data/bindex_network/bindex_edges_1304.csv"),
                Path("data/bindex_network/rACC_399_fullSTRING.csv"),
                Path("9606.protein.info.v12.0.txt.gz"),
                Path("9606.protein.links.v12.0.txt.gz"),
            )
        },
    }
    return LeaveOneSeedOutResult(
        rows=tuple(rows),
        variant_summary=variant_summary,
        seed_summary=tuple(seed_summary),
        metrics=metrics,
    )


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("Cannot write an empty CSV")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: (
                        ""
                        if value is None
                        else f"{value:.12g}"
                        if isinstance(value, float)
                        else value
                    )
                    for key, value in row.items()
                }
            )


def _format_optional(value: Any, decimals: int = 1) -> str:
    return "NA" if value is None else f"{float(value):.{decimals}f}"


def _render_report(result: LeaveOneSeedOutResult) -> str:
    m = result.metrics
    top = ", ".join(
        f"{row['seed']} ({row['max_abs_shift']})"
        for row in m["top_influential_seeds"]
    )
    lines = [
        "# Full 45-seed leave-one-out influence audit",
        "",
        f"- Implementation: `{m['analysis_version']}`.",
        f"- Protocol SHA-256: `{m['protocol_sha256']}`.",
        "- Classification: reviewer-requested, result-known post-hoc descriptive analysis.",
        f"- Seeds: {m['seed_count']}; locked drugs: {m['drug_count']}; variants: {m['variant_count']}; deterministic runs: {m['run_count']}.",
        "",
        "## Headline result",
        "",
        f"Across variants, the minimum ADRS-rank Spearman correlation was {m['minimum_rho']:.4f} and the minimum Top-20 Jaccard overlap was {m['minimum_top20_jaccard']:.4f}. RB1 produced a worst-case shift of {m['rb1_max_abs_shift']} ranks and ranked {m['rb1_worst_shift_rank']}/45 by maximum single-drug movement and {m['rb1_exposed_shift_rank']}/45 by maximum directly exposed-drug movement.",
        "",
        f"{m['zero_exposure_seed_count']} of 45 seeds had no direct association to a drug in the locked universe. Largest worst-case shifts: {top}.",
        "",
        "## Variant-level summary",
        "",
        "| Variant | Minimum rho | Median rho | Minimum Top-20 Jaccard | Largest absolute rank shift | Seeds causing any >=10-rank shift |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for variant in VARIANTS:
        row = result.variant_summary[variant]
        lines.append(
            f"| {variant} | {row['minimum_rho']:.4f} | {row['median_rho']:.4f} | {row['minimum_top20_jaccard']:.4f} | {row['largest_abs_shift']} | {row['seeds_with_any_shift_ge_10']} |"
        )
    lines.extend(
        [
            "",
            "## Seed-level worst case across four variants",
            "",
            "| Seed | Normalized weight | Directly exposed drugs | Minimum rho | Minimum Top-20 Jaccard | Maximum absolute shift | Maximum number shifting >=10 | Maximum exposed shift | Maximum unexposed shift |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in result.seed_summary:
        lines.append(
            f"| {row['seed']} | {row['normalized_weight']:.5f} | {row['exposed_n']} | {row['minimum_rho']:.4f} | {row['minimum_top20_jaccard']:.4f} | {row['max_abs_shift']} | {row['max_n_shift_ge_10']} | {row['max_exposed_shift']} | {row['max_unexposed_shift']} |"
        )
    lines.extend(
        [
            "",
            "## Complete seed-by-variant results",
            "",
            "| Seed | Variant | rho | Top-20 Jaccard | Median abs. shift | Max abs. shift | n >=5 | n >=10 | exp n | exp median | exp max | unexp median | unexp max |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in result.rows:
        lines.append(
            f"| {row['seed']} | {row['variant']} | {row['rho']:.4f} | {row['top20_jaccard']:.4f} | {row['median_abs_shift']:.1f} | {row['max_abs_shift']} | {row['n_shift_ge_5']} | {row['n_shift_ge_10']} | {row['exposed_n']} | {_format_optional(row['exposed_median_abs_shift'])} | {_format_optional(row['exposed_max_abs_shift'], 0)} | {row['unexposed_median_abs_shift']:.1f} | {row['unexposed_max_abs_shift']} |"
        )
    convergence = m["convergence"]
    lines.extend(
        [
            "",
            "## Convergence and reproducibility",
            "",
            f"- Python: `{m['environment']['python']}`; NumPy: `{m['environment']['numpy']}`; SciPy: `{m['environment']['scipy']}`.",
            f"- Column-stochastic batch: {convergence['column_iterations']} iterations; final maximum L1 delta `{convergence['column_final_max_l1_delta']:.3e}`.",
            f"- Symmetric batch: {convergence['symmetric_iterations']} iterations; final maximum L1 delta `{convergence['symmetric_final_max_l1_delta']:.3e}`.",
            f"- Uniform-reference run: {convergence['uniform_iterations']} iterations; final maximum L1 delta `{convergence['uniform_final_max_l1_delta']:.3e}`.",
            "",
            "### Input SHA-256",
            "",
        ]
    )
    for path, digest in m["input_sha256"].items():
        lines.append(f"- `{path}`: `{digest}`")
    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "This scan measures deterministic ranking sensitivity when one curated seed is removed. It does not test whether an omitted seed is biologically correct and does not provide external efficacy validation.",
            "",
        ]
    )
    return "\n".join(lines)


def write_outputs(
    project_root: Path,
    result: LeaveOneSeedOutResult,
    output_dir: Path | None = None,
) -> tuple[Path, Path, Path, Path]:
    root = project_root.resolve()
    target = (root / OUTPUT_DIR) if output_dir is None else output_dir.resolve()
    target.mkdir(parents=True, exist_ok=True)
    rows_path = target / "leave_one_seed_out_summary.csv"
    seed_path = target / "leave_one_seed_out_seed_summary.csv"
    metrics_path = target / "leave_one_seed_out_metrics.json"
    report_path = target / "leave_one_seed_out_audit.md"
    _write_csv(rows_path, result.rows)
    _write_csv(seed_path, result.seed_summary)
    metrics_payload = {
        **result.metrics,
        "variant_summary": result.variant_summary,
    }
    metrics_path.write_text(
        json.dumps(metrics_payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    report_path.write_text(_render_report(result), encoding="utf-8")
    return rows_path, seed_path, metrics_path, report_path


def run(project_root: Path, output_dir: Path | None = None) -> LeaveOneSeedOutResult:
    result = run_analysis(project_root)
    write_outputs(project_root, result, output_dir)
    return result


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    result = run(args.project_root, args.output_dir)
    print(json.dumps(result.metrics, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
