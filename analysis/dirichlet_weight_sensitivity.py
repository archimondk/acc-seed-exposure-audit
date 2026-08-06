"""Post-hoc Dirichlet sensitivity of ACC disease-component weights.

The frozen primary model is not modified. This module varies only the five
active disease-biology component coefficients while retaining the 45 seed
genes, the primary network estimator, the locked activity residual and the
locked 108-drug universe.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import stats

from analysis.acc_primary_pipeline import compute_primary_analysis, load_inputs
from analysis.method_strengthening import (
    C_ACC_PSEUDO_COUNT,
    RWR_RESTART,
    THERAPY_ONLY_SEEDS,
    build_association_matrix,
    compute_c_acc_matrix,
    load_string_graph,
    minmax_columns,
    percentile_columns,
    random_walk_with_restart,
)


ANALYSIS_VERSION = "dirichlet-component-weight-sensitivity-v1"
ACTIVE_COMPONENTS = ("G", "R", "P", "L", "S")
COMPONENT_COLUMNS = (
    "Genomic_driver",
    "Recurrence_score",
    "Core_pathway_score",
    "Lineage_or_biomarker",
    "Prognostic_or_subtype",
)
BASE_COMPONENT_WEIGHTS = np.asarray([0.30, 0.20, 0.20, 0.10, 0.05])
DIRICHLET_ALPHA = np.ones(len(ACTIVE_COMPONENTS), dtype=float)
DIRICHLET_DRAWS = 1000
RNG_SEED = 20260729
PRIMARY_DRUG_COUNT = 108
DISEASE_SEED_COUNT = 45
CDK46_DRUGS = ("Abemaciclib", "Palbociclib", "Ribociclib")
TOP_K = 20
PNG_DPI = 1000
TARGET_WIDTH_MM = 170.0

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
GREY = "#A7A9AC"
DARK = "#202124"


@dataclass(frozen=True)
class DirichletSensitivityResult:
    drug_names: tuple[str, ...]
    component_draws: np.ndarray
    c_acc_percentiles: np.ndarray
    adrs_scores: np.ndarray
    ranks: np.ndarray
    draw_rows: tuple[dict[str, Any], ...]
    drug_rows: tuple[dict[str, Any], ...]
    summary: dict[str, Any]


def _read_csv(path: Path, required: Iterable[str]) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Required input is missing: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"Input has no data rows: {path}")
    missing = set(required) - set(rows[0])
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
    return rows


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_seed_component_matrix(path: Path) -> tuple[tuple[str, ...], np.ndarray]:
    rows = _read_csv(path, ("Gene", *COMPONENT_COLUMNS))
    parsed: list[tuple[str, list[float]]] = []
    seen: set[str] = set()
    for row in rows:
        gene = row["Gene"].strip()
        if not gene:
            raise ValueError("Seed gene names must be non-empty")
        if gene in seen:
            raise ValueError(f"Duplicate seed gene: {gene}")
        seen.add(gene)
        if gene in THERAPY_ONLY_SEEDS:
            continue
        values = [float(row[column]) for column in COMPONENT_COLUMNS]
        if not all(math.isfinite(value) and value >= 0 for value in values):
            raise ValueError(f"Invalid component value for {gene}")
        if not any(value > 0 for value in values):
            raise ValueError(f"Retained seed has no positive component: {gene}")
        parsed.append((gene, values))
    parsed.sort(key=lambda item: item[0])
    genes = tuple(item[0] for item in parsed)
    matrix = np.asarray([item[1] for item in parsed], dtype=float)
    if len(genes) != DISEASE_SEED_COUNT:
        raise ValueError(
            f"Expected {DISEASE_SEED_COUNT} disease-only seeds, got {len(genes)}"
        )
    return genes, matrix


def sample_component_weights(
    n_draws: int = DIRICHLET_DRAWS,
    rng_seed: int = RNG_SEED,
) -> np.ndarray:
    if n_draws <= 0:
        raise ValueError("Dirichlet draw count must be positive")
    rng = np.random.default_rng(rng_seed)
    draws = rng.dirichlet(DIRICHLET_ALPHA, size=n_draws)
    if draws.shape != (n_draws, len(ACTIVE_COMPONENTS)):
        raise RuntimeError("Unexpected Dirichlet draw shape")
    if not np.all(np.isfinite(draws)) or not np.all(draws > 0):
        raise ValueError("Dirichlet draws must be finite and strictly positive")
    if not np.allclose(draws.sum(axis=1), 1.0, atol=1e-12):
        raise ValueError("Dirichlet draws do not sum to one")
    return draws


def compute_seed_weight_matrix(
    seed_component_matrix: np.ndarray,
    component_draws: np.ndarray,
) -> np.ndarray:
    components = np.asarray(seed_component_matrix, dtype=float)
    draws = np.asarray(component_draws, dtype=float)
    if components.ndim != 2 or draws.ndim != 2:
        raise ValueError("Component matrices must be two-dimensional")
    if components.shape[1] != len(ACTIVE_COMPONENTS):
        raise ValueError("Seed component matrix has the wrong number of columns")
    if draws.shape[1] != len(ACTIVE_COMPONENTS):
        raise ValueError("Dirichlet matrix has the wrong number of columns")
    if not np.all(np.isfinite(components)) or np.any(components < 0):
        raise ValueError("Seed components must be finite and non-negative")
    if not np.all(np.isfinite(draws)) or np.any(draws < 0):
        raise ValueError("Component draws must be finite and non-negative")
    raw = components @ draws.T
    totals = raw.sum(axis=0)
    if np.any(raw <= 0) or np.any(totals <= 0):
        raise ValueError("Every retained seed must have positive restart weight")
    normalized = raw / totals
    if not np.allclose(normalized.sum(axis=0), 1.0, atol=1e-12):
        raise ValueError("Seed restart weights do not sum to one")
    return normalized


def ordinal_rank_columns(
    scores: np.ndarray,
    drug_names: Sequence[str],
) -> np.ndarray:
    values = np.asarray(scores, dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    if values.ndim != 2 or values.shape[0] != len(drug_names):
        raise ValueError("Score matrix and drug names do not align")
    if not np.all(np.isfinite(values)):
        raise ValueError("Drug scores must be finite")
    names = np.asarray([str(name) for name in drug_names], dtype=object)
    ranks = np.empty(values.shape, dtype=np.int16)
    target = np.arange(1, values.shape[0] + 1, dtype=np.int16)
    for column in range(values.shape[1]):
        order = np.lexsort((names, -values[:, column]))
        ranks[order, column] = target
    expected = np.arange(1, values.shape[0] + 1)
    for column in range(values.shape[1]):
        if not np.array_equal(np.sort(ranks[:, column]), expected):
            raise ValueError(f"Draw {column + 1} is not a complete rank permutation")
    return ranks


def rank_spearman_columns(
    ranks: np.ndarray,
    baseline_ranks: np.ndarray,
) -> np.ndarray:
    matrix = np.asarray(ranks, dtype=float)
    baseline = np.asarray(baseline_ranks, dtype=float)
    if matrix.ndim != 2 or baseline.shape != (matrix.shape[0],):
        raise ValueError("Rank matrices do not align")
    centered_baseline = baseline - baseline.mean()
    centered = matrix - matrix.mean(axis=0, keepdims=True)
    denominator = np.sqrt(
        np.sum(centered_baseline**2) * np.sum(centered**2, axis=0)
    )
    if np.any(denominator <= 0):
        raise ValueError("Rank variance must be positive")
    return (centered_baseline @ centered) / denominator


def top_k_jaccard_columns(
    ranks: np.ndarray,
    baseline_ranks: np.ndarray,
    k: int = TOP_K,
) -> np.ndarray:
    matrix = np.asarray(ranks)
    baseline = np.asarray(baseline_ranks)
    if matrix.ndim != 2 or baseline.shape != (matrix.shape[0],):
        raise ValueError("Rank matrices do not align")
    if not 0 < k < matrix.shape[0]:
        raise ValueError("Top-k must lie between 1 and n_drugs - 1")
    baseline_set = baseline <= k
    draw_sets = matrix <= k
    intersection = np.sum(draw_sets & baseline_set[:, None], axis=0)
    union = np.sum(draw_sets | baseline_set[:, None], axis=0)
    return intersection / union


def _distribution_summary(values: np.ndarray) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    if array.size == 0 or not np.all(np.isfinite(array)):
        raise ValueError("Summary values must be non-empty and finite")
    return {
        "minimum": float(np.min(array)),
        "q05": float(np.quantile(array, 0.05)),
        "median": float(np.median(array)),
        "q95": float(np.quantile(array, 0.95)),
        "maximum": float(np.max(array)),
    }


def _restart_from_seed_matrix(
    n_nodes: int,
    seed_indices: np.ndarray,
    seed_weights: np.ndarray,
) -> np.ndarray:
    weights = np.asarray(seed_weights, dtype=float)
    restart = np.zeros((n_nodes, weights.shape[1]), dtype=float)
    restart[seed_indices, :] = weights
    if not np.all(np.isfinite(restart)):
        raise ValueError("Restart matrix contains non-finite values")
    if not np.allclose(restart.sum(axis=0), 1.0, atol=1e-12):
        raise ValueError("Restart vectors do not sum to one")
    return restart


def run_analysis(
    project_root: Path,
    n_draws: int = DIRICHLET_DRAWS,
    batch_size: int = 64,
    rng_seed: int = RNG_SEED,
) -> DirichletSensitivityResult:
    if batch_size <= 0:
        raise ValueError("Batch size must be positive")
    started = time.perf_counter()
    root = project_root.resolve()
    inputs = load_inputs(root)
    primary = compute_primary_analysis(inputs)
    primary_by_drug = {
        str(row["drug"]): row for row in primary.primary_rows
    }
    drug_names = tuple(sorted(primary_by_drug))
    if len(drug_names) != PRIMARY_DRUG_COUNT:
        raise ValueError(
            f"Expected {PRIMARY_DRUG_COUNT} primary drugs, got {len(drug_names)}"
        )
    baseline_ranks = np.asarray(
        [int(primary_by_drug[drug]["rank_comp"]) for drug in drug_names],
        dtype=np.int16,
    )
    residual_pct = np.asarray(
        [float(primary_by_drug[drug]["residual_pct"]) for drug in drug_names],
        dtype=float,
    )
    locked_c_acc_pct = np.asarray(
        [float(primary_by_drug[drug]["C_ACC_pct"]) for drug in drug_names],
        dtype=float,
    )

    seed_path = root / "data" / "ACC_P0.5C_gene_weights_v1.csv"
    seed_names, seed_components = load_seed_component_matrix(seed_path)
    component_draws = sample_component_weights(n_draws, rng_seed)
    seed_weight_draws = compute_seed_weight_matrix(
        seed_components,
        component_draws,
    )

    all_drugs = tuple(sorted(inputs.associations))
    associated_genes = tuple(sorted(set().union(*inputs.associations.values())))
    graph = load_string_graph(
        root / "9606.protein.info.v12.0.txt.gz",
        root / "9606.protein.links.v12.0.txt.gz",
        required_nodes=set(associated_genes) | set(seed_names),
    )
    seed_indices = np.asarray(
        [graph.node_index[gene] for gene in seed_names],
        dtype=int,
    )
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
        [all_drug_index[drug] for drug in drug_names],
        dtype=int,
    )

    baseline_components = BASE_COMPONENT_WEIGHTS / BASE_COMPONENT_WEIGHTS.sum()
    baseline_seed_weights = compute_seed_weight_matrix(
        seed_components,
        baseline_components[None, :],
    )
    baseline_restart = _restart_from_seed_matrix(
        len(graph.node_names),
        seed_indices,
        baseline_seed_weights,
    )
    baseline_propagated, baseline_iterations, baseline_delta = (
        random_walk_with_restart(
            graph,
            baseline_restart,
            alpha=RWR_RESTART,
        )
    )
    baseline_r_acc = minmax_columns(
        baseline_propagated[gene_indices, :]
    )[:, 0]
    frozen_r_acc = np.asarray(
        [inputs.r_acc[gene] for gene in associated_genes],
        dtype=float,
    )
    baseline_r_acc_max_abs = float(
        np.max(np.abs(baseline_r_acc - frozen_r_acc))
    )
    baseline_r_acc_rho = float(
        stats.spearmanr(baseline_r_acc, frozen_r_acc).statistic
    )
    if (
        baseline_r_acc_max_abs > 1.1e-6
        or baseline_r_acc_rho < 0.999999
    ):
        raise ValueError(
            "Baseline r_ACC failed the frozen reproduction gate: "
            f"max_abs={baseline_r_acc_max_abs}, rho={baseline_r_acc_rho}"
        )
    baseline_c_acc_all = compute_c_acc_matrix(
        association_matrix,
        baseline_r_acc,
        pseudo_count=C_ACC_PSEUDO_COUNT,
    )[:, 0]
    baseline_c_acc_pct = percentile_columns(
        baseline_c_acc_all[primary_indices]
    )[:, 0]
    baseline_c_acc_rho = float(
        stats.spearmanr(baseline_c_acc_pct, locked_c_acc_pct).statistic
    )
    if baseline_c_acc_rho < 0.999999:
        raise ValueError("Baseline C_ACC percentile did not reproduce")

    c_acc_percentiles = np.empty(
        (len(drug_names), n_draws),
        dtype=np.float32,
    )
    adrs_scores = np.empty_like(c_acc_percentiles)
    ranks = np.empty((len(drug_names), n_draws), dtype=np.int16)
    iteration_counts: list[int] = []
    final_deltas: list[float] = []
    for start in range(0, n_draws, batch_size):
        stop = min(start + batch_size, n_draws)
        restart = _restart_from_seed_matrix(
            len(graph.node_names),
            seed_indices,
            seed_weight_draws[:, start:stop],
        )
        propagated, iterations, final_delta = random_walk_with_restart(
            graph,
            restart,
            alpha=RWR_RESTART,
        )
        iteration_counts.extend([iterations] * (stop - start))
        final_deltas.extend([final_delta] * (stop - start))
        r_acc = minmax_columns(propagated[gene_indices, :])
        c_acc_all = compute_c_acc_matrix(
            association_matrix,
            r_acc,
            pseudo_count=C_ACC_PSEUDO_COUNT,
        )
        c_acc_pct = percentile_columns(
            c_acc_all[primary_indices, :]
        )
        adrs = 0.50 * c_acc_pct + 0.50 * residual_pct[:, None]
        c_acc_percentiles[:, start:stop] = c_acc_pct.astype(np.float32)
        adrs_scores[:, start:stop] = adrs.astype(np.float32)
        ranks[:, start:stop] = ordinal_rank_columns(
            adrs,
            drug_names,
        )

    spearman = rank_spearman_columns(ranks, baseline_ranks)
    top20_jaccard = top_k_jaccard_columns(
        ranks,
        baseline_ranks,
        k=TOP_K,
    )
    drug_index = {drug: index for index, drug in enumerate(drug_names)}
    cdk_indices = np.asarray(
        [drug_index[drug] for drug in CDK46_DRUGS],
        dtype=int,
    )
    cdk_mean_rank = ranks[cdk_indices, :].mean(axis=0)

    draw_rows: list[dict[str, Any]] = []
    for draw in range(n_draws):
        row: dict[str, Any] = {
            "draw": draw + 1,
            **{
                f"weight_{component}": float(component_draws[draw, index])
                for index, component in enumerate(ACTIVE_COMPONENTS)
            },
            "spearman_vs_locked_ADRS_rank": float(spearman[draw]),
            "top20_jaccard_vs_locked": float(top20_jaccard[draw]),
            "CDK46_mean_rank": float(cdk_mean_rank[draw]),
        }
        for drug in CDK46_DRUGS:
            row[f"{drug}_rank"] = int(ranks[drug_index[drug], draw])
        draw_rows.append(row)

    top10_probability = np.mean(ranks <= 10, axis=1)
    top20_probability = np.mean(ranks <= 20, axis=1)
    drug_rows: list[dict[str, Any]] = []
    for index, drug in enumerate(drug_names):
        drug_ranks = ranks[index, :].astype(float)
        drug_rows.append(
            {
                "drug": drug,
                "locked_rank": int(baseline_ranks[index]),
                "rank_median": float(np.median(drug_ranks)),
                "rank_q05": float(np.quantile(drug_ranks, 0.05)),
                "rank_q25": float(np.quantile(drug_ranks, 0.25)),
                "rank_q75": float(np.quantile(drug_ranks, 0.75)),
                "rank_q95": float(np.quantile(drug_ranks, 0.95)),
                "rank_sd": float(np.std(drug_ranks, ddof=1)),
                "prob_top10": float(top10_probability[index]),
                "prob_top20": float(top20_probability[index]),
            }
        )
    drug_rows.sort(key=lambda row: (row["locked_rank"], row["drug"]))

    runtime_seconds = time.perf_counter() - started
    summary = {
        "analysis_version": ANALYSIS_VERSION,
        "protocol_id": "dirichlet_component_weight_sensitivity_v1",
        "post_hoc": True,
        "frozen_verdict_revised": False,
        "active_components": list(ACTIVE_COMPONENTS),
        "therapeutic_component_included": False,
        "therapy_only_seeds_included": False,
        "disease_seed_count": len(seed_names),
        "primary_drug_count": len(drug_names),
        "dirichlet_alpha": DIRICHLET_ALPHA.tolist(),
        "draws": n_draws,
        "rng_seed": rng_seed,
        "restart_probability": RWR_RESTART,
        "c_acc_pseudo_count": C_ACC_PSEUDO_COUNT,
        "baseline_component_proportions": {
            component: float(baseline_components[index])
            for index, component in enumerate(ACTIVE_COMPONENTS)
        },
        "baseline_reproduction": {
            "rACC_max_abs_difference": baseline_r_acc_max_abs,
            "rACC_spearman": baseline_r_acc_rho,
            "C_ACC_percentile_spearman": baseline_c_acc_rho,
            "RWR_iterations": baseline_iterations,
            "RWR_final_max_L1_delta": baseline_delta,
        },
        "ADRS_rank_spearman_vs_locked": _distribution_summary(spearman),
        "top20_jaccard_vs_locked": _distribution_summary(top20_jaccard),
        "CDK46_mean_rank": _distribution_summary(cdk_mean_rank),
        "CDK46_drugs": {
            row["drug"]: row
            for row in drug_rows
            if row["drug"] in CDK46_DRUGS
        },
        "n_drugs_prob_top20_ge_0_80": int(
            np.sum(top20_probability >= 0.80)
        ),
        "n_drugs_prob_top20_le_0_20": int(
            np.sum(top20_probability <= 0.20)
        ),
        "quality_control": {
            "component_draw_count": int(component_draws.shape[0]),
            "component_simplex_max_abs_error": float(
                np.max(np.abs(component_draws.sum(axis=1) - 1.0))
            ),
            "restart_simplex_max_abs_error": float(
                np.max(np.abs(seed_weight_draws.sum(axis=0) - 1.0))
            ),
            "all_component_weights_positive": bool(
                np.all(component_draws > 0)
            ),
            "all_seed_weights_positive": bool(
                np.all(seed_weight_draws > 0)
            ),
            "all_rank_columns_complete": bool(
                all(
                    np.array_equal(
                        np.sort(ranks[:, column]),
                        np.arange(1, len(drug_names) + 1),
                    )
                    for column in range(n_draws)
                )
            ),
            "max_RWR_iterations": int(max(iteration_counts)),
            "max_RWR_final_L1_delta": float(max(final_deltas)),
        },
        "runtime_seconds": runtime_seconds,
    }
    return DirichletSensitivityResult(
        drug_names=drug_names,
        component_draws=component_draws,
        c_acc_percentiles=c_acc_percentiles,
        adrs_scores=adrs_scores,
        ranks=ranks,
        draw_rows=tuple(draw_rows),
        drug_rows=tuple(drug_rows),
        summary=summary,
    )


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write an empty CSV: {path}")
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


def _draw_long_rows(
    result: DirichletSensitivityResult,
) -> Iterable[dict[str, Any]]:
    for draw in range(result.ranks.shape[1]):
        for drug_index, drug in enumerate(result.drug_names):
            yield {
                "draw": draw + 1,
                "drug": drug,
                "C_ACC_pct": float(
                    result.c_acc_percentiles[drug_index, draw]
                ),
                "ADRS_comp": float(result.adrs_scores[drug_index, draw]),
                "rank": int(result.ranks[drug_index, draw]),
            }


def _write_long_csv(path: Path, rows: Iterable[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["draw", "drug", "C_ACC_pct", "ADRS_comp", "rank"]
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        count = 0
        for row in rows:
            writer.writerow(
                {
                    key: f"{value:.12g}" if isinstance(value, float) else value
                    for key, value in row.items()
                }
            )
            count += 1
    if count == 0:
        raise ValueError("Cannot write an empty draw-level CSV")


def make_figure(
    result: DirichletSensitivityResult,
    output_prefix: Path,
) -> tuple[Path, Path, Path]:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    draw_spearman = np.asarray(
        [
            float(row["spearman_vs_locked_ADRS_rank"])
            for row in result.draw_rows
        ]
    )
    draw_jaccard = np.asarray(
        [float(row["top20_jaccard_vs_locked"]) for row in result.draw_rows]
    )
    by_drug = {str(row["drug"]): row for row in result.drug_rows}
    top15 = [
        str(row["drug"])
        for row in result.drug_rows
        if int(row["locked_rank"]) <= 15
    ]
    drug_index = {drug: index for index, drug in enumerate(result.drug_names)}

    width_inches = TARGET_WIDTH_MM / 25.4
    fig = plt.figure(figsize=(width_inches, 7.6))
    grid = fig.add_gridspec(3, 1, height_ratios=(1.0, 1.35, 1.0), hspace=0.62)

    ax_a = fig.add_subplot(grid[0, 0])
    ax_a.scatter(
        draw_spearman,
        draw_jaccard,
        s=10,
        alpha=0.24,
        color=BLUE,
        linewidths=0,
        rasterized=True,
    )
    ax_a.scatter([1.0], [1.0], marker="D", s=38, color=ORANGE, zorder=5)
    ax_a.set_xlabel("Spearman correlation with locked ADRS rank")
    ax_a.set_ylabel("Top-20 Jaccard overlap")
    ax_a.set_title("a  Global rank concordance across 1,000 weight draws", loc="left")
    ax_a.grid(True, color="#E6E6E6", linewidth=0.6)
    ax_a.text(
        0.02,
        0.05,
        (
            f"median rho = {np.median(draw_spearman):.3f}\n"
            f"median Jaccard = {np.median(draw_jaccard):.3f}"
        ),
        transform=ax_a.transAxes,
        fontsize=8,
        va="bottom",
        color=DARK,
    )

    ax_b = fig.add_subplot(grid[1, 0])
    x = np.arange(len(top15))
    locked = np.asarray([by_drug[drug]["locked_rank"] for drug in top15])
    median = np.asarray([by_drug[drug]["rank_median"] for drug in top15])
    q05 = np.asarray([by_drug[drug]["rank_q05"] for drug in top15])
    q95 = np.asarray([by_drug[drug]["rank_q95"] for drug in top15])
    ax_b.vlines(x, q05, q95, color=GREY, linewidth=2.0, zorder=1)
    ax_b.scatter(x, median, color=BLUE, s=24, label="Draw median", zorder=3)
    ax_b.scatter(
        x,
        locked,
        color=ORANGE,
        marker="D",
        s=24,
        label="Locked rank",
        zorder=4,
    )
    ax_b.set_xticks(x, top15, rotation=50, ha="right")
    ax_b.set_ylabel("ADRS rank (lower is better)")
    ax_b.set_title(
        "b  Rank intervals for the locked Top 15 (5th–95th percentiles)",
        loc="left",
    )
    ax_b.invert_yaxis()
    ax_b.grid(True, axis="y", color="#E6E6E6", linewidth=0.6)
    ax_b.legend(frameon=False, ncol=2, fontsize=8, loc="upper right")

    ax_c = fig.add_subplot(grid[2, 0])
    cdk_data = [
        result.ranks[drug_index[drug], :].astype(float)
        for drug in CDK46_DRUGS
    ]
    violin = ax_c.violinplot(
        cdk_data,
        showmeans=False,
        showmedians=True,
        showextrema=True,
    )
    for body in violin["bodies"]:
        body.set_facecolor(GREEN)
        body.set_edgecolor(DARK)
        body.set_alpha(0.55)
    for key in ("cmedians", "cmins", "cmaxes", "cbars"):
        violin[key].set_color(DARK)
        violin[key].set_linewidth(0.9)
    locked_cdk = [by_drug[drug]["locked_rank"] for drug in CDK46_DRUGS]
    ax_c.scatter(
        np.arange(1, len(CDK46_DRUGS) + 1),
        locked_cdk,
        marker="D",
        color=ORANGE,
        s=30,
        zorder=4,
    )
    ax_c.set_xticks(np.arange(1, 4), CDK46_DRUGS)
    ax_c.set_ylabel("ADRS rank (lower is better)")
    ax_c.set_title("c  CDK4/6-drug rank distributions", loc="left")
    ax_c.invert_yaxis()
    ax_c.grid(True, axis="y", color="#E6E6E6", linewidth=0.6)

    for axis in (ax_a, ax_b, ax_c):
        axis.tick_params(labelsize=8)
        axis.xaxis.label.set_size(8)
        axis.yaxis.label.set_size(8)
        axis.title.set_size(9)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    fig.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.10)
    png = output_prefix.with_suffix(".png")
    pdf = output_prefix.with_suffix(".pdf")
    svg = output_prefix.with_suffix(".svg")
    fig.savefig(png, dpi=PNG_DPI, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    fig.savefig(svg, facecolor="white")
    plt.close(fig)
    return png, pdf, svg


def write_outputs(
    project_root: Path,
    result: DirichletSensitivityResult,
) -> dict[str, Path]:
    root = project_root.resolve()
    output_dir = root / "results" / "dirichlet_weight_sensitivity"
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "component_weights": output_dir / "component_weight_draws.csv",
        "drug_rank_draws": output_dir / "drug_rank_draws.csv",
        "drug_summary": output_dir / "drug_rank_summary.csv",
        "draw_summary": output_dir / "draw_summary.csv",
        "summary": output_dir / "dirichlet_weight_sensitivity_summary.json",
        "manifest": output_dir / "run_manifest.md",
    }
    component_rows = [
        {
            "draw": draw + 1,
            **{
                f"weight_{component}": float(
                    result.component_draws[draw, index]
                )
                for index, component in enumerate(ACTIVE_COMPONENTS)
            },
        }
        for draw in range(result.component_draws.shape[0])
    ]
    _write_csv(paths["component_weights"], component_rows)
    _write_long_csv(paths["drug_rank_draws"], _draw_long_rows(result))
    _write_csv(paths["drug_summary"], result.drug_rows)
    _write_csv(paths["draw_summary"], result.draw_rows)
    with paths["summary"].open("w", encoding="utf-8") as stream:
        json.dump(result.summary, stream, indent=2, ensure_ascii=False)
        stream.write("\n")

    figure_prefix = (
        root / "figures" / "revision" / "FigS3_dirichlet_weight_sensitivity"
    )
    png, pdf, svg = make_figure(result, figure_prefix)
    paths.update({"figure_png": png, "figure_pdf": pdf, "figure_svg": svg})

    input_paths = (
        root / "experiments" / "dirichlet_component_weight_sensitivity_protocol_v1.md",
        root / "experiments" / "DIRICHLET_WEIGHT_SENSITIVITY_FREEZE.txt",
        root / "data" / "ACC_P0.5C_gene_weights_v1.csv",
        root / "data" / "bindex_network" / "bindex_edges_1304.csv",
        root / "data" / "bindex_network" / "rACC_399_fullSTRING.csv",
        root / "data" / "bindex_network" / "Sactivity_124_v1.csv",
        root / "data" / "bindex_network" / "NCI60_potency_124.csv",
        root / "9606.protein.info.v12.0.txt.gz",
        root / "9606.protein.links.v12.0.txt.gz",
    )
    input_lines = "\n".join(
        f"- `{path.relative_to(root).as_posix()}`: `{_sha256(path)}`"
        for path in input_paths
    )
    output_lines = "\n".join(
        f"- `{path.relative_to(root).as_posix()}`: `{_sha256(path)}`"
        for key, path in paths.items()
        if key != "manifest"
    )
    summary = result.summary
    rho = summary["ADRS_rank_spearman_vs_locked"]
    jaccard = summary["top20_jaccard_vs_locked"]
    manifest = f"""# Dirichlet component-weight sensitivity run manifest

## Frozen design

- Analysis version: `{ANALYSIS_VERSION}`
- Protocol: `dirichlet_component_weight_sensitivity_v1`
- Post-hoc status: `true`
- Active components: `{", ".join(ACTIVE_COMPONENTS)}`
- Fixed disease-only seeds: `{summary["disease_seed_count"]}`
- Draws: `{summary["draws"]}`
- Dirichlet alpha: `{summary["dirichlet_alpha"]}`
- RNG seed: `{summary["rng_seed"]}`
- Primary drug universe: `{summary["primary_drug_count"]}`
- Frozen leakage verdict revised: `false`

## Key descriptive results

- ADRS-rank Spearman versus locked ranking:
  median `{rho["median"]:.6f}`, 5th–95th percentile
  `{rho["q05"]:.6f}`–`{rho["q95"]:.6f}`.
- Top-20 Jaccard versus locked ranking:
  median `{jaccard["median"]:.6f}`, 5th–95th percentile
  `{jaccard["q05"]:.6f}`–`{jaccard["q95"]:.6f}`.
- Drugs with Top-20 probability >=0.80:
  `{summary["n_drugs_prob_top20_ge_0_80"]}`.

## Quality control

- Baseline r_ACC maximum absolute difference:
  `{summary["baseline_reproduction"]["rACC_max_abs_difference"]:.12g}`.
- Baseline r_ACC Spearman:
  `{summary["baseline_reproduction"]["rACC_spearman"]:.12g}`.
- Complete rank permutations:
  `{str(summary["quality_control"]["all_rank_columns_complete"]).lower()}`.
- Maximum RWR iterations:
  `{summary["quality_control"]["max_RWR_iterations"]}`.
- Maximum final RWR L1 delta:
  `{summary["quality_control"]["max_RWR_final_L1_delta"]:.12g}`.

## Runtime environment

- Python: `{platform.python_version()}`
- NumPy: `{np.__version__}`
- SciPy: `{scipy.__version__}`
- Platform: `{platform.platform()}`
- Logical CPUs: `{os.cpu_count()}`
- Wall-clock seconds: `{summary["runtime_seconds"]:.3f}`

## Input SHA-256

{input_lines}

## Output SHA-256

{output_lines}
"""
    paths["manifest"].write_text(manifest, encoding="utf-8")
    return paths


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument("--draws", type=int, default=DIRICHLET_DRAWS)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--rng-seed", type=int, default=RNG_SEED)
    args = parser.parse_args(argv)
    result = run_analysis(
        args.project_root,
        n_draws=args.draws,
        batch_size=args.batch_size,
        rng_seed=args.rng_seed,
    )
    paths = write_outputs(args.project_root, result)
    print(
        json.dumps(
            {
                "status": "ok",
                "draws": result.summary["draws"],
                "runtime_seconds": result.summary["runtime_seconds"],
                "summary": str(paths["summary"]),
                "manifest": str(paths["manifest"]),
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
