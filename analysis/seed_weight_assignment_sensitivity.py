"""Post-hoc sensitivity to ACC seed-weight magnitude and gene assignment.

Protocol Amendment 3 retains the fixed 45-seed membership and locked primary
estimator. W1 replaces the restart weights with 1/45. W2 independently
permutes the normalized baseline weights across the same genes.
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
from analysis.dirichlet_weight_sensitivity import (
    BASE_COMPONENT_WEIGHTS,
    CDK46_DRUGS,
    DISEASE_SEED_COUNT,
    PRIMARY_DRUG_COUNT,
    TOP_K,
    compute_seed_weight_matrix,
    load_seed_component_matrix,
    ordinal_rank_columns,
    rank_spearman_columns,
    top_k_jaccard_columns,
)
from analysis.method_strengthening import (
    C_ACC_PSEUDO_COUNT,
    RWR_RESTART,
    build_association_matrix,
    compute_c_acc_matrix,
    load_string_graph,
    minmax_columns,
    percentile_columns,
    random_walk_with_restart,
)


ANALYSIS_VERSION = "seed-weight-assignment-sensitivity-v1"
PERMUTATION_DRAWS = 1000
RNG_SEED = 20260729
PNG_DPI = 1000
TARGET_WIDTH_MM = 170.0

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PURPLE = "#7A5195"
GREY = "#A7A9AC"
DARK = "#202124"


@dataclass(frozen=True)
class SeedWeightAssignmentResult:
    drug_names: tuple[str, ...]
    baseline_seed_weights: np.ndarray
    uniform_c_acc_percentile: np.ndarray
    uniform_adrs: np.ndarray
    uniform_ranks: np.ndarray
    permutation_c_acc_percentiles: np.ndarray
    permutation_adrs_scores: np.ndarray
    permutation_ranks: np.ndarray
    uniform_rows: tuple[dict[str, Any], ...]
    draw_rows: tuple[dict[str, Any], ...]
    drug_rows: tuple[dict[str, Any], ...]
    dirichlet_spearman: np.ndarray
    dirichlet_jaccard: np.ndarray
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


def _interval_comparison(
    first: Mapping[str, float],
    second: Mapping[str, float],
) -> dict[str, Any]:
    lower = max(float(first["q05"]), float(second["q05"]))
    upper = min(float(first["q95"]), float(second["q95"]))
    union_lower = min(float(first["q05"]), float(second["q05"]))
    union_upper = max(float(first["q95"]), float(second["q95"]))
    overlap = max(0.0, upper - lower)
    union = union_upper - union_lower
    return {
        "q05_q95_intervals_overlap": bool(upper >= lower),
        "overlap_width": overlap,
        "union_width": union,
        "overlap_over_union": float(overlap / union) if union > 0 else 1.0,
        "median_difference_W2_minus_Dirichlet": (
            float(first["median"]) - float(second["median"])
        ),
    }


def make_uniform_seed_weights(n_seeds: int) -> np.ndarray:
    if n_seeds <= 0:
        raise ValueError("Seed count must be positive")
    return np.full((n_seeds, 1), 1.0 / n_seeds, dtype=float)


def make_permuted_seed_weights(
    baseline_weights: np.ndarray,
    n_draws: int = PERMUTATION_DRAWS,
    rng_seed: int = RNG_SEED,
) -> np.ndarray:
    weights = np.asarray(baseline_weights, dtype=float)
    if weights.ndim != 1 or weights.size < 2:
        raise ValueError("Baseline weights must be a one-dimensional vector")
    if not np.all(np.isfinite(weights)) or np.any(weights <= 0):
        raise ValueError("Baseline weights must be finite and strictly positive")
    if not math.isclose(float(weights.sum()), 1.0, abs_tol=1e-12):
        raise ValueError("Baseline weights must sum to one")
    if n_draws <= 0:
        raise ValueError("Permutation draw count must be positive")
    rng = np.random.default_rng(rng_seed)
    draws = np.column_stack(
        [rng.permutation(weights) for _ in range(n_draws)]
    )
    if draws.shape != (weights.size, n_draws):
        raise RuntimeError("Unexpected permutation matrix shape")
    sorted_baseline = np.sort(weights)
    if not all(
        np.array_equal(np.sort(draws[:, column]), sorted_baseline)
        for column in range(n_draws)
    ):
        raise RuntimeError("A W2 draw did not preserve the weight multiset")
    return draws


def _restart_matrix(
    n_nodes: int,
    seed_indices: np.ndarray,
    seed_weights: np.ndarray,
) -> np.ndarray:
    weights = np.asarray(seed_weights, dtype=float)
    if weights.ndim != 2 or weights.shape[0] != seed_indices.size:
        raise ValueError("Seed-weight matrix does not align with seed indices")
    restart = np.zeros((n_nodes, weights.shape[1]), dtype=float)
    restart[seed_indices, :] = weights
    if not np.all(np.isfinite(restart)):
        raise ValueError("Restart matrix contains non-finite values")
    if not np.allclose(restart.sum(axis=0), 1.0, atol=1e-12):
        raise ValueError("Restart vectors do not sum to one")
    return restart


def run_analysis(
    project_root: Path,
    n_draws: int = PERMUTATION_DRAWS,
    batch_size: int = 64,
    rng_seed: int = RNG_SEED,
) -> SeedWeightAssignmentResult:
    if batch_size <= 0:
        raise ValueError("Batch size must be positive")
    started = time.perf_counter()
    root = project_root.resolve()
    inputs = load_inputs(root)
    primary = compute_primary_analysis(inputs)
    primary_by_drug = {str(row["drug"]): row for row in primary.primary_rows}
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
    raw_baseline_weights = seed_components @ BASE_COMPONENT_WEIGHTS
    baseline_seed_weights = compute_seed_weight_matrix(
        seed_components,
        BASE_COMPONENT_WEIGHTS[None, :],
    )[:, 0]
    uniform_weights = make_uniform_seed_weights(len(seed_names))
    permutation_weights = make_permuted_seed_weights(
        baseline_seed_weights,
        n_draws=n_draws,
        rng_seed=rng_seed,
    )
    arm_weights = np.column_stack(
        [uniform_weights[:, 0], permutation_weights]
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

    baseline_restart = _restart_matrix(
        len(graph.node_names),
        seed_indices,
        baseline_seed_weights[:, None],
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
    if baseline_r_acc_max_abs > 1.1e-6 or baseline_r_acc_rho < 0.999999:
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

    n_arms = n_draws + 1
    c_acc_percentiles = np.empty(
        (len(drug_names), n_arms),
        dtype=np.float32,
    )
    adrs_scores = np.empty_like(c_acc_percentiles)
    ranks = np.empty((len(drug_names), n_arms), dtype=np.int16)
    iteration_counts: list[int] = []
    final_deltas: list[float] = []
    for start in range(0, n_arms, batch_size):
        stop = min(start + batch_size, n_arms)
        restart = _restart_matrix(
            len(graph.node_names),
            seed_indices,
            arm_weights[:, start:stop],
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
        c_acc_pct = percentile_columns(c_acc_all[primary_indices, :])
        adrs = 0.50 * c_acc_pct + 0.50 * residual_pct[:, None]
        c_acc_percentiles[:, start:stop] = c_acc_pct.astype(np.float32)
        adrs_scores[:, start:stop] = adrs.astype(np.float32)
        ranks[:, start:stop] = ordinal_rank_columns(adrs, drug_names)

    uniform_c_acc = c_acc_percentiles[:, 0].copy()
    uniform_adrs = adrs_scores[:, 0].copy()
    uniform_ranks = ranks[:, 0].copy()
    permutation_c_acc = c_acc_percentiles[:, 1:].copy()
    permutation_adrs = adrs_scores[:, 1:].copy()
    permutation_ranks = ranks[:, 1:].copy()

    uniform_spearman = float(
        rank_spearman_columns(
            uniform_ranks[:, None],
            baseline_ranks,
        )[0]
    )
    uniform_jaccard = float(
        top_k_jaccard_columns(
            uniform_ranks[:, None],
            baseline_ranks,
            k=TOP_K,
        )[0]
    )
    uniform_abs_change = np.abs(
        uniform_ranks.astype(float) - baseline_ranks.astype(float)
    )

    permutation_spearman_locked = rank_spearman_columns(
        permutation_ranks,
        baseline_ranks,
    )
    permutation_jaccard_locked = top_k_jaccard_columns(
        permutation_ranks,
        baseline_ranks,
        k=TOP_K,
    )
    permutation_spearman_uniform = rank_spearman_columns(
        permutation_ranks,
        uniform_ranks,
    )
    permutation_jaccard_uniform = top_k_jaccard_columns(
        permutation_ranks,
        uniform_ranks,
        k=TOP_K,
    )
    absolute_rank_change = np.abs(
        permutation_ranks.astype(float) - baseline_ranks[:, None]
    )
    mean_absolute_rank_change = absolute_rank_change.mean(axis=0)
    max_absolute_rank_change = absolute_rank_change.max(axis=0)

    dirichlet_rows = _read_csv(
        root
        / "results"
        / "dirichlet_weight_sensitivity"
        / "draw_summary.csv",
        (
            "spearman_vs_locked_ADRS_rank",
            "top20_jaccard_vs_locked",
        ),
    )
    dirichlet_spearman = np.asarray(
        [
            float(row["spearman_vs_locked_ADRS_rank"])
            for row in dirichlet_rows
        ],
        dtype=float,
    )
    dirichlet_jaccard = np.asarray(
        [float(row["top20_jaccard_vs_locked"]) for row in dirichlet_rows],
        dtype=float,
    )

    drug_index = {drug: index for index, drug in enumerate(drug_names)}
    uniform_rows: list[dict[str, Any]] = []
    for index, drug in enumerate(drug_names):
        uniform_rows.append(
            {
                "drug": drug,
                "locked_rank": int(baseline_ranks[index]),
                "uniform_rank": int(uniform_ranks[index]),
                "uniform_rank_change": int(
                    uniform_ranks[index] - baseline_ranks[index]
                ),
                "uniform_C_ACC_pct": float(uniform_c_acc[index]),
                "uniform_ADRS_comp": float(uniform_adrs[index]),
            }
        )
    uniform_rows.sort(key=lambda row: (row["locked_rank"], row["drug"]))

    draw_rows: list[dict[str, Any]] = []
    for draw in range(n_draws):
        draw_rows.append(
            {
                "draw": draw + 1,
                "spearman_vs_locked_ADRS_rank": float(
                    permutation_spearman_locked[draw]
                ),
                "top20_jaccard_vs_locked": float(
                    permutation_jaccard_locked[draw]
                ),
                "spearman_vs_uniform_ADRS_rank": float(
                    permutation_spearman_uniform[draw]
                ),
                "top20_jaccard_vs_uniform": float(
                    permutation_jaccard_uniform[draw]
                ),
                "mean_absolute_rank_change_vs_locked": float(
                    mean_absolute_rank_change[draw]
                ),
                "maximum_absolute_rank_change_vs_locked": float(
                    max_absolute_rank_change[draw]
                ),
                **{
                    f"{drug}_rank": int(
                        permutation_ranks[drug_index[drug], draw]
                    )
                    for drug in CDK46_DRUGS
                },
            }
        )

    top10_probability = np.mean(permutation_ranks <= 10, axis=1)
    top20_probability = np.mean(permutation_ranks <= 20, axis=1)
    drug_rows: list[dict[str, Any]] = []
    for index, drug in enumerate(drug_names):
        drug_ranks = permutation_ranks[index, :].astype(float)
        drug_rows.append(
            {
                "drug": drug,
                "locked_rank": int(baseline_ranks[index]),
                "uniform_rank": int(uniform_ranks[index]),
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

    permutation_spearman_summary = _distribution_summary(
        permutation_spearman_locked
    )
    permutation_jaccard_summary = _distribution_summary(
        permutation_jaccard_locked
    )
    dirichlet_spearman_summary = _distribution_summary(dirichlet_spearman)
    dirichlet_jaccard_summary = _distribution_summary(dirichlet_jaccard)
    runtime_seconds = time.perf_counter() - started
    summary = {
        "analysis_version": ANALYSIS_VERSION,
        "protocol_id": "seed_weight_assignment_sensitivity_v1",
        "post_hoc": True,
        "frozen_verdict_revised": False,
        "disease_seed_count": len(seed_names),
        "primary_drug_count": len(drug_names),
        "permutation_draws": n_draws,
        "rng_seed": rng_seed,
        "restart_probability": RWR_RESTART,
        "c_acc_pseudo_count": C_ACC_PSEUDO_COUNT,
        "baseline_seed_weights": {
            "raw_minimum": float(np.min(raw_baseline_weights)),
            "raw_maximum": float(np.max(raw_baseline_weights)),
            "normalized_minimum": float(np.min(baseline_seed_weights)),
            "normalized_maximum": float(np.max(baseline_seed_weights)),
            "maximum_minimum_ratio": float(
                np.max(baseline_seed_weights)
                / np.min(baseline_seed_weights)
            ),
            "unique_raw_weights": int(np.unique(raw_baseline_weights).size),
            "normalized_weight_cv": float(
                np.std(baseline_seed_weights, ddof=1)
                / np.mean(baseline_seed_weights)
            ),
        },
        "baseline_reproduction": {
            "rACC_max_abs_difference": baseline_r_acc_max_abs,
            "rACC_spearman": baseline_r_acc_rho,
            "C_ACC_percentile_spearman": baseline_c_acc_rho,
            "RWR_iterations": baseline_iterations,
            "RWR_final_max_L1_delta": baseline_delta,
        },
        "W1_uniform": {
            "ADRS_rank_spearman_vs_locked": uniform_spearman,
            "top20_jaccard_vs_locked": uniform_jaccard,
            "mean_absolute_rank_change_vs_locked": float(
                np.mean(uniform_abs_change)
            ),
            "maximum_absolute_rank_change_vs_locked": float(
                np.max(uniform_abs_change)
            ),
            "number_of_changed_ranks": int(
                np.sum(uniform_ranks != baseline_ranks)
            ),
            "number_of_top20_membership_changes": int(
                np.sum(
                    (uniform_ranks <= TOP_K)
                    != (baseline_ranks <= TOP_K)
                )
            ),
            "CDK46_drugs": {
                drug: {
                    "locked_rank": int(baseline_ranks[drug_index[drug]]),
                    "uniform_rank": int(uniform_ranks[drug_index[drug]]),
                }
                for drug in CDK46_DRUGS
            },
        },
        "W2_permuted": {
            "ADRS_rank_spearman_vs_locked": permutation_spearman_summary,
            "top20_jaccard_vs_locked": permutation_jaccard_summary,
            "ADRS_rank_spearman_vs_uniform": _distribution_summary(
                permutation_spearman_uniform
            ),
            "top20_jaccard_vs_uniform": _distribution_summary(
                permutation_jaccard_uniform
            ),
            "mean_absolute_rank_change_vs_locked": _distribution_summary(
                mean_absolute_rank_change
            ),
            "maximum_absolute_rank_change_vs_locked": _distribution_summary(
                max_absolute_rank_change
            ),
            "CDK46_drugs": {
                row["drug"]: row
                for row in drug_rows
                if row["drug"] in CDK46_DRUGS
            },
        },
        "comparison_with_Amendment_2": {
            "Dirichlet_ADRS_rank_spearman_vs_locked": (
                dirichlet_spearman_summary
            ),
            "Dirichlet_top20_jaccard_vs_locked": (
                dirichlet_jaccard_summary
            ),
            "spearman_interval_comparison": _interval_comparison(
                permutation_spearman_summary,
                dirichlet_spearman_summary,
            ),
            "top20_jaccard_interval_comparison": _interval_comparison(
                permutation_jaccard_summary,
                dirichlet_jaccard_summary,
            ),
        },
        "quality_control": {
            "uniform_simplex_max_abs_error": float(
                abs(float(uniform_weights.sum()) - 1.0)
            ),
            "permutation_simplex_max_abs_error": float(
                np.max(np.abs(permutation_weights.sum(axis=0) - 1.0))
            ),
            "all_permutations_preserve_weights": bool(
                all(
                    np.array_equal(
                        np.sort(permutation_weights[:, column]),
                        np.sort(baseline_seed_weights),
                    )
                    for column in range(n_draws)
                )
            ),
            "all_seed_weights_positive": bool(np.all(arm_weights > 0)),
            "all_rank_columns_complete": bool(
                all(
                    np.array_equal(
                        np.sort(ranks[:, column]),
                        np.arange(1, len(drug_names) + 1),
                    )
                    for column in range(n_arms)
                )
            ),
            "max_RWR_iterations": int(max(iteration_counts)),
            "max_RWR_final_L1_delta": float(max(final_deltas)),
        },
        "runtime_seconds": runtime_seconds,
    }
    return SeedWeightAssignmentResult(
        drug_names=drug_names,
        baseline_seed_weights=baseline_seed_weights,
        uniform_c_acc_percentile=uniform_c_acc,
        uniform_adrs=uniform_adrs,
        uniform_ranks=uniform_ranks,
        permutation_c_acc_percentiles=permutation_c_acc,
        permutation_adrs_scores=permutation_adrs,
        permutation_ranks=permutation_ranks,
        uniform_rows=tuple(uniform_rows),
        draw_rows=tuple(draw_rows),
        drug_rows=tuple(drug_rows),
        dirichlet_spearman=dirichlet_spearman,
        dirichlet_jaccard=dirichlet_jaccard,
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


def _permutation_long_rows(
    result: SeedWeightAssignmentResult,
) -> Iterable[dict[str, Any]]:
    for draw in range(result.permutation_ranks.shape[1]):
        for drug_index, drug in enumerate(result.drug_names):
            yield {
                "draw": draw + 1,
                "drug": drug,
                "C_ACC_pct": float(
                    result.permutation_c_acc_percentiles[drug_index, draw]
                ),
                "ADRS_comp": float(
                    result.permutation_adrs_scores[drug_index, draw]
                ),
                "rank": int(result.permutation_ranks[drug_index, draw]),
            }


def _write_long_csv(
    path: Path,
    rows: Iterable[Mapping[str, Any]],
) -> None:
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
    result: SeedWeightAssignmentResult,
    output_prefix: Path,
) -> tuple[Path, Path, Path]:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    summary = result.summary
    baseline_ranks = np.asarray(
        [int(row["locked_rank"]) for row in result.uniform_rows],
        dtype=float,
    )
    uniform_ranks = np.asarray(
        [int(row["uniform_rank"]) for row in result.uniform_rows],
        dtype=float,
    )
    permutation_spearman = np.asarray(
        [
            float(row["spearman_vs_locked_ADRS_rank"])
            for row in result.draw_rows
        ]
    )
    permutation_jaccard = np.asarray(
        [float(row["top20_jaccard_vs_locked"]) for row in result.draw_rows]
    )
    by_drug = {str(row["drug"]): row for row in result.drug_rows}
    top15 = [
        str(row["drug"])
        for row in result.drug_rows
        if int(row["locked_rank"]) <= 15
    ]

    width_inches = TARGET_WIDTH_MM / 25.4
    fig = plt.figure(figsize=(width_inches, 7.0))
    grid = fig.add_gridspec(2, 2, height_ratios=(1.0, 1.18), hspace=0.56, wspace=0.34)

    ax_a = fig.add_subplot(grid[0, 0])
    ax_a.scatter(
        baseline_ranks,
        uniform_ranks,
        s=18,
        alpha=0.70,
        color=BLUE,
        linewidths=0,
    )
    ax_a.plot([1, PRIMARY_DRUG_COUNT], [1, PRIMARY_DRUG_COUNT], color=GREY, lw=1)
    ax_a.set_xlabel("Locked ADRS rank")
    ax_a.set_ylabel("W1 uniform-weight rank")
    ax_a.set_title("a  Membership-only restart", loc="left")
    ax_a.text(
        0.04,
        0.95,
        (
            f"rho = {summary['W1_uniform']['ADRS_rank_spearman_vs_locked']:.4f}\n"
            f"Top-20 Jaccard = {summary['W1_uniform']['top20_jaccard_vs_locked']:.3f}"
        ),
        transform=ax_a.transAxes,
        va="top",
        fontsize=8,
    )

    ax_b = fig.add_subplot(grid[0, 1])
    violin_b = ax_b.violinplot(
        [result.dirichlet_spearman, permutation_spearman],
        showmedians=True,
        showextrema=True,
    )
    for body, color in zip(
        violin_b["bodies"],
        (GREEN, PURPLE),
        strict=True,
    ):
        body.set_facecolor(color)
        body.set_edgecolor(DARK)
        body.set_alpha(0.55)
    for key in ("cmedians", "cmins", "cmaxes", "cbars"):
        violin_b[key].set_color(DARK)
        violin_b[key].set_linewidth(0.9)
    ax_b.set_xticks([1, 2], ["Amendment 2\nDirichlet", "W2\npermuted"])
    ax_b.set_ylabel("Spearman correlation with locked rank")
    ax_b.set_title("b  Global rank concordance", loc="left")

    ax_c = fig.add_subplot(grid[1, 0])
    values = np.unique(
        np.concatenate([result.dirichlet_jaccard, permutation_jaccard])
    )
    width = 0.012
    for offset, data, color, label in (
        (-width / 2, result.dirichlet_jaccard, GREEN, "Amendment 2 Dirichlet"),
        (width / 2, permutation_jaccard, PURPLE, "W2 permuted"),
    ):
        counts = np.asarray([np.mean(data == value) for value in values])
        ax_c.bar(
            values + offset,
            counts,
            width=width,
            color=color,
            alpha=0.65,
            label=label,
        )
    ax_c.set_xlabel("Top-20 Jaccard overlap with locked ranking")
    ax_c.set_ylabel("Proportion of draws")
    ax_c.set_title("c  Top-20 membership concordance", loc="left")
    ax_c.legend(frameon=False, fontsize=7, loc="upper left")

    ax_d = fig.add_subplot(grid[1, 1])
    x = np.arange(len(top15))
    locked = np.asarray([by_drug[drug]["locked_rank"] for drug in top15])
    uniform = np.asarray([by_drug[drug]["uniform_rank"] for drug in top15])
    median = np.asarray([by_drug[drug]["rank_median"] for drug in top15])
    q05 = np.asarray([by_drug[drug]["rank_q05"] for drug in top15])
    q95 = np.asarray([by_drug[drug]["rank_q95"] for drug in top15])
    ax_d.vlines(x, q05, q95, color=GREY, linewidth=1.8, zorder=1)
    ax_d.scatter(x, median, color=PURPLE, s=19, label="W2 median", zorder=3)
    ax_d.scatter(x, uniform, color=BLUE, s=19, label="W1 uniform", zorder=4)
    ax_d.scatter(
        x,
        locked,
        color=ORANGE,
        marker="D",
        s=19,
        label="Locked",
        zorder=5,
    )
    ax_d.set_xticks(x, top15, rotation=55, ha="right")
    ax_d.set_ylabel("ADRS rank (lower is better)")
    ax_d.set_title("d  Locked Top 15 under W1/W2", loc="left")
    ax_d.invert_yaxis()
    ax_d.legend(frameon=False, fontsize=7, ncol=1, loc="lower right")

    for axis in (ax_a, ax_b, ax_c, ax_d):
        axis.grid(True, color="#E6E6E6", linewidth=0.6)
        axis.tick_params(labelsize=7.5)
        axis.xaxis.label.set_size(8)
        axis.yaxis.label.set_size(8)
        axis.title.set_size(9)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    fig.subplots_adjust(left=0.10, right=0.98, top=0.97, bottom=0.13)
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
    result: SeedWeightAssignmentResult,
) -> dict[str, Path]:
    root = project_root.resolve()
    output_dir = root / "results" / "seed_weight_assignment_sensitivity"
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "uniform_ranks": output_dir / "uniform_drug_ranks.csv",
        "permutation_draw_summary": (
            output_dir / "permutation_draw_summary.csv"
        ),
        "permutation_drug_rank_draws": (
            output_dir / "permutation_drug_rank_draws.csv"
        ),
        "permutation_drug_rank_summary": (
            output_dir / "permutation_drug_rank_summary.csv"
        ),
        "summary": output_dir / "seed_weight_assignment_summary.json",
        "manifest": output_dir / "run_manifest.md",
    }
    _write_csv(paths["uniform_ranks"], result.uniform_rows)
    _write_csv(paths["permutation_draw_summary"], result.draw_rows)
    _write_long_csv(
        paths["permutation_drug_rank_draws"],
        _permutation_long_rows(result),
    )
    _write_csv(paths["permutation_drug_rank_summary"], result.drug_rows)
    with paths["summary"].open("w", encoding="utf-8") as stream:
        json.dump(result.summary, stream, indent=2, ensure_ascii=False)
        stream.write("\n")

    figure_prefix = (
        root
        / "figures"
        / "revision"
        / "FigS4_seed_weight_assignment_sensitivity"
    )
    png, pdf, svg = make_figure(result, figure_prefix)
    paths.update({"figure_png": png, "figure_pdf": pdf, "figure_svg": svg})

    input_paths = (
        root
        / "experiments"
        / "seed_weight_assignment_sensitivity_protocol_v1.md",
        root
        / "experiments"
        / "SEED_WEIGHT_ASSIGNMENT_SENSITIVITY_FREEZE.txt",
        root / "data" / "ACC_P0.5C_gene_weights_v1.csv",
        root / "data" / "bindex_network" / "bindex_edges_1304.csv",
        root / "data" / "bindex_network" / "rACC_399_fullSTRING.csv",
        root / "data" / "bindex_network" / "Sactivity_124_v1.csv",
        root / "data" / "bindex_network" / "NCI60_potency_124.csv",
        root
        / "results"
        / "dirichlet_weight_sensitivity"
        / "draw_summary.csv",
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
    w1 = summary["W1_uniform"]
    w2_rho = summary["W2_permuted"]["ADRS_rank_spearman_vs_locked"]
    w2_jaccard = summary["W2_permuted"]["top20_jaccard_vs_locked"]
    manifest = f"""# Seed-weight assignment sensitivity run manifest

## Frozen design

- Analysis version: `{ANALYSIS_VERSION}`
- Protocol: `seed_weight_assignment_sensitivity_v1`
- Post-hoc status: `true`
- Fixed disease-only seeds: `{summary["disease_seed_count"]}`
- W1: uniform restart weight `1/45`
- W2 permutations: `{summary["permutation_draws"]}`
- RNG seed: `{summary["rng_seed"]}`
- Primary drug universe: `{summary["primary_drug_count"]}`
- Frozen leakage verdict revised: `false`

## Key descriptive results

- W1 ADRS-rank Spearman versus locked:
  `{w1["ADRS_rank_spearman_vs_locked"]:.6f}`.
- W1 Top-20 Jaccard versus locked:
  `{w1["top20_jaccard_vs_locked"]:.6f}`.
- W2 ADRS-rank Spearman versus locked:
  median `{w2_rho["median"]:.6f}`, 5th–95th percentile
  `{w2_rho["q05"]:.6f}`–`{w2_rho["q95"]:.6f}`.
- W2 Top-20 Jaccard versus locked:
  median `{w2_jaccard["median"]:.6f}`, 5th–95th percentile
  `{w2_jaccard["q05"]:.6f}`–`{w2_jaccard["q95"]:.6f}`.

## Quality control

- Baseline r_ACC maximum absolute difference:
  `{summary["baseline_reproduction"]["rACC_max_abs_difference"]:.12g}`.
- Baseline r_ACC Spearman:
  `{summary["baseline_reproduction"]["rACC_spearman"]:.12g}`.
- Every W2 draw preserves the baseline weight multiset:
  `{str(summary["quality_control"]["all_permutations_preserve_weights"]).lower()}`.
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
    parser.add_argument(
        "--draws",
        type=int,
        default=PERMUTATION_DRAWS,
    )
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
                "permutation_draws": result.summary["permutation_draws"],
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
