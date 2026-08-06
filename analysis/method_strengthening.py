"""Reviewer-priority baselines, centrality diagnostics and random-seed nulls.

The frozen C1 primary score is never modified here.  This module asks whether
that ranking is reducible to simple baselines and whether the ACC network
component exceeds what is expected from degree-matched random STRING seeds.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
import platform
import sys
from array import array
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import scipy
from scipy import sparse, stats
from scipy.stats import rankdata

from analysis.acc_primary_pipeline import (
    compute_primary_analysis,
    load_inputs,
)
from analysis.mechanism_enrichment import benjamini_hochberg


ANALYSIS_VERSION = "method-strengthening-v2-reviewer-null-resolution"
STRING_THRESHOLD = 400
RWR_RESTART = 0.4
PAGERANK_DAMPING = 0.85
C_ACC_PSEUDO_COUNT = 3.0
NULL_DRAWS = 10_000
RNG_SEED = 20260727
REQUESTED_DEGREE_BINS = 10
CDK46_DRUGS = ("Abemaciclib", "Palbociclib", "Ribociclib")
THERAPY_ONLY_SEEDS = frozenset({"MGMT", "SLFN11", "ABCB1", "SOAT1", "UBA1"})


@dataclass(frozen=True)
class StringGraph:
    node_names: tuple[str, ...]
    node_index: Mapping[str, int]
    adjacency: sparse.csr_matrix
    transition: sparse.csr_matrix
    degree: np.ndarray
    strength: np.ndarray
    dangling: np.ndarray
    source_edge_rows: int


@dataclass(frozen=True)
class StrengtheningResult:
    baseline_rows: tuple[dict[str, Any], ...]
    centrality_gene_rows: tuple[dict[str, Any], ...]
    centrality_drug_rows: tuple[dict[str, Any], ...]
    matched_seed_rows: tuple[dict[str, Any], ...]
    null_rows: tuple[dict[str, Any], ...]
    metrics: dict[str, Any]


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


def load_disease_seed_weights(path: Path) -> dict[str, float]:
    rows = _read_csv(
        path,
        (
            "Gene",
            "Genomic_driver",
            "Recurrence_score",
            "Core_pathway_score",
            "Lineage_or_biomarker",
            "Prognostic_or_subtype",
        ),
    )
    weights: dict[str, float] = {}
    for row in rows:
        gene = row["Gene"].strip()
        if gene in THERAPY_ONLY_SEEDS:
            continue
        weight = (
            0.30 * float(row["Genomic_driver"])
            + 0.20 * float(row["Recurrence_score"])
            + 0.20 * float(row["Core_pathway_score"])
            + 0.10 * float(row["Lineage_or_biomarker"])
            + 0.05 * float(row["Prognostic_or_subtype"])
        )
        if weight > 0:
            weights[gene] = weight
    if not weights:
        raise ValueError("Disease seed set is empty")
    return weights


def load_string_graph(
    protein_info_path: Path,
    protein_links_path: Path,
    required_nodes: Iterable[str] = (),
    threshold: int = STRING_THRESHOLD,
) -> StringGraph:
    if not protein_info_path.is_file() or not protein_links_path.is_file():
        raise FileNotFoundError("STRING protein-info or protein-links file is missing")

    protein_to_symbol: dict[str, str] = {}
    with gzip.open(protein_info_path, "rt", encoding="utf-8") as stream:
        header = next(stream).rstrip("\n").split("\t")
        try:
            protein_index = header.index("#string_protein_id")
            symbol_index = header.index("preferred_name")
        except ValueError as error:
            raise ValueError("Unexpected STRING protein-info header") from error
        for line in stream:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(protein_index, symbol_index):
                continue
            protein_to_symbol[parts[protein_index]] = parts[symbol_index]

    node_names = tuple(
        sorted(set(protein_to_symbol.values()).union(str(x) for x in required_nodes))
    )
    node_index = {name: index for index, name in enumerate(node_names)}
    rows = array("i")
    columns = array("i")
    values = array("d")
    source_edge_rows = 0
    with gzip.open(protein_links_path, "rt", encoding="utf-8") as stream:
        next(stream)
        for line in stream:
            protein_a, protein_b, score_text = line.split()
            score = int(score_text)
            if score < threshold:
                continue
            gene_a = protein_to_symbol.get(protein_a)
            gene_b = protein_to_symbol.get(protein_b)
            if not gene_a or not gene_b or gene_a == gene_b:
                continue
            index_a = node_index[gene_a]
            index_b = node_index[gene_b]
            weight = score / 1000.0
            rows.extend((index_a, index_b))
            columns.extend((index_b, index_a))
            values.extend((weight, weight))
            source_edge_rows += 1

    shape = (len(node_names), len(node_names))
    adjacency = sparse.coo_matrix(
        (
            np.frombuffer(values, dtype=np.float64),
            (
                np.frombuffer(rows, dtype=np.int32),
                np.frombuffer(columns, dtype=np.int32),
            ),
        ),
        shape=shape,
        dtype=np.float64,
    ).tocsr()
    adjacency.sum_duplicates()
    adjacency.eliminate_zeros()
    degree = np.diff(adjacency.indptr).astype(float)
    strength = np.asarray(adjacency.sum(axis=0)).ravel()
    inverse_strength = np.zeros_like(strength)
    nonzero = strength > 0
    inverse_strength[nonzero] = 1.0 / strength[nonzero]
    transition = (adjacency @ sparse.diags(inverse_strength)).tocsr()
    return StringGraph(
        node_names=node_names,
        node_index=node_index,
        adjacency=adjacency,
        transition=transition,
        degree=degree,
        strength=strength,
        dangling=~nonzero,
        source_edge_rows=source_edge_rows,
    )


def build_restart_matrix(
    graph: StringGraph,
    columns: Sequence[Mapping[str, float]],
) -> np.ndarray:
    restart = np.zeros((len(graph.node_names), len(columns)), dtype=float)
    for column_index, weights in enumerate(columns):
        for gene, weight in weights.items():
            if gene not in graph.node_index:
                raise ValueError(f"Seed is absent from STRING node index: {gene}")
            restart[graph.node_index[gene], column_index] = float(weight)
    totals = restart.sum(axis=0)
    if np.any(totals <= 0):
        raise ValueError("Every restart column must contain positive mass")
    restart /= totals
    return restart


def random_walk_with_restart(
    graph: StringGraph,
    restart: np.ndarray,
    alpha: float = RWR_RESTART,
    tolerance: float = 1e-10,
    max_iterations: int = 500,
) -> tuple[np.ndarray, int, float]:
    if restart.ndim == 1:
        restart = restart[:, None]
    if restart.shape[0] != len(graph.node_names):
        raise ValueError("Restart matrix does not match STRING graph size")
    if not 0 < alpha <= 1:
        raise ValueError("RWR alpha must lie in (0, 1]")

    propagated = restart.copy()
    final_delta = math.inf
    for iteration in range(1, max_iterations + 1):
        dangling_mass = propagated[graph.dangling, :].sum(axis=0)
        updated = (1.0 - alpha) * (
            graph.transition @ propagated + restart * dangling_mass
        ) + alpha * restart
        final_delta = float(
            np.max(np.sum(np.abs(updated - propagated), axis=0))
        )
        propagated = updated
        if final_delta < tolerance:
            break
    else:
        raise RuntimeError(
            f"RWR failed to converge in {max_iterations} iterations; "
            f"last max L1 delta={final_delta}"
        )
    if not np.allclose(propagated.sum(axis=0), 1.0, atol=1e-9):
        raise ValueError("RWR probability mass is not conserved")
    return propagated, iteration, final_delta


def weighted_pagerank(
    graph: StringGraph,
    damping: float = PAGERANK_DAMPING,
    tolerance: float = 1e-12,
    max_iterations: int = 1000,
) -> tuple[np.ndarray, int, float]:
    if not 0 < damping < 1:
        raise ValueError("PageRank damping must lie in (0, 1)")
    n_nodes = len(graph.node_names)
    uniform = np.full(n_nodes, 1.0 / n_nodes, dtype=float)
    pagerank = uniform.copy()
    final_delta = math.inf
    for iteration in range(1, max_iterations + 1):
        dangling_mass = float(pagerank[graph.dangling].sum())
        updated = damping * (
            graph.transition @ pagerank + uniform * dangling_mass
        ) + (1.0 - damping) * uniform
        final_delta = float(np.abs(updated - pagerank).sum())
        pagerank = updated
        if final_delta < tolerance:
            break
    else:
        raise RuntimeError("PageRank failed to converge")
    return pagerank, iteration, final_delta


def minmax_columns(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    minimum = values.min(axis=0)
    span = values.max(axis=0) - minimum
    if np.any(span <= 0):
        raise ValueError("Min-max scaling requires non-constant columns")
    return (values - minimum) / span


def percentile_columns(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    if values.shape[0] < 2:
        raise ValueError("At least two rows are required for percentile ranks")
    return (rankdata(values, axis=0, method="average") - 1.0) / (
        values.shape[0] - 1.0
    )


def build_association_matrix(
    drug_names: Sequence[str],
    gene_names: Sequence[str],
    associations: Mapping[str, Iterable[str]],
) -> sparse.csr_matrix:
    gene_index = {gene: index for index, gene in enumerate(gene_names)}
    rows: list[int] = []
    columns: list[int] = []
    for drug_index, drug in enumerate(drug_names):
        for gene in associations[drug]:
            if gene not in gene_index:
                raise ValueError(f"Association gene is outside the gene matrix: {gene}")
            rows.append(drug_index)
            columns.append(gene_index[gene])
    values = np.ones(len(rows), dtype=float)
    return sparse.csr_matrix(
        (values, (rows, columns)),
        shape=(len(drug_names), len(gene_names)),
    )


def compute_c_acc_matrix(
    association_matrix: sparse.csr_matrix,
    propagated_gene_values: np.ndarray,
    pseudo_count: float = C_ACC_PSEUDO_COUNT,
) -> np.ndarray:
    propagated = np.asarray(propagated_gene_values, dtype=float)
    if propagated.ndim == 1:
        propagated = propagated[:, None]
    if association_matrix.shape[1] != propagated.shape[0]:
        raise ValueError("Association and propagated-gene matrices do not align")
    n_associations = np.asarray(association_matrix.sum(axis=1)).ravel()
    if np.any(n_associations <= 0):
        raise ValueError("Every drug must have at least one associated gene")
    gene_edge_counts = np.asarray(association_matrix.sum(axis=0)).ravel()
    background = (gene_edge_counts @ propagated) / gene_edge_counts.sum()
    target_sums = association_matrix @ propagated
    return (target_sums + pseudo_count * background[None, :]) / (
        n_associations[:, None] + pseudo_count
    )


def _degree_bins(
    degrees: np.ndarray,
    n_bins: int,
) -> tuple[np.ndarray, np.ndarray]:
    transformed = np.log2(np.asarray(degrees, dtype=float) + 1.0)
    edges = np.unique(
        np.quantile(transformed, np.linspace(0.0, 1.0, n_bins + 1))
    )
    if len(edges) < 2:
        raise ValueError("Degree distribution cannot be binned")
    labels = np.digitize(transformed, edges[1:-1], right=True)
    return labels, edges


def generate_degree_matched_seed_sets(
    seed_weights: Mapping[str, float],
    node_names: Sequence[str],
    degrees: np.ndarray,
    n_draws: int,
    rng_seed: int,
    n_bins: int = REQUESTED_DEGREE_BINS,
) -> tuple[list[dict[str, Any]], np.ndarray]:
    if n_draws < 1:
        raise ValueError("n_draws must be positive")
    if len(node_names) != len(degrees):
        raise ValueError("Node names and degree vector do not align")
    node_index = {name: index for index, name in enumerate(node_names)}
    missing = set(seed_weights) - set(node_index)
    if missing:
        raise ValueError(f"Seeds missing from degree universe: {sorted(missing)}")
    true_seed_set = set(seed_weights)

    selected_labels: np.ndarray | None = None
    selected_edges: np.ndarray | None = None
    selected_pools: dict[int, np.ndarray] | None = None
    for candidate_bin_count in range(n_bins, 0, -1):
        labels, edges = _degree_bins(degrees, candidate_bin_count)
        pools: dict[int, np.ndarray] = {}
        feasible = True
        seed_bin_ids = {
            int(labels[node_index[seed]]) for seed in seed_weights
        }
        for bin_id in sorted(seed_bin_ids):
            seed_count = sum(
                labels[node_index[seed]] == bin_id for seed in seed_weights
            )
            candidates = np.asarray(
                [
                    index
                    for index, name in enumerate(node_names)
                    if labels[index] == bin_id and name not in true_seed_set
                ],
                dtype=int,
            )
            if len(candidates) < seed_count:
                feasible = False
                break
            pools[int(bin_id)] = candidates
        if feasible:
            selected_labels = labels
            selected_edges = edges
            selected_pools = pools
            break
    if selected_labels is None or selected_edges is None or selected_pools is None:
        raise ValueError("No feasible degree-bin matching scheme")

    seeds_by_bin: dict[int, list[str]] = {}
    for seed in sorted(seed_weights):
        bin_id = int(selected_labels[node_index[seed]])
        seeds_by_bin.setdefault(bin_id, []).append(seed)

    rng = np.random.default_rng(rng_seed)
    rows: list[dict[str, Any]] = []
    for replicate in range(n_draws):
        used: set[str] = set()
        for bin_id, seeds in sorted(seeds_by_bin.items()):
            chosen_indices = rng.choice(
                selected_pools[bin_id],
                size=len(seeds),
                replace=False,
            )
            for matched_seed, chosen_index in zip(
                seeds, chosen_indices, strict=True
            ):
                null_seed = node_names[int(chosen_index)]
                if null_seed in used:
                    raise AssertionError("Null seed was sampled twice in one draw")
                used.add(null_seed)
                rows.append(
                    {
                        "replicate": replicate,
                        "matched_seed": matched_seed,
                        "null_seed": null_seed,
                        "weight": float(seed_weights[matched_seed]),
                        "matched_seed_degree": int(
                            degrees[node_index[matched_seed]]
                        ),
                        "null_seed_degree": int(degrees[int(chosen_index)]),
                        "matched_seed_degree_bin": bin_id,
                        "null_seed_degree_bin": int(
                            selected_labels[int(chosen_index)]
                        ),
                    }
                )
    return rows, selected_edges


def empirical_upper_p(observed: float, null_values: np.ndarray) -> float:
    null = np.asarray(null_values, dtype=float)
    if null.ndim != 1 or len(null) < 1 or not np.isfinite(null).all():
        raise ValueError("Null values must be a non-empty finite vector")
    exceedances = int(np.count_nonzero(null >= observed))
    return (exceedances + 1.0) / (len(null) + 1.0)


def exact_lower_tail_rank_sum(
    all_ranks: np.ndarray,
    class_indices: Sequence[int],
) -> float:
    ranks = np.asarray(all_ranks, dtype=float)
    if ranks.ndim != 1:
        raise ValueError("Ranks must be one-dimensional")
    class_indices = tuple(class_indices)
    if not class_indices:
        raise ValueError("Class indices cannot be empty")
    observed = float(ranks[list(class_indices)].sum())
    favorable = 0
    total = 0
    for combination in itertools.combinations(ranks, len(class_indices)):
        total += 1
        favorable += sum(combination) <= observed + 1e-12
    return favorable / total


def partial_spearman(
    x: np.ndarray,
    y: np.ndarray,
    control: np.ndarray,
) -> float:
    x_rank = rankdata(np.asarray(x, dtype=float), method="average")
    y_rank = rankdata(np.asarray(y, dtype=float), method="average")
    z_rank = rankdata(np.asarray(control, dtype=float), method="average")
    if not (len(x_rank) == len(y_rank) == len(z_rank)):
        raise ValueError("Partial-correlation vectors must have equal length")
    design = np.column_stack((np.ones(len(z_rank)), z_rank))
    x_residual = x_rank - design @ np.linalg.lstsq(
        design, x_rank, rcond=None
    )[0]
    y_residual = y_rank - design @ np.linalg.lstsq(
        design, y_rank, rcond=None
    )[0]
    if np.std(x_residual) == 0 or np.std(y_residual) == 0:
        raise ValueError("Partial Spearman residual is constant")
    return float(np.corrcoef(x_residual, y_residual)[0, 1])


def _bootstrap_spearman(
    x: np.ndarray,
    y: np.ndarray,
    n_resamples: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    values: list[float] = []
    n = len(x)
    for _ in range(n_resamples):
        sample = rng.integers(0, n, size=n)
        rho = stats.spearmanr(x[sample], y[sample]).statistic
        if np.isfinite(rho):
            values.append(float(rho))
    if len(values) < n_resamples * 0.99:
        raise ValueError("Too many undefined bootstrap correlations")
    return tuple(float(value) for value in np.percentile(values, [2.5, 97.5]))


def _correlation_record(
    x: np.ndarray,
    y: np.ndarray,
    rng: np.random.Generator,
    n_bootstrap: int = 2000,
) -> dict[str, Any]:
    result = stats.spearmanr(x, y)
    ci_low, ci_high = _bootstrap_spearman(x, y, n_bootstrap, rng)
    return {
        "n": len(x),
        "rho": float(result.statistic),
        "p_two_sided_descriptive": float(result.pvalue),
        "bootstrap_ci_95": [ci_low, ci_high],
        "bootstrap_resamples": n_bootstrap,
    }


def _rank_scores(
    drug_names: Sequence[str],
    scores: Mapping[str, float],
) -> tuple[np.ndarray, list[str]]:
    values = np.asarray([scores[drug] for drug in drug_names], dtype=float)
    ranks = rankdata(-values, method="average")
    deterministic_order = sorted(drug_names, key=lambda drug: (-scores[drug], drug))
    return ranks, deterministic_order


def build_baseline_comparison(
    primary_rows: Sequence[Mapping[str, Any]],
    external_scores: Mapping[str, float],
    associations: Mapping[str, Iterable[str]],
    seed_weights: Mapping[str, float],
    random_null_mean: Mapping[str, float],
    cdk_null_p: float,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    drug_names = tuple(str(row["drug"]) for row in primary_rows)
    row_by_drug = {str(row["drug"]): row for row in primary_rows}
    seed_set = set(seed_weights)
    score_maps: list[tuple[str, str, dict[str, float]]] = [
        (
            "ADRS_comp",
            "locked_primary",
            {drug: float(row_by_drug[drug]["ADRS_comp"]) for drug in drug_names},
        ),
        (
            "raw_MIPE_potency",
            "simple_baseline",
            {
                drug: float(row_by_drug[drug]["MIPE_potency_pct"])
                for drug in drug_names
            },
        ),
        (
            "residual_alone",
            "model_component",
            {
                drug: float(row_by_drug[drug]["residual_pct"])
                for drug in drug_names
            },
        ),
        (
            "C_ACC_alone",
            "model_component",
            {
                drug: float(row_by_drug[drug]["C_ACC_pct"])
                for drug in drug_names
            },
        ),
        (
            "association_count",
            "simple_baseline",
            {drug: float(len(tuple(associations[drug]))) for drug in drug_names},
        ),
        (
            "direct_seed_overlap_fraction",
            "simple_baseline",
            {
                drug: len(set(associations[drug]) & seed_set)
                / len(tuple(associations[drug]))
                for drug in drug_names
            },
        ),
        (
            "S_external",
            "exploratory_non_independent",
            {drug: float(external_scores[drug]) for drug in drug_names},
        ),
        (
            "degree_matched_random_seed_mean",
            "random_network_baseline",
            {drug: float(random_null_mean[drug]) for drug in drug_names},
        ),
    ]
    primary_scores = score_maps[0][2]
    locked_primary_ranks = np.asarray(
        [int(row_by_drug[drug]["rank_comp"]) for drug in drug_names],
        dtype=float,
    )
    primary_top20 = {
        drug
        for drug, rank in zip(drug_names, locked_primary_ranks, strict=True)
        if rank <= 20
    }
    cdk_indices = tuple(drug_names.index(drug) for drug in CDK46_DRUGS)
    rows: list[dict[str, Any]] = []
    for name, role, scores in score_maps:
        if name == "ADRS_comp":
            ranks = locked_primary_ranks
            order = [
                drug
                for drug, _ in sorted(
                    zip(drug_names, ranks, strict=True),
                    key=lambda item: item[1],
                )
            ]
        else:
            ranks, order = _rank_scores(drug_names, scores)
        top20 = set(order[:20])
        overlap = len(primary_top20 & top20)
        p_value = exact_lower_tail_rank_sum(ranks, cdk_indices)
        cdk_ranks = [float(ranks[index]) for index in cdk_indices]
        rows.append(
            {
                "ranking": name,
                "role": role,
                "n_universe": len(drug_names),
                "spearman_vs_ADRS_comp": float(
                    stats.spearmanr(
                        [primary_scores[drug] for drug in drug_names],
                        [scores[drug] for drug in drug_names],
                    ).statistic
                ),
                "top20_overlap_n": overlap,
                "top20_jaccard": overlap / len(primary_top20 | top20),
                "cdk46_ranks": "; ".join(f"{rank:g}" for rank in cdk_ranks),
                "cdk46_mean_rank": float(np.mean(cdk_ranks)),
                "cdk46_p_exact": p_value,
            }
        )
    q_values = benjamini_hochberg([row["cdk46_p_exact"] for row in rows])
    for row, q_value in zip(rows, q_values, strict=True):
        row["cdk46_q_bh_across_rankings"] = q_value

    substitutes = [
        row
        for row in rows[1:]
        if abs(row["spearman_vs_ADRS_comp"]) >= 0.90
        and row["top20_jaccard"] >= 0.80
    ]
    adrs_p = rows[0]["cdk46_p_exact"]
    h3_better_than_all = all(
        adrs_p < row["cdk46_p_exact"] for row in rows[1:]
    )
    hypothesis = {
        "H1_nonredundancy": {
            "status": "retired_descriptive_only",
            "substituting_baselines": [row["ranking"] for row in substitutes],
            "rule": (
                "retired as a hypothesis test because ADRS_comp is an equal-weight "
                "combination of the two ranked components; baseline similarity is "
                "reported only as a composition audit"
            ),
        },
        "composition_audit": {
            "status": "descriptive",
            "substituting_baselines": [row["ranking"] for row in substitutes],
            "rule": (
                "report Spearman correlation and Top-20 overlap without an "
                "inferential non-redundancy claim"
            ),
        },
        "H3_CDK46_robustness": {
            "status": (
                "supported"
                if h3_better_than_all and cdk_null_p < 0.05
                else "not_supported"
            ),
            "ADRS_p_lower_than_all_baselines": h3_better_than_all,
            "degree_matched_seed_group_p": cdk_null_p,
            "rule": "ADRS exact P lower than every baseline and null P<0.05",
        },
    }
    return rows, hypothesis


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


def _rebuild_observed_network_scores(
    graph: StringGraph,
    seed_weights: Mapping[str, float],
    associated_genes: Sequence[str],
) -> tuple[np.ndarray, int, float]:
    restart = build_restart_matrix(graph, (seed_weights,))
    propagated, iterations, delta = random_walk_with_restart(graph, restart)
    indices = [graph.node_index[gene] for gene in associated_genes]
    return minmax_columns(propagated[indices, :])[:, 0], iterations, delta


def run_analysis(
    project_root: Path,
    n_null: int = NULL_DRAWS,
    null_batch_size: int = 64,
) -> StrengtheningResult:
    root = project_root.resolve()
    inputs = load_inputs(root)
    primary_result = compute_primary_analysis(inputs)
    primary_rows = list(primary_result.primary_rows)
    primary_drugs = tuple(str(row["drug"]) for row in primary_rows)
    all_drugs = tuple(sorted(inputs.associations))
    associated_genes = tuple(sorted(set().union(*inputs.associations.values())))

    seed_path = root / "data" / "ACC_P0.5C_gene_weights_v1.csv"
    protein_info_path = root / "9606.protein.info.v12.0.txt.gz"
    protein_links_path = root / "9606.protein.links.v12.0.txt.gz"
    frozen_r_acc_path = root / "data" / "bindex_network" / "rACC_399_fullSTRING.csv"
    seed_weights = load_disease_seed_weights(seed_path)
    graph = load_string_graph(
        protein_info_path,
        protein_links_path,
        required_nodes=set(associated_genes) | set(seed_weights),
    )

    observed_r_acc, observed_iterations, observed_delta = (
        _rebuild_observed_network_scores(
            graph,
            seed_weights,
            associated_genes,
        )
    )
    frozen_rows = _read_csv(frozen_r_acc_path, ("gene", "rACC_full"))
    frozen_r_acc_map = {
        row["gene"].strip(): float(row["rACC_full"]) for row in frozen_rows
    }
    frozen_r_acc = np.asarray(
        [frozen_r_acc_map[gene] for gene in associated_genes],
        dtype=float,
    )
    r_acc_max_abs_difference = float(
        np.max(np.abs(observed_r_acc - frozen_r_acc))
    )
    r_acc_spearman = float(
        stats.spearmanr(observed_r_acc, frozen_r_acc).statistic
    )
    if r_acc_max_abs_difference > 1.1e-6 or r_acc_spearman < 0.999999:
        raise ValueError(
            "Recomputed full-STRING r_ACC failed the frozen-result gate: "
            f"max_abs={r_acc_max_abs_difference}, rho={r_acc_spearman}"
        )

    association_matrix = build_association_matrix(
        all_drugs,
        associated_genes,
        inputs.associations,
    )
    observed_c_acc_all = compute_c_acc_matrix(
        association_matrix,
        observed_r_acc,
    )[:, 0]
    all_drug_index = {drug: index for index, drug in enumerate(all_drugs)}
    primary_indices = np.asarray(
        [all_drug_index[drug] for drug in primary_drugs],
        dtype=int,
    )
    observed_c_acc_pct = percentile_columns(
        observed_c_acc_all[primary_indices]
    )[:, 0]

    matched_seed_rows, degree_bin_edges = generate_degree_matched_seed_sets(
        seed_weights,
        graph.node_names,
        graph.degree,
        n_draws=n_null,
        rng_seed=RNG_SEED,
        n_bins=REQUESTED_DEGREE_BINS,
    )
    rows_by_replicate: dict[int, list[dict[str, Any]]] = {}
    for row in matched_seed_rows:
        rows_by_replicate.setdefault(int(row["replicate"]), []).append(row)
    null_c_acc_pct = np.empty((len(primary_drugs), n_null), dtype=float)
    null_iteration_counts: list[int] = []
    null_final_deltas: list[float] = []
    for batch_start in range(0, n_null, null_batch_size):
        batch_replicates = list(
            range(batch_start, min(batch_start + null_batch_size, n_null))
        )
        batch_weights = []
        for replicate in batch_replicates:
            batch_weights.append(
                {
                    row["null_seed"]: float(row["weight"])
                    for row in rows_by_replicate[replicate]
                }
            )
        restart = build_restart_matrix(graph, batch_weights)
        propagated, iterations, delta = random_walk_with_restart(graph, restart)
        null_iteration_counts.extend([iterations] * len(batch_replicates))
        null_final_deltas.extend([delta] * len(batch_replicates))
        gene_indices = [graph.node_index[gene] for gene in associated_genes]
        null_r_acc = minmax_columns(propagated[gene_indices, :])
        null_c_acc_all = compute_c_acc_matrix(
            association_matrix,
            null_r_acc,
        )
        null_c_acc_pct[:, batch_replicates] = percentile_columns(
            null_c_acc_all[primary_indices, :]
        )

    null_mean = null_c_acc_pct.mean(axis=1)
    null_std = null_c_acc_pct.std(axis=1, ddof=1)
    empirical_p = np.asarray(
        [
            empirical_upper_p(observed_c_acc_pct[index], null_c_acc_pct[index])
            for index in range(len(primary_drugs))
        ],
        dtype=float,
    )
    q_values = np.asarray(benjamini_hochberg(empirical_p.tolist()), dtype=float)
    z_values = np.divide(
        observed_c_acc_pct - null_mean,
        null_std,
        out=np.full_like(null_std, np.nan),
        where=null_std > 0,
    )
    null_rows = [
        {
            "drug": drug,
            "observed_C_ACC_pct": float(observed_c_acc_pct[index]),
            "null_mean_C_ACC_pct": float(null_mean[index]),
            "null_sd_C_ACC_pct": float(null_std[index]),
            "z_degree_matched": float(z_values[index]),
            "empirical_p_upper": float(empirical_p[index]),
            "q_bh_108": float(q_values[index]),
            "null_draws": n_null,
            "monte_carlo_se": float(
                math.sqrt(
                    empirical_p[index]
                    * (1.0 - empirical_p[index])
                    / (n_null + 1)
                )
            ),
        }
        for index, drug in enumerate(primary_drugs)
    ]
    null_rows.sort(key=lambda row: (row["empirical_p_upper"], row["drug"]))
    primary_index = {drug: index for index, drug in enumerate(primary_drugs)}
    cdk_indices = np.asarray(
        [primary_index[drug] for drug in CDK46_DRUGS],
        dtype=int,
    )
    observed_cdk_mean = float(observed_c_acc_pct[cdk_indices].mean())
    null_cdk_means = null_c_acc_pct[cdk_indices, :].mean(axis=0)
    cdk_null_p = empirical_upper_p(observed_cdk_mean, null_cdk_means)

    random_null_mean = {
        drug: float(null_mean[index]) for index, drug in enumerate(primary_drugs)
    }
    baseline_rows, hypothesis = build_baseline_comparison(
        primary_rows,
        inputs.external_score,
        inputs.associations,
        seed_weights,
        random_null_mean,
        cdk_null_p,
    )

    pagerank, pagerank_iterations, pagerank_delta = weighted_pagerank(graph)
    gene_graph_indices = np.asarray(
        [graph.node_index[gene] for gene in associated_genes],
        dtype=int,
    )
    gene_degree = graph.degree[gene_graph_indices]
    gene_strength = graph.strength[gene_graph_indices]
    gene_pagerank = pagerank[gene_graph_indices]
    bootstrap_rng = np.random.default_rng(RNG_SEED)
    centrality_correlations = {
        "gene_rACC_vs_degree": _correlation_record(
            observed_r_acc, gene_degree, bootstrap_rng
        ),
        "gene_rACC_vs_strength": _correlation_record(
            observed_r_acc, gene_strength, bootstrap_rng
        ),
        "gene_rACC_vs_PageRank": _correlation_record(
            observed_r_acc, gene_pagerank, bootstrap_rng
        ),
    }
    centrality_correlations[
        "gene_rACC_vs_PageRank_partial_log_degree"
    ] = {
        "n": len(associated_genes),
        "rho": partial_spearman(
            observed_r_acc,
            gene_pagerank,
            np.log2(gene_degree + 1.0),
        ),
        "interpretation": "rank residual correlation controlling log2(degree+1)",
    }
    centrality_gene_rows = [
        {
            "gene": gene,
            "r_ACC": float(observed_r_acc[index]),
            "is_ACC_seed": gene in seed_weights,
            "STRING_degree": int(gene_degree[index]),
            "STRING_strength": float(gene_strength[index]),
            "STRING_PageRank": float(gene_pagerank[index]),
            "log2_degree_plus1": float(np.log2(gene_degree[index] + 1.0)),
        }
        for index, gene in enumerate(associated_genes)
    ]

    mean_degree: dict[str, float] = {}
    mean_strength: dict[str, float] = {}
    mean_pagerank: dict[str, float] = {}
    for drug in primary_drugs:
        indices = [graph.node_index[gene] for gene in inputs.associations[drug]]
        mean_degree[drug] = float(graph.degree[indices].mean())
        mean_strength[drug] = float(graph.strength[indices].mean())
        mean_pagerank[drug] = float(pagerank[indices].mean())
    primary_row_map = {str(row["drug"]): row for row in primary_rows}
    centrality_drug_rows = [
        {
            "drug": drug,
            "n_assoc": int(primary_row_map[drug]["n_assoc"]),
            "mean_target_STRING_degree": mean_degree[drug],
            "mean_target_STRING_strength": mean_strength[drug],
            "mean_target_STRING_PageRank": mean_pagerank[drug],
            "C_ACC_pct": float(primary_row_map[drug]["C_ACC_pct"]),
            "ADRS_comp": float(primary_row_map[drug]["ADRS_comp"]),
        }
        for drug in primary_drugs
    ]
    drug_vectors = {
        "n_assoc": np.asarray(
            [primary_row_map[drug]["n_assoc"] for drug in primary_drugs],
            dtype=float,
        ),
        "mean_degree": np.asarray(
            [mean_degree[drug] for drug in primary_drugs], dtype=float
        ),
        "mean_strength": np.asarray(
            [mean_strength[drug] for drug in primary_drugs], dtype=float
        ),
        "mean_PageRank": np.asarray(
            [mean_pagerank[drug] for drug in primary_drugs], dtype=float
        ),
    }
    c_acc_vector = np.asarray(
        [primary_row_map[drug]["C_ACC_pct"] for drug in primary_drugs],
        dtype=float,
    )
    adrs_vector = np.asarray(
        [primary_row_map[drug]["ADRS_comp"] for drug in primary_drugs],
        dtype=float,
    )
    for name, vector in drug_vectors.items():
        centrality_correlations[f"drug_C_ACC_vs_{name}"] = _correlation_record(
            c_acc_vector, vector, bootstrap_rng
        )
        centrality_correlations[f"drug_ADRS_comp_vs_{name}"] = (
            _correlation_record(adrs_vector, vector, bootstrap_rng)
        )

    degree_rho = abs(
        centrality_correlations["gene_rACC_vs_degree"]["rho"]
    )
    pagerank_rho = abs(
        centrality_correlations["gene_rACC_vs_PageRank"]["rho"]
    )
    bh_alpha = 0.05
    bh_family_size = len(primary_drugs)
    minimum_null_draws_for_bh = math.floor(bh_family_size / bh_alpha)
    bh_q_minimum_possible = min(
        1.0,
        bh_family_size / (n_null + 1.0),
    )
    bh_resolution_adequate = bh_q_minimum_possible < bh_alpha
    any_drug_q = bool(np.any(q_values < bh_alpha))
    centrality_below_threshold = degree_rho < 0.70 and pagerank_rho < 0.70
    null_signal = any_drug_q or cdk_null_p < 0.05
    if centrality_below_threshold and null_signal:
        h2_status = "supported"
    elif centrality_below_threshold or null_signal:
        h2_status = "partially_supported"
    else:
        h2_status = "not_supported"
    hypothesis["H2_disease_context_beyond_centrality"] = {
        "status": h2_status,
        "abs_rho_rACC_degree": degree_rho,
        "abs_rho_rACC_PageRank": pagerank_rho,
        "n_drugs_q_lt_0_05": int(np.count_nonzero(q_values < 0.05)),
        "CDK46_degree_matched_seed_p": cdk_null_p,
        "rule": (
            "both centrality |rho|<0.70 and at least one drug q<0.05 or "
            "CDK4/6 null P<0.05"
        ),
    }

    metrics: dict[str, Any] = {
        "analysis_version": ANALYSIS_VERSION,
        "primary_universe_n": len(primary_drugs),
        "associated_gene_n": len(associated_genes),
        "all_drug_n": len(all_drugs),
        "ACC_seed_n": len(seed_weights),
        "STRING_threshold": STRING_THRESHOLD,
        "STRING_symbol_nodes": len(graph.node_names),
        "STRING_source_edge_rows": graph.source_edge_rows,
        "STRING_aggregated_undirected_edges": int(graph.adjacency.nnz // 2),
        "RWR_restart_probability": RWR_RESTART,
        "RWR_observed_iterations": observed_iterations,
        "RWR_observed_final_max_L1_delta": observed_delta,
        "rACC_frozen_max_abs_difference": r_acc_max_abs_difference,
        "rACC_frozen_spearman": r_acc_spearman,
        "PageRank_damping": PAGERANK_DAMPING,
        "PageRank_iterations": pagerank_iterations,
        "PageRank_final_L1_delta": pagerank_delta,
        "degree_bins_requested": REQUESTED_DEGREE_BINS,
        "degree_bins_effective": len(degree_bin_edges) - 1,
        "degree_bin_edges_log2_degree_plus1": degree_bin_edges.tolist(),
        "null_draws": n_null,
        "null_rng_seed": RNG_SEED,
        "null_batch_size": null_batch_size,
        "null_RWR_iterations_min": min(null_iteration_counts),
        "null_RWR_iterations_max": max(null_iteration_counts),
        "null_RWR_final_delta_max": max(null_final_deltas),
        "empirical_p_minimum_resolution": 1.0 / (n_null + 1),
        "BH_family_size": bh_family_size,
        "BH_alpha": bh_alpha,
        "BH_q_minimum_possible": bh_q_minimum_possible,
        "minimum_null_draws_for_any_BH_q_lt_0_05": (
            minimum_null_draws_for_bh
        ),
        "BH_resolution_adequate": bh_resolution_adequate,
        "CDK46_observed_mean_C_ACC_pct": observed_cdk_mean,
        "CDK46_degree_matched_null_mean": float(null_cdk_means.mean()),
        "CDK46_degree_matched_null_sd": float(
            null_cdk_means.std(ddof=1)
        ),
        "CDK46_degree_matched_empirical_p": cdk_null_p,
        "centrality_correlations": centrality_correlations,
        "hypothesis_decisions": hypothesis,
        "clinical_performance_gate": {
            "auc_estimable": False,
            "strict_positive_n": 2,
            "strict_negative_n": 0,
            "reason": (
                "C4 strict drug-specific clinical subset has no unambiguous "
                "negative comparator"
            ),
        },
    }
    return StrengtheningResult(
        baseline_rows=tuple(baseline_rows),
        centrality_gene_rows=tuple(centrality_gene_rows),
        centrality_drug_rows=tuple(centrality_drug_rows),
        matched_seed_rows=tuple(matched_seed_rows),
        null_rows=tuple(null_rows),
        metrics=metrics,
    )


def _write_claim_table(path: Path, result: StrengtheningResult) -> None:
    metrics = result.metrics
    decisions = metrics["hypothesis_decisions"]
    adrs = next(
        row for row in result.baseline_rows if row["ranking"] == "ADRS_comp"
    )
    strongest = max(
        (row for row in result.baseline_rows if row["ranking"] != "ADRS_comp"),
        key=lambda row: abs(row["spearman_vs_ADRS_comp"]),
    )
    centrality = metrics["centrality_correlations"]
    path.write_text(
        "\n".join(
            [
                "# Method-strengthening claim–evidence table",
                "",
                "| Claim | Direct evidence | Statistical support | Boundary |",
                "|---|---|---|---|",
                (
                    "| ADRS_comp composition audit "
                    f"({decisions['composition_audit']['status']}) "
                    f"| strongest single component/baseline: {strongest['ranking']}, "
                    f"rho={strongest['spearman_vs_ADRS_comp']:.3f}, "
                    f"Top-20 Jaccard={strongest['top20_jaccard']:.3f} "
                    "| Pre-frozen joint threshold: abs(rho)≥0.90 and Jaccard≥0.80 "
                    "| Descriptive only; former H1 is retired because the "
                    "equal-weight score definition makes component similarity "
                    "unsuitable as an inferential non-redundancy claim |"
                ),
                (
                    "| r_ACC beyond centrality "
                    f"({decisions['H2_disease_context_beyond_centrality']['status']}) "
                    f"| rho with degree="
                    f"{centrality['gene_rACC_vs_degree']['rho']:.3f}; "
                    f"rho with PageRank="
                    f"{centrality['gene_rACC_vs_PageRank']['rho']:.3f} "
                    f"| CDK4/6 degree-matched seed P="
                    f"{metrics['CDK46_degree_matched_empirical_p']:.4f}; "
                    f"drug q<0.05 n="
                    f"{decisions['H2_disease_context_beyond_centrality']['n_drugs_q_lt_0_05']} "
                    "| Network nodes are dependent; centrality correlation is descriptive |"
                ),
                (
                    "| CDK4/6 ranking robustness "
                    f"({decisions['H3_CDK46_robustness']['status']}) "
                    f"| ADRS ranks {adrs['cdk46_ranks']}, exact P="
                    f"{adrs['cdk46_p_exact']:.4f} "
                    f"| degree-matched seed group P="
                    f"{metrics['CDK46_degree_matched_empirical_p']:.4f} "
                    "| A non-significant trend cannot establish efficacy |"
                ),
                (
                    "| Clinical predictive performance is not estimable "
                    "| strict clinical subset: 2 positive, 0 negative "
                    "| ROC-AUC, PR-AUC and calibration not computed "
                    "| No pseudo-negative or mixed-evidence benchmark is introduced |"
                ),
                "",
            ]
        ),
        encoding="utf-8",
    )


def _write_report(path: Path, result: StrengtheningResult) -> None:
    metrics = result.metrics
    decisions = metrics["hypothesis_decisions"]
    centrality = metrics["centrality_correlations"]
    baseline_lines = [
        (
            f"| {row['ranking']} | {row['spearman_vs_ADRS_comp']:.3f} | "
            f"{row['top20_jaccard']:.3f} | {row['cdk46_mean_rank']:.2f} | "
            f"{row['cdk46_p_exact']:.4f} | "
            f"{row['cdk46_q_bh_across_rankings']:.4f} |"
        )
        for row in result.baseline_rows
    ]
    significant_drugs = [
        row for row in result.null_rows if row["q_bh_108"] < 0.05
    ]
    significant_text = (
        ", ".join(
            f"{row['drug']} (Z={row['z_degree_matched']:.2f}, "
            f"q={row['q_bh_108']:.4f})"
            for row in significant_drugs
        )
        if significant_drugs
        else "None"
    )
    path.write_text(
        "\n".join(
            [
                "# Method strengthening report",
                "",
                "## Design lock",
                "",
                f"- Version: `{metrics['analysis_version']}`.",
                f"- Primary universe: {metrics['primary_universe_n']} drugs.",
                f"- Full STRING: {metrics['STRING_symbol_nodes']:,} symbols and "
                f"{metrics['STRING_aggregated_undirected_edges']:,} aggregated "
                f"undirected edges at combined score ≥{metrics['STRING_threshold']}.",
                f"- Disease-only ACC seeds: {metrics['ACC_seed_n']}.",
                f"- Degree-matched null: {metrics['null_draws']:,} draws; "
                f"RNG seed {metrics['null_rng_seed']}; empirical P floor "
                f"{metrics['empirical_p_minimum_resolution']:.6f}.",
                f"- Best-case BH q floor across "
                f"{metrics['BH_family_size']} drugs: "
                f"{metrics['BH_q_minimum_possible']:.6f}; resolution adequate "
                f"for q<0.05: {metrics['BH_resolution_adequate']}.",
                "- Clinical AUC/PR-AUC/calibration: not estimable (2 strict "
                "positives, 0 strict negatives).",
                "",
                "## 1. Simple baselines",
                "",
                "| Ranking | Spearman vs ADRS | Top-20 Jaccard | CDK4/6 mean rank | Exact P | BH q across rankings |",
                "|---|---:|---:|---:|---:|---:|",
                *baseline_lines,
                "",
                "Former H1 decision: **retired**. Correlation and Top-20 overlap "
                "are retained only as a descriptive composition audit because "
                "ADRS_comp is an equal-weight sum of two ranked components.",
                "",
                "## 2. STRING centrality",
                "",
                f"- `r_ACC` vs degree: rho="
                f"{centrality['gene_rACC_vs_degree']['rho']:.3f}, 95% bootstrap "
                f"CI [{centrality['gene_rACC_vs_degree']['bootstrap_ci_95'][0]:.3f}, "
                f"{centrality['gene_rACC_vs_degree']['bootstrap_ci_95'][1]:.3f}].",
                f"- `r_ACC` vs strength: rho="
                f"{centrality['gene_rACC_vs_strength']['rho']:.3f}.",
                f"- `r_ACC` vs PageRank: rho="
                f"{centrality['gene_rACC_vs_PageRank']['rho']:.3f}, 95% bootstrap "
                f"CI [{centrality['gene_rACC_vs_PageRank']['bootstrap_ci_95'][0]:.3f}, "
                f"{centrality['gene_rACC_vs_PageRank']['bootstrap_ci_95'][1]:.3f}].",
                f"- Partial Spearman with PageRank controlling log-degree: rho="
                f"{centrality['gene_rACC_vs_PageRank_partial_log_degree']['rho']:.3f}.",
                "",
                "Bootstrap intervals are descriptive because network nodes are not "
                "independent observational units.",
                "",
                "## 3. Degree-matched random-seed null",
                "",
                f"- Drugs with BH q<0.05: {len(significant_drugs)}.",
                f"- Significant drugs: {significant_text}.",
                f"- CDK4/6 observed mean C_ACC percentile: "
                f"{metrics['CDK46_observed_mean_C_ACC_pct']:.3f}.",
                f"- CDK4/6 null mean ± SD: "
                f"{metrics['CDK46_degree_matched_null_mean']:.3f} ± "
                f"{metrics['CDK46_degree_matched_null_sd']:.3f}; empirical "
                f"P={metrics['CDK46_degree_matched_empirical_p']:.4f}.",
                f"- H2 decision: **"
                f"{decisions['H2_disease_context_beyond_centrality']['status']}**.",
                f"- H3 decision: **"
                f"{decisions['H3_CDK46_robustness']['status']}**.",
                "",
                "## 4. Quality control",
                "",
                f"- Recomputed vs frozen r_ACC: maximum absolute difference "
                f"{metrics['rACC_frozen_max_abs_difference']:.3g}; Spearman "
                f"{metrics['rACC_frozen_spearman']:.9f}.",
                f"- Effective degree bins: {metrics['degree_bins_effective']} "
                f"(requested {metrics['degree_bins_requested']}).",
                f"- Null RWR iterations: {metrics['null_RWR_iterations_min']}–"
                f"{metrics['null_RWR_iterations_max']}; maximum final L1 delta "
                f"{metrics['null_RWR_final_delta_max']:.3g}.",
                "",
                "## 5. Interpretation boundary",
                "",
                "These analyses can show non-redundant ranking structure and quantify "
                "network/seed dependence. They cannot establish drug efficacy, causal "
                "targets or prospective generalization. `S_external` remains a "
                "non-independent exploratory baseline because the same literature "
                "informed evidence-aware reprioritization.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def write_outputs(
    project_root: Path,
    result: StrengtheningResult,
    output_dir: Path,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "baselines": output_dir / "baseline_comparison_primary108.csv",
        "centrality_gene": output_dir / "centrality_gene399.csv",
        "centrality_drug": output_dir / "centrality_drug108.csv",
        "matched_seeds": output_dir / "degree_matched_seed_sets.csv",
        "null_summary": output_dir / "random_seed_null_primary108.csv",
        "metrics": output_dir / "method_strengthening_metrics.json",
        "claims": output_dir / "claim_evidence_table.md",
        "report": output_dir / "method_strengthening_report.md",
        "manifest": output_dir / "run_manifest.md",
    }
    _write_csv(paths["baselines"], result.baseline_rows)
    _write_csv(paths["centrality_gene"], result.centrality_gene_rows)
    _write_csv(paths["centrality_drug"], result.centrality_drug_rows)
    _write_csv(paths["matched_seeds"], result.matched_seed_rows)
    _write_csv(paths["null_summary"], result.null_rows)
    paths["metrics"].write_text(
        json.dumps(result.metrics, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    _write_claim_table(paths["claims"], result)
    _write_report(paths["report"], result)

    input_paths = (
        project_root / "data" / "bindex_network" / "bindex_edges_1304.csv",
        project_root / "data" / "bindex_network" / "rACC_399_fullSTRING.csv",
        project_root / "data" / "bindex_network" / "Sactivity_124_v1.csv",
        project_root / "data" / "bindex_network" / "NCI60_potency_124.csv",
        project_root / "data" / "bindex_network" / "S_external_curated.csv",
        project_root / "data" / "ACC_P0.5C_gene_weights_v1.csv",
        project_root / "9606.protein.info.v12.0.txt.gz",
        project_root / "9606.protein.links.v12.0.txt.gz",
    )
    input_lines = [
        f"- `{path.relative_to(project_root).as_posix()}`: `{_sha256(path)}`"
        for path in input_paths
    ]
    output_lines = [
        f"- `{path.name}`: `{_sha256(path)}`"
        for name, path in paths.items()
        if name != "manifest"
    ]
    paths["manifest"].write_text(
        "\n".join(
            [
                "# Method-strengthening run manifest",
                "",
                f"- Analysis version: `{ANALYSIS_VERSION}`",
                "- Command: `python -m analysis.method_strengthening --project-root .`",
                f"- Python: `{platform.python_version()}`",
                f"- NumPy: `{np.__version__}`",
                f"- SciPy: `{scipy.__version__}`",
                f"- Platform: `{platform.platform()}`",
                f"- RNG seed: `{result.metrics['null_rng_seed']}`",
                f"- Null draws: `{result.metrics['null_draws']}`",
                "",
                "## Input SHA-256",
                "",
                *input_lines,
                "",
                "## Output SHA-256",
                "",
                *output_lines,
                "",
                "## Guardrails",
                "",
                "- The C1 primary score is read and compared, not refitted.",
                "- No AUC, PR-AUC or calibration is computed without negative "
                "clinical comparators.",
                "- Random-seed empirical P values use add-one correction and BH "
                "adjustment across 108 drugs.",
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
) -> StrengtheningResult:
    root = project_root.resolve()
    result = run_analysis(root, n_null=n_null, null_batch_size=null_batch_size)
    target = (
        output_dir.resolve()
        if output_dir is not None
        else root / "results" / "method_strengthening"
    )
    write_outputs(root, result, target)
    return result


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
    result = run(
        args.project_root,
        output_dir=args.output_dir,
        n_null=args.null_draws,
        null_batch_size=args.null_batch_size,
    )
    decisions = result.metrics["hypothesis_decisions"]
    print(
        json.dumps(
            {
                "status": "ok",
                "analysis_version": ANALYSIS_VERSION,
                "null_draws": result.metrics["null_draws"],
                "H1": decisions["H1_nonredundancy"]["status"],
                "H2": decisions[
                    "H2_disease_context_beyond_centrality"
                ]["status"],
                "H3": decisions["H3_CDK46_robustness"]["status"],
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
