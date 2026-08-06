import csv
import json
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

from analysis.method_strengthening import (
    build_baseline_comparison,
    compute_c_acc_matrix,
    empirical_upper_p,
    exact_lower_tail_rank_sum,
    generate_degree_matched_seed_sets,
    load_disease_seed_weights,
    partial_spearman,
)
from analysis.acc_primary_pipeline import compute_primary_analysis, load_inputs


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_empirical_upper_p_uses_add_one_correction() -> None:
    null = np.asarray([0.1, 0.3, 0.4])

    assert empirical_upper_p(0.3, null) == pytest.approx(0.75)
    assert empirical_upper_p(1.0, null) == pytest.approx(0.25)


def test_exact_lower_tail_rank_sum_accepts_average_ranks() -> None:
    ranks = np.asarray([1.0, 2.0, 3.0, 4.0])

    result = exact_lower_tail_rank_sum(ranks, (0, 1))

    assert result == pytest.approx(1 / 6)


def test_degree_matched_seed_sets_preserve_bins_and_exclude_true_seeds() -> None:
    names = ("seed_low", "candidate_low", "candidate_low_2",
             "seed_high", "candidate_high", "candidate_high_2")
    degrees = np.asarray([1, 1, 2, 10, 9, 10], dtype=float)
    seed_weights = {"seed_low": 0.7, "seed_high": 0.3}

    rows, _ = generate_degree_matched_seed_sets(
        seed_weights,
        names,
        degrees,
        n_draws=8,
        rng_seed=20260727,
        n_bins=2,
    )

    assert len(rows) == 16
    assert {row["matched_seed"] for row in rows} == set(seed_weights)
    assert not ({row["null_seed"] for row in rows} & set(seed_weights))
    assert all(
        row["matched_seed_degree_bin"] == row["null_seed_degree_bin"]
        for row in rows
    )
    assert all(row["weight"] == seed_weights[row["matched_seed"]] for row in rows)


def test_compute_c_acc_matrix_matches_hand_calculation() -> None:
    association_matrix = sparse.csr_matrix(
        np.asarray([[1, 0], [0, 1], [1, 1]], dtype=float)
    )
    propagated = np.asarray([[0.2], [0.8]], dtype=float)

    result = compute_c_acc_matrix(
        association_matrix,
        propagated,
        pseudo_count=3.0,
    )

    assert result[:, 0] == pytest.approx([0.425, 0.575, 0.5])


def test_partial_spearman_is_symmetric_and_bounded() -> None:
    x = np.asarray([1, 2, 4, 3, 6, 5, 8, 7], dtype=float)
    y = np.asarray([2, 1, 3, 5, 4, 7, 6, 8], dtype=float)
    control = np.arange(1, 9, dtype=float)

    xy = partial_spearman(x, y, control)
    yx = partial_spearman(y, x, control)

    assert np.isfinite(xy)
    assert -1.0 <= xy <= 1.0
    assert xy == pytest.approx(yx)


def test_primary_baseline_uses_locked_c1_ordinal_ranks() -> None:
    inputs = load_inputs(PROJECT_ROOT)
    primary = compute_primary_analysis(inputs)
    seed_weights = load_disease_seed_weights(
        PROJECT_ROOT / "data" / "ACC_P0.5C_gene_weights_v1.csv"
    )
    random_null_mean = {
        row["drug"]: float(row["C_ACC_pct"]) for row in primary.primary_rows
    }

    rows, _ = build_baseline_comparison(
        primary.primary_rows,
        inputs.external_score,
        inputs.associations,
        seed_weights,
        random_null_mean,
        cdk_null_p=0.5,
    )

    adrs = next(row for row in rows if row["ranking"] == "ADRS_comp")
    assert adrs["cdk46_ranks"] == "8; 26; 51"
    assert adrs["cdk46_p_exact"] == pytest.approx(0.07635827504457375)


def test_formal_method_strengthening_outputs_lock_negative_results() -> None:
    output_dir = PROJECT_ROOT / "results" / "method_strengthening"
    metrics = json.loads(
        (output_dir / "method_strengthening_metrics.json").read_text(
            encoding="utf-8"
        )
    )

    assert metrics["null_draws"] == 10_000
    assert metrics["BH_resolution_adequate"] is True
    assert metrics["BH_q_minimum_possible"] < 0.05
    assert metrics["minimum_null_draws_for_any_BH_q_lt_0_05"] == 2160
    assert "runtime_seconds" not in metrics
    assert metrics["rACC_frozen_max_abs_difference"] <= 1.1e-6
    assert metrics["rACC_frozen_spearman"] >= 0.999999
    assert metrics["clinical_performance_gate"]["auc_estimable"] is False
    decisions = metrics["hypothesis_decisions"]
    assert (
        decisions["H1_nonredundancy"]["status"]
        == "retired_descriptive_only"
    )
    assert (
        decisions["H2_disease_context_beyond_centrality"]["status"]
        == "partially_supported"
    )
    assert decisions["H3_CDK46_robustness"]["status"] == "not_supported"
    assert metrics["CDK46_degree_matched_empirical_p"] == pytest.approx(
        0.28387161283871615
    )


def test_formal_degree_matching_and_empirical_p_integrity() -> None:
    output_dir = PROJECT_ROOT / "results" / "method_strengthening"
    with (output_dir / "degree_matched_seed_sets.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        matched_rows = list(csv.DictReader(stream))
    with (output_dir / "random_seed_null_primary108.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        null_rows = list(csv.DictReader(stream))

    assert len(matched_rows) == 450_000
    assert len({int(row["replicate"]) for row in matched_rows}) == 10_000
    assert all(
        row["matched_seed_degree_bin"] == row["null_seed_degree_bin"]
        for row in matched_rows
    )
    assert len(null_rows) == 108
    assert min(
        float(row["empirical_p_upper"]) for row in null_rows
    ) == pytest.approx(1 / 10001)
    assert sum(float(row["q_bh_108"]) < 0.05 for row in null_rows) == 2
