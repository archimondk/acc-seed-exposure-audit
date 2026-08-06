from pathlib import Path

import csv
import hashlib
import json
import numpy as np
import pytest

from analysis.dirichlet_weight_sensitivity import (
    ACTIVE_COMPONENTS,
    compute_seed_weight_matrix,
    ordinal_rank_columns,
    sample_component_weights,
    top_k_jaccard_columns,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_component_draws_are_deterministic_positive_simplex() -> None:
    first = sample_component_weights(n_draws=8, rng_seed=20260729)
    second = sample_component_weights(n_draws=8, rng_seed=20260729)

    assert first.shape == (8, len(ACTIVE_COMPONENTS))
    assert np.array_equal(first, second)
    assert np.all(np.isfinite(first))
    assert np.all(first > 0)
    assert first.sum(axis=1) == pytest.approx(np.ones(8), abs=1e-12)


def test_seed_weight_matrix_applies_components_then_normalizes() -> None:
    component_matrix = np.asarray(
        [
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0],
        ]
    )
    component_draws = np.asarray(
        [[0.25, 0.75, 1e-9, 1e-9, 1e-9]]
    )

    result = compute_seed_weight_matrix(component_matrix, component_draws)

    assert result.shape == (2, 1)
    assert result[:, 0] == pytest.approx([0.25, 0.75], abs=1e-8)
    assert result.sum(axis=0) == pytest.approx([1.0], abs=1e-12)


def test_ordinal_rank_columns_breaks_ties_by_drug_name() -> None:
    names = ("DrugB", "DrugA", "DrugC")
    scores = np.asarray(
        [
            [1.0, 0.0],
            [1.0, 2.0],
            [0.0, 1.0],
        ]
    )

    ranks = ordinal_rank_columns(scores, names)

    assert ranks[:, 0].tolist() == [2, 1, 3]
    assert ranks[:, 1].tolist() == [3, 1, 2]
    assert sorted(ranks[:, 0].tolist()) == [1, 2, 3]
    assert sorted(ranks[:, 1].tolist()) == [1, 2, 3]


def test_top_k_jaccard_columns_matches_hand_calculation() -> None:
    baseline = np.asarray([1, 2, 3, 4])
    draws = np.asarray(
        [
            [1, 1],
            [2, 3],
            [3, 2],
            [4, 4],
        ]
    )

    result = top_k_jaccard_columns(draws, baseline, k=2)

    assert result.tolist() == pytest.approx([1.0, 1.0 / 3.0])


def test_frozen_protocol_hash_and_formal_outputs() -> None:
    protocol = (
        PROJECT_ROOT
        / "experiments"
        / "dirichlet_component_weight_sensitivity_protocol_v1.md"
    )
    freeze = (
        PROJECT_ROOT
        / "experiments"
        / "DIRICHLET_WEIGHT_SENSITIVITY_FREEZE.txt"
    )
    protocol_hash = hashlib.sha256(protocol.read_bytes()).hexdigest()

    assert protocol_hash in freeze.read_text(encoding="utf-8")

    output_dir = PROJECT_ROOT / "results" / "dirichlet_weight_sensitivity"
    summary = json.loads(
        (output_dir / "dirichlet_weight_sensitivity_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert summary["draws"] == 1000
    assert summary["rng_seed"] == 20260729
    assert summary["active_components"] == list(ACTIVE_COMPONENTS)
    assert summary["disease_seed_count"] == 45
    assert summary["therapeutic_component_included"] is False
    assert summary["frozen_verdict_revised"] is False
    assert summary["quality_control"]["all_rank_columns_complete"] is True
    assert summary["baseline_reproduction"]["rACC_max_abs_difference"] <= 1.1e-6
    assert summary["ADRS_rank_spearman_vs_locked"]["q05"] > 0.99
    assert summary["top20_jaccard_vs_locked"]["median"] == pytest.approx(1.0)
    assert summary["CDK46_drugs"]["Abemaciclib"]["prob_top20"] == pytest.approx(1.0)

    expected_rows = {
        "component_weight_draws.csv": 1000,
        "drug_rank_draws.csv": 108_000,
        "drug_rank_summary.csv": 108,
        "draw_summary.csv": 1000,
    }
    for name, expected in expected_rows.items():
        with (output_dir / name).open(
            "r", encoding="utf-8-sig", newline=""
        ) as stream:
            assert sum(1 for _ in csv.DictReader(stream)) == expected

    for suffix in ("png", "pdf", "svg"):
        assert (
            PROJECT_ROOT
            / "figures"
            / "revision"
            / f"FigS3_dirichlet_weight_sensitivity.{suffix}"
        ).is_file()
