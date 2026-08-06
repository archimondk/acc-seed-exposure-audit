from pathlib import Path

import csv
import hashlib
import json

import numpy as np
import pytest

from analysis.seed_weight_assignment_sensitivity import (
    DISEASE_SEED_COUNT,
    make_permuted_seed_weights,
    make_uniform_seed_weights,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_uniform_weights_are_exact_membership_only_simplex() -> None:
    weights = make_uniform_seed_weights(DISEASE_SEED_COUNT)

    assert weights.shape == (DISEASE_SEED_COUNT, 1)
    assert weights[:, 0] == pytest.approx(
        np.full(DISEASE_SEED_COUNT, 1.0 / DISEASE_SEED_COUNT),
        abs=1e-15,
    )
    assert weights[:, 0].sum() == pytest.approx(1.0, abs=1e-15)


def test_permutations_are_deterministic_and_preserve_weight_multiset() -> None:
    baseline = np.asarray([0.1, 0.2, 0.3, 0.4])

    first = make_permuted_seed_weights(baseline, n_draws=12, rng_seed=17)
    second = make_permuted_seed_weights(baseline, n_draws=12, rng_seed=17)

    assert np.array_equal(first, second)
    assert first.shape == (4, 12)
    assert first.sum(axis=0) == pytest.approx(np.ones(12), abs=1e-15)
    for column in range(first.shape[1]):
        assert np.array_equal(
            np.sort(first[:, column]),
            np.sort(baseline),
        )


def test_permutation_inputs_are_validated() -> None:
    with pytest.raises(ValueError, match="positive"):
        make_permuted_seed_weights(
            np.asarray([0.0, 0.5, 0.5]),
            n_draws=2,
            rng_seed=1,
        )
    with pytest.raises(ValueError, match="sum to one"):
        make_permuted_seed_weights(
            np.asarray([0.2, 0.2, 0.2]),
            n_draws=2,
            rng_seed=1,
        )


def test_frozen_protocol_hash_and_formal_outputs() -> None:
    protocol = (
        PROJECT_ROOT
        / "experiments"
        / "seed_weight_assignment_sensitivity_protocol_v1.md"
    )
    freeze = (
        PROJECT_ROOT
        / "experiments"
        / "SEED_WEIGHT_ASSIGNMENT_SENSITIVITY_FREEZE.txt"
    )
    protocol_hash = hashlib.sha256(protocol.read_bytes()).hexdigest()

    assert protocol_hash in freeze.read_text(encoding="utf-8")

    output_dir = PROJECT_ROOT / "results" / "seed_weight_assignment_sensitivity"
    summary = json.loads(
        (output_dir / "seed_weight_assignment_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert summary["permutation_draws"] == 1000
    assert summary["rng_seed"] == 20260729
    assert summary["disease_seed_count"] == 45
    assert summary["frozen_verdict_revised"] is False
    assert summary["quality_control"]["all_permutations_preserve_weights"] is True
    assert summary["quality_control"]["all_rank_columns_complete"] is True

    expected_rows = {
        "uniform_drug_ranks.csv": 108,
        "permutation_draw_summary.csv": 1000,
        "permutation_drug_rank_draws.csv": 108_000,
        "permutation_drug_rank_summary.csv": 108,
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
            / f"FigS4_seed_weight_assignment_sensitivity.{suffix}"
        ).is_file()
