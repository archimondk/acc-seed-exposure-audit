from __future__ import annotations

from pathlib import Path

import pytest

from analysis.leakage_audit import (
    build_arm_seed_sets,
    evaluate_leakage_verdict,
    verify_frozen_protocol,
)


ROOT = Path(__file__).resolve().parents[1]


def test_frozen_protocol_matches_external_hash_and_locked_thresholds() -> None:
    record = verify_frozen_protocol(ROOT)

    assert record["protocol_id"] == "leakage_audit_v1"
    assert record["protocol_sha256"] == (
        "36c9638ade80bb761f6e8481889575b1b80feb460bbd51c9d5d68617a4155e85"
    )
    assert record["l2_threshold"] == pytest.approx(1.5)
    assert record["l3_threshold"] == pytest.approx(2.0)
    assert record["interventional_results_observed_at_freeze"] is False


def test_arm_construction_is_symmetric_and_uses_fixed_b2_anchor() -> None:
    arms = build_arm_seed_sets(ROOT)

    assert set(arms) == {"A1", "A2", "B1", "B2", "B2_lo", "B2_hi"}
    assert len(arms["A1"]) == 45
    assert len(arms["A2"]) == 44
    assert len(arms["B1"]) == 24
    assert len(arms["B2"]) == 25
    assert "RB1" in arms["A1"] and "RB1" not in arms["A2"]
    assert "RB1" not in arms["B1"] and "RB1" in arms["B2"]

    breast_median = sorted(arms["B1"].values())[len(arms["B1"]) // 2 - 1 :][
        :2
    ]
    breast_median_value = sum(breast_median) / 2.0
    assert arms["B2"]["RB1"] == pytest.approx(breast_median_value)
    assert arms["B2_lo"]["RB1"] == pytest.approx(0.5 * breast_median_value)
    assert arms["B2_hi"]["RB1"] == pytest.approx(1.5 * breast_median_value)


def test_leakage_verdict_requires_all_frozen_rules_and_honors_falsification() -> None:
    passing = {
        variant: {
            "z_A1": 3.4,
            "z_A2": 1.1,
            "z_B1": 0.2,
            "z_B2": 2.0,
            "ribociclib_delta_acc": 0.1,
            "ribociclib_delta_breast": 0.2,
        }
        for variant in (
            "column_minmax",
            "column_gene_rank",
            "uniform_ratio_gene_rank",
            "symmetric_gene_rank",
        )
    }
    verdict = evaluate_leakage_verdict(passing)
    assert verdict["status"] == "LEAKAGE_SUPPORTED"
    assert verdict["criteria"]["L4"]["n_variants_passing"] == 4

    falsified = {name: dict(values) for name, values in passing.items()}
    falsified["column_minmax"]["z_A2"] = 2.1
    verdict = evaluate_leakage_verdict(falsified)
    assert verdict["status"] == "FALSIFIED_ACC_SIGNAL_SURVIVES"
    assert verdict["criteria"]["F1"]["passed"] is True
