from __future__ import annotations

import csv
from pathlib import Path

import pytest

from analysis.positive_control import (
    classify_erpos_her2neg,
    evaluate_success_criteria,
    load_positive_control_seed_weights,
)


def test_classify_erpos_her2neg_uses_fish_then_ihc_fallback() -> None:
    assert classify_erpos_her2neg("Positive", "Equivocal", "Negative")
    assert classify_erpos_her2neg(
        "Positive", "Negative", "[Not Available]"
    )
    assert not classify_erpos_her2neg("Negative", "Negative", "Negative")
    assert not classify_erpos_her2neg("Positive", "Negative", "Positive")
    assert not classify_erpos_her2neg(
        "Positive", "Equivocal", "[Not Available]"
    )


def test_seed_loader_rejects_direct_targets_and_outcome_fields(
    tmp_path: Path,
) -> None:
    path = tmp_path / "seeds.csv"
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "gene",
                "genomic_driver",
                "recurrence",
                "core_pathway",
                "lineage_biomarker",
                "prognostic_subtype",
                "raw_weight",
                "include_primary",
                "exclusion_reason",
                "source_id",
                "source_version",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "gene": "CDK4",
                "genomic_driver": 1,
                "recurrence": 0.5,
                "core_pathway": 0,
                "lineage_biomarker": 0,
                "prognostic_subtype": 0,
                "raw_weight": 0.4,
                "include_primary": "yes",
                "exclusion_reason": "",
                "source_id": "test",
                "source_version": "1",
            }
        )
    with pytest.raises(ValueError, match="direct positive-control target"):
        load_positive_control_seed_weights(path)

    path.write_text(
        "gene,raw_weight,include_primary,drug_rank\n"
        "PIK3CA,1,yes,1\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="outcome-derived"):
        load_positive_control_seed_weights(path)


def test_success_gate_requires_all_four_prespecified_criteria() -> None:
    passing = evaluate_success_criteria(
        primary_group_p=0.01,
        primary_group_q=0.04,
        top_quartile_flags={
            "Abemaciclib": True,
            "Palbociclib": True,
            "Ribociclib": False,
        },
        concordant_variants=3,
    )
    assert passing["status"] == "pass"
    assert all(passing["criteria"].values())

    partial = evaluate_success_criteria(
        primary_group_p=0.01,
        primary_group_q=0.08,
        top_quartile_flags={
            "Abemaciclib": True,
            "Palbociclib": True,
            "Ribociclib": False,
        },
        concordant_variants=3,
    )
    assert partial["status"] == "partial_recovery"
    assert not partial["criteria"]["primary_group_q_lt_0_05"]

