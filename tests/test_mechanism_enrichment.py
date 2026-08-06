from pathlib import Path

import pytest

from analysis.mechanism_enrichment import (
    benjamini_hochberg,
    compute_mechanism_enrichment,
    exact_lower_tail_rank_sum_p,
    load_analysis_inputs,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_exact_rank_sum_null_uses_the_primary_universe() -> None:
    assert exact_lower_tail_rank_sum_p(108, 3, 85) == pytest.approx(
        0.07635827504457375
    )


def test_bh_adjustment_is_monotone_in_p_value_order() -> None:
    adjusted = benjamini_hochberg([0.04, 0.001, 0.03, 0.2])
    assert adjusted == pytest.approx([0.0533333333, 0.004, 0.0533333333, 0.2])


def test_cdk46_result_is_reproducible_from_c1_ranking() -> None:
    primary_rows, metadata = load_analysis_inputs(PROJECT_ROOT)
    result = compute_mechanism_enrichment(primary_rows, metadata)

    assert result.n_universe == 108
    assert len(result.class_rows) == 10
    assert all(row["n_universe"] == 108 for row in result.class_rows)

    cdk = next(row for row in result.class_rows if row["mechanism_class"] == "CDK4/6")
    assert cdk["k"] == 3
    assert cdk["members"] == "Abemaciclib; Palbociclib; Ribociclib"
    assert cdk["member_ranks"] == "8; 26; 51"
    assert cdk["rank_sum"] == 85
    assert cdk["mean_rank"] == pytest.approx(28.3333333333)
    assert cdk["p_exact"] == pytest.approx(0.07635827504457375)
    assert cdk["q_bh"] == pytest.approx(0.37110183715720657)


def test_all_observed_members_are_inside_the_same_null_universe() -> None:
    primary_rows, metadata = load_analysis_inputs(PROJECT_ROOT)
    result = compute_mechanism_enrichment(primary_rows, metadata)
    primary_drugs = {row["drug"] for row in primary_rows}

    assert len(result.member_rows) == 108
    assert {row["drug"] for row in result.member_rows} == primary_drugs
    assert all(row["in_primary_universe"] for row in result.member_rows)
