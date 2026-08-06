from pathlib import Path

from analysis.evidence_label_audit import (
    build_audit_rows,
    evaluate_benchmark,
    load_evidence_labels,
    load_primary_scores,
    load_second_review_summary,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_evidence_table_has_complete_unique_rows() -> None:
    rows = load_evidence_labels(
        PROJECT_ROOT
        / "data"
        / "evidence"
        / "evidence_labels_v3_adjudicated.csv"
    )

    assert len(rows) == 19
    assert len({row["drug"] for row in rows}) == 19
    assert sum(row["legacy_benchmark_member"] for row in rows) == 14
    assert all(row["source_verification"] == "primary_source_checked" for row in rows)
    assert all(row["independent_second_review"] != "pending" for row in rows)


def test_strict_binary_clinical_benchmark_is_not_estimable() -> None:
    evidence = load_evidence_labels(
        PROJECT_ROOT
        / "data"
        / "evidence"
        / "evidence_labels_v3_adjudicated.csv"
    )
    second_review = load_second_review_summary(
        PROJECT_ROOT
        / "results"
        / "evidence_audit"
        / "second_reviewer_agreement.json"
    )
    metrics = evaluate_benchmark(evidence, second_review)

    assert metrics["legacy_n"] == 14
    assert metrics["legacy_positive_n"] == 10
    assert metrics["legacy_negative_n"] == 4
    assert metrics["strict_candidate_drugs"] == ["Cabozantinib", "Mitotane"]
    assert metrics["strict_positive_n"] == 2
    assert metrics["strict_negative_n"] == 0
    assert metrics["auc_estimable"] is False
    assert metrics["auc"] is None
    assert metrics["second_reviewer_status"] == "completed_and_adjudicated"
    assert metrics["strict_binary_agreement_n"] == 19
    assert metrics["strict_binary_agreement_total"] == 19
    assert metrics["strict_binary_cohen_kappa"] == 1.0


def test_legacy_figure_rows_use_frozen_primary_ranks() -> None:
    evidence = load_evidence_labels(
        PROJECT_ROOT
        / "data"
        / "evidence"
        / "evidence_labels_v3_adjudicated.csv"
    )
    primary = load_primary_scores(
        PROJECT_ROOT
        / "results"
        / "primary_analysis"
        / "ADRS_comp_primary_108.csv"
    )
    rows = build_audit_rows(evidence, primary, legacy_only=True)
    rank_by_drug = {row["drug"]: row["rank_comp"] for row in rows}

    assert len(rows) == 14
    assert rank_by_drug["Mitotane"] == 17
    assert rank_by_drug["Cabozantinib"] == 82
    assert rank_by_drug["Sunitinib"] == 69
    assert rank_by_drug["Cisplatin"] == 103


def test_external_labels_are_never_claimed_as_independent_validation() -> None:
    rows = load_evidence_labels(
        PROJECT_ROOT
        / "data"
        / "evidence"
        / "evidence_labels_v3_adjudicated.csv"
    )

    reused = [row for row in rows if row["used_in_external_score"]]
    assert len(reused) == 19
    assert not any(row["eligible_as_independent_validation"] for row in reused)
