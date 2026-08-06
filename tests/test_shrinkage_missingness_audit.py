import pytest

from analysis.shrinkage_missingness_audit import (
    classify_missing_drug,
    top_k_jaccard_from_ranks,
)


def test_missing_drug_classification_is_explicit() -> None:
    assert classify_missing_drug("Talazoparib") == ("PARP", "targeted")
    assert classify_missing_drug("Ifosfamide") == (
        "Alkylator/platinum",
        "broad_cytotoxic",
    )


def test_top_k_jaccard_from_ranks() -> None:
    first = {"a": 1, "b": 2, "c": 3}
    second = {"a": 1, "b": 3, "c": 2}

    assert top_k_jaccard_from_ranks(first, second, 2) == pytest.approx(1 / 3)
