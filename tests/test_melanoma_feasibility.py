"""Synthetic guards: coverage, same-codon ambiguity, and participant duplication."""
import pandas as pd
from analysis.melanoma_feasibility import classify_braf, select_one_sample, normalize_name, external_model_mask


def test_coverage_is_required_for_reference_group():
    assert classify_braf([], sequenced=False) == "not_profiled"
    assert classify_braf([], sequenced=True) == "no_reported_BRAF_variant"
    assert classify_braf([float("nan")], sequenced=True) == "other_BRAF_variant"


def test_exact_variant_not_prefix_match():
    assert classify_braf(["p.V600E"], True) == "V600E"
    assert classify_braf(["V600K"], True) == "V600K"
    assert classify_braf(["V600R"], True) == "other_BRAF_variant"
    assert classify_braf(["V600_K601delinsE"], True) == "ambiguous_V600"


def test_split_codon_is_not_silently_called_e_or_k():
    assert classify_braf(["p.V600E", "p.V600M"], True) == "ambiguous_V600"
    assert classify_braf(["p.V600E", "p.V600="], True) == "V600E"
    assert classify_braf(["p.V600E", "p.S467L"], True) == "V600E"


def test_patient_dedup_is_deterministic_and_scenario_specific():
    frame = pd.DataFrame({"PATIENT_ID": ["P1", "P1", "P2"],
                          "SAMPLE_ID": ["P1-06", "P1-01", "P2-06"],
                          "SAMPLE_TYPE": ["Metastasis", "Primary", "Metastasis"]})
    assert select_one_sample(frame, "all").SAMPLE_ID.tolist() == ["P1-01", "P2-06"]
    assert select_one_sample(frame, "metastatic").SAMPLE_ID.tolist() == ["P1-06", "P2-06"]


def test_name_normalization_does_not_merge_related_lines():
    assert normalize_name("SK-MEL-28") == normalize_name("SK_MEL_28")
    assert normalize_name("MDA-MB-435") != normalize_name("MDA-MB-435S")


def test_overlap_exclusion_survives_nullable_object_merge():
    frame = pd.DataFrame({"QC_eligible": [True, True, True, False],
                          "NCI60_overlap": pd.Series([True, False, None, False], dtype=object)})
    assert external_model_mask(frame).tolist() == [False, True, False, False]
