from pathlib import Path

import pytest

from scripts.reproduce import (
    EXPECTED_FIGURES,
    REQUIRED_INPUTS,
    SCIENTIFIC_OUTPUTS,
    compare_reproduction,
    copy_input_manifest,
    resolve_run_root,
    sha256_file,
    validate_frozen_inputs,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_required_inputs_are_portable_and_complete() -> None:
    expected = {
        "data/bindex_network/bindex_edges_1304.csv",
        "data/bindex_network/rACC_399_fullSTRING.csv",
        "data/bindex_network/Sactivity_124_v1.csv",
        "data/bindex_network/NCI60_potency_124.csv",
        "data/bindex_network/S_external_curated.csv",
        "data/ACC_MIPE5_ranked_by_mean_ZAUC.csv",
        "data/ACC_P0.5C_gene_weights_v1.csv",
        "data/evidence/evidence_labels_v2.csv",
        "independent_review/claude_blinded_evidence_review_v1/claude_review_locked.csv",
        "independent_review/cabozantinib_b02_rerate/claude_b02_rerated.csv",
        "9606.protein.info.v12.0.txt.gz",
        "9606.protein.links.v12.0.txt.gz",
        "experiments/amendment4_leave_one_seed_out_protocol_v1.md",
        "experiments/amendment5_seed_excluded_scoring_protocol_v1.md",
    }

    assert set(REQUIRED_INPUTS) == expected
    assert all(not Path(path).is_absolute() for path in REQUIRED_INPUTS)
    assert all((PROJECT_ROOT / path).is_file() for path in REQUIRED_INPUTS)


def test_method_strengthening_outputs_are_in_reproduction_contract() -> None:
    expected_scientific = {
        "results/method_strengthening/baseline_comparison_primary108.csv",
        "results/method_strengthening/centrality_gene399.csv",
        "results/method_strengthening/centrality_drug108.csv",
        "results/method_strengthening/degree_matched_seed_sets.csv",
        "results/method_strengthening/random_seed_null_primary108.csv",
        "results/method_strengthening/method_strengthening_metrics.json",
        "figure_data/revision/Fig2d_random_seed_top12_primary108.csv",
    }
    expected_figures = {
        f"figures/revision/Fig2_method_strengthening_primary108.{suffix}"
        for suffix in ("pdf", "svg", "png")
    }

    assert expected_scientific <= set(SCIENTIFIC_OUTPUTS)
    assert expected_figures <= set(EXPECTED_FIGURES)


def test_major_revision_sensitivity_outputs_are_in_reproduction_contract() -> None:
    expected = {
        "results/normalization_sensitivity/normalization_sensitivity_summary.csv",
        "results/normalization_sensitivity/normalization_sensitivity_null_primary108.csv",
        "results/normalization_sensitivity/normalization_sensitivity_metrics.json",
        "results/reviewer_minor_audits/shrinkage_k_sensitivity_summary.csv",
        "results/reviewer_minor_audits/shrinkage_k_sensitivity_drug124.csv",
        "results/reviewer_minor_audits/MIPE_missing16_audit.csv",
        "results/reviewer_minor_audits/shrinkage_missingness_metrics.json",
    }

    assert expected <= set(SCIENTIFIC_OUTPUTS)


def test_rev13_posthoc_outputs_are_in_reproduction_contract() -> None:
    expected = {
        "results/leave_one_seed_out/leave_one_seed_out_summary.csv",
        "results/leave_one_seed_out/leave_one_seed_out_seed_summary.csv",
        "results/leave_one_seed_out/leave_one_seed_out_metrics.json",
        "results/leave_one_seed_out/leave_one_seed_out_audit.md",
        "results/seed_excluded_scoring/seed_excluded_scores.csv",
        "results/seed_excluded_scoring/seed_excluded_metrics.json",
        "results/seed_excluded_scoring/seed_excluded_scoring_audit.md",
    }

    assert expected <= set(SCIENTIFIC_OUTPUTS)


def test_adjudicated_evidence_outputs_are_in_reproduction_contract() -> None:
    expected = {
        "data/evidence/evidence_labels_v3_adjudicated.csv",
        "results/evidence_audit/second_reviewer_agreement.json",
        "results/evidence_audit/second_reviewer_agreement.md",
        "results/evidence_audit/second_reviewer_disagreements.csv",
        "results/evidence_audit/second_reviewer_adjudication_final.csv",
    }

    assert expected <= set(SCIENTIFIC_OUTPUTS)
    assert len(REQUIRED_INPUTS) == 14
    assert len(SCIENTIFIC_OUTPUTS) == 38


def test_frozen_input_validation_detects_tampering(tmp_path: Path) -> None:
    input_path = tmp_path / "input.csv"
    input_path.write_text("a,b\n1,2\n", encoding="utf-8")
    manifest = {"input.csv": sha256_file(input_path)}

    validate_frozen_inputs(tmp_path, manifest)
    input_path.write_text("a,b\n1,3\n", encoding="utf-8")

    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        validate_frozen_inputs(tmp_path, manifest)


def test_reproduction_comparison_rejects_changed_scientific_output(
    tmp_path: Path,
) -> None:
    authoritative = tmp_path / "authoritative"
    reproduced = tmp_path / "reproduced"
    authoritative.mkdir()
    reproduced.mkdir()
    (authoritative / "table.csv").write_text("x\n1\n", encoding="utf-8")
    (reproduced / "table.csv").write_text("x\n2\n", encoding="utf-8")

    comparison = compare_reproduction(
        authoritative,
        reproduced,
        ("table.csv",),
    )

    assert comparison["all_match"] is False
    assert comparison["files"]["table.csv"]["match"] is False


def test_primary_pipeline_no_longer_emits_retired_auc() -> None:
    from analysis.acc_primary_pipeline import compute_primary_analysis, load_inputs

    result = compute_primary_analysis(load_inputs(PROJECT_ROOT))

    assert "legacy_mixed_evidence_auc" not in result.metrics
    assert "legacy_benchmark_warning" not in result.metrics


def test_relative_run_directory_is_anchored_to_project_root(tmp_path: Path) -> None:
    project_root = tmp_path / "project"
    project_root.mkdir()

    resolved = resolve_run_root(project_root, Path("repro_outputs/check"))

    assert resolved == (project_root / "repro_outputs/check").resolve()


def test_input_manifest_is_copied_into_isolated_run(tmp_path: Path) -> None:
    source = tmp_path / "source_manifest.json"
    source.write_text('{"sha256": {}}\n', encoding="utf-8")
    run_root = tmp_path / "run"

    copied = copy_input_manifest(source, run_root)

    assert copied == run_root / "reproducibility" / "input_manifest.json"
    assert copied.read_bytes() == source.read_bytes()
