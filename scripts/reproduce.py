"""Run the frozen core workflow plus rev13 Amendments 4-5 in isolation."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib
import numpy as np
import scipy
from PIL import __version__ as pillow_version

from analysis.acc_primary_pipeline import run_pipeline
from analysis.evidence_label_audit import generate_c4_outputs
from analysis.leave_one_seed_out import run as run_leave_one_seed_out
from analysis.mechanism_enrichment import run as run_mechanism_enrichment
from analysis.method_strengthening import run as run_method_strengthening
from analysis.method_strengthening_figure import (
    generate_method_strengthening_figure,
)
from analysis.normalization_sensitivity import run as run_normalization_sensitivity
from analysis.revision_figures import generate_c3_figures
from analysis.second_reviewer_agreement import analyze as analyze_second_review
from analysis.seed_excluded_scoring import run as run_seed_excluded_scoring
from analysis.shrinkage_missingness_audit import (
    run as run_shrinkage_missingness_audit,
)


RUNNER_VERSION = "rev13-reproduce-v5-amendments-4-5"

REQUIRED_INPUTS = (
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
)

SCIENTIFIC_OUTPUTS = (
    "results/primary_analysis/ADRS_comp_primary_108.csv",
    "results/primary_analysis/ADRS_evidence_informed_108.csv",
    "results/primary_analysis/ADRS_context_only_16.csv",
    "results/primary_analysis/primary_metrics.json",
    "results/mechanism_enrichment/mechanism_enrichment_primary108.csv",
    "results/mechanism_enrichment/mechanism_members_primary108.csv",
    "figure_data/revision/Fig3_component_correlation_primary108.csv",
    "figure_data/revision/Fig5a_weight_scan_primary108.csv",
    "figure_data/revision/Fig5b_CDK46_exact_null_primary108.csv",
    "figure_data/revision/C3_figure_stats.json",
    "figure_data/revision/Fig4_evidence_audit_primary108.csv",
    "results/evidence_audit/evidence_audit_metrics.json",
    "data/evidence/evidence_labels_v3_adjudicated.csv",
    "results/evidence_audit/second_reviewer_agreement.json",
    "results/evidence_audit/second_reviewer_agreement.md",
    "results/evidence_audit/second_reviewer_disagreements.csv",
    "results/evidence_audit/second_reviewer_adjudication_final.csv",
    "results/method_strengthening/baseline_comparison_primary108.csv",
    "results/method_strengthening/centrality_gene399.csv",
    "results/method_strengthening/centrality_drug108.csv",
    "results/method_strengthening/degree_matched_seed_sets.csv",
    "results/method_strengthening/random_seed_null_primary108.csv",
    "results/method_strengthening/method_strengthening_metrics.json",
    "figure_data/revision/Fig2d_random_seed_top12_primary108.csv",
    "results/normalization_sensitivity/normalization_sensitivity_summary.csv",
    "results/normalization_sensitivity/normalization_sensitivity_null_primary108.csv",
    "results/normalization_sensitivity/normalization_sensitivity_metrics.json",
    "results/reviewer_minor_audits/shrinkage_k_sensitivity_summary.csv",
    "results/reviewer_minor_audits/shrinkage_k_sensitivity_drug124.csv",
    "results/reviewer_minor_audits/MIPE_missing16_audit.csv",
    "results/reviewer_minor_audits/shrinkage_missingness_metrics.json",
    "results/leave_one_seed_out/leave_one_seed_out_summary.csv",
    "results/leave_one_seed_out/leave_one_seed_out_seed_summary.csv",
    "results/leave_one_seed_out/leave_one_seed_out_metrics.json",
    "results/leave_one_seed_out/leave_one_seed_out_audit.md",
    "results/seed_excluded_scoring/seed_excluded_scores.csv",
    "results/seed_excluded_scoring/seed_excluded_metrics.json",
    "results/seed_excluded_scoring/seed_excluded_scoring_audit.md",
)

EXPECTED_FIGURES = (
    "figures/revision/Fig2_method_strengthening_primary108.pdf",
    "figures/revision/Fig2_method_strengthening_primary108.svg",
    "figures/revision/Fig2_method_strengthening_primary108.png",
    "figures/revision/Fig3_weight_stability_CDK46_primary108.pdf",
    "figures/revision/Fig3_weight_stability_CDK46_primary108.svg",
    "figures/revision/Fig3_weight_stability_CDK46_primary108.png",
    "figures/revision/Fig4_evidence_audit_primary108.pdf",
    "figures/revision/Fig4_evidence_audit_primary108.svg",
    "figures/revision/Fig4_evidence_audit_primary108.png",
)

CODE_FILES = (
    "analysis/acc_primary_pipeline.py",
    "analysis/mechanism_enrichment.py",
    "analysis/revision_figures.py",
    "analysis/evidence_label_audit.py",
    "analysis/second_reviewer_agreement.py",
    "analysis/method_strengthening.py",
    "analysis/method_strengthening_figure.py",
    "analysis/normalization_sensitivity.py",
    "analysis/shrinkage_missingness_audit.py",
    "analysis/leave_one_seed_out.py",
    "analysis/seed_excluded_scoring.py",
    "scripts/reproduce.py",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_run_root(project_root: Path, requested: Path) -> Path:
    """Resolve a user-supplied run directory relative to the project root."""
    project_root = project_root.resolve()
    return (
        requested.resolve()
        if requested.is_absolute()
        else (project_root / requested).resolve()
    )


def load_frozen_input_manifest(path: Path) -> dict[str, str]:
    if not path.is_file():
        raise FileNotFoundError(f"Frozen input manifest is missing: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    hashes = payload.get("sha256")
    if not isinstance(hashes, dict):
        raise ValueError("Input manifest must contain an object named 'sha256'")
    normalized = {str(key): str(value) for key, value in hashes.items()}
    if set(normalized) != set(REQUIRED_INPUTS):
        missing = sorted(set(REQUIRED_INPUTS) - set(normalized))
        extra = sorted(set(normalized) - set(REQUIRED_INPUTS))
        raise ValueError(
            f"Input manifest inventory mismatch; missing={missing}, extra={extra}"
        )
    return normalized


def validate_frozen_inputs(
    project_root: Path,
    expected_hashes: Mapping[str, str],
) -> dict[str, str]:
    actual: dict[str, str] = {}
    for relative in expected_hashes:
        path = project_root / relative
        if not path.is_file():
            raise FileNotFoundError(f"Required frozen input is missing: {path}")
        digest = sha256_file(path)
        actual[relative] = digest
        if digest != expected_hashes[relative]:
            raise ValueError(
                f"SHA-256 mismatch for {relative}: "
                f"expected {expected_hashes[relative]}, got {digest}"
            )
    return actual


def _assert_empty_destination(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(
            f"Run directory must be absent or empty; refusing to overwrite: {path}"
        )
    path.mkdir(parents=True, exist_ok=True)


def materialize_inputs(
    project_root: Path,
    run_root: Path,
) -> dict[str, str]:
    copied: dict[str, str] = {}
    for relative in REQUIRED_INPUTS:
        source = project_root / relative
        destination = run_root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        copied[relative] = sha256_file(destination)
    return copied


def copy_input_manifest(manifest_path: Path, run_root: Path) -> Path:
    """Store the exact frozen-input inventory beside the isolated run."""
    destination = run_root / "reproducibility" / "input_manifest.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(manifest_path, destination)
    return destination


def compare_reproduction(
    authoritative_root: Path,
    reproduced_root: Path,
    relative_paths: Sequence[str] = SCIENTIFIC_OUTPUTS,
) -> dict[str, Any]:
    files: dict[str, dict[str, Any]] = {}
    for relative in relative_paths:
        authoritative = authoritative_root / relative
        reproduced = reproduced_root / relative
        authoritative_hash = (
            sha256_file(authoritative) if authoritative.is_file() else None
        )
        reproduced_hash = sha256_file(reproduced) if reproduced.is_file() else None
        files[relative] = {
            "authoritative_exists": authoritative.is_file(),
            "reproduced_exists": reproduced.is_file(),
            "authoritative_sha256": authoritative_hash,
            "reproduced_sha256": reproduced_hash,
            "match": (
                authoritative_hash is not None
                and reproduced_hash is not None
                and authoritative_hash == reproduced_hash
            ),
        }
    return {
        "all_match": all(item["match"] for item in files.values()),
        "files": files,
    }


def _environment() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor() or "not reported by operating system",
        "logical_cpu_cores": str(os.cpu_count() or "not reported"),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "pillow": pillow_version,
    }


def _collect_output_hashes(run_root: Path) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for folder in ("results", "figure_data", "figures", "projects", "databases"):
        root = run_root / folder
        if not root.is_dir():
            continue
        for path in sorted(item for item in root.rglob("*") if item.is_file()):
            relative = path.relative_to(run_root).as_posix()
            if relative.startswith("results/reproducibility/"):
                continue
            hashes[relative] = sha256_file(path)
    return hashes


def _write_run_records(
    project_root: Path,
    run_root: Path,
    payload: Mapping[str, Any],
) -> tuple[Path, Path, Path]:
    output_dir = run_root / "results" / "reproducibility"
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "run_manifest.json"
    markdown_path = output_dir / "run_manifest.md"
    patches_path = output_dir / "PATCHES.md"
    json_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    comparison = payload["comparison"]
    mismatch_lines = [
        f"- `{path}`"
        for path, item in comparison["files"].items()
        if not item["match"]
    ]
    run_display = (
        run_root.relative_to(project_root).as_posix()
        if run_root.is_relative_to(project_root)
        else str(run_root)
    )
    markdown_path.write_text(
        "\n".join(
            [
                "# Method-strengthening reproducibility run manifest",
                "",
                f"- Status: **{payload['status']}**",
                f"- Runner: `{payload['runner_version']}`",
                f"- Run directory: `{run_display}`",
                f"- UTC completed: `{payload['completed_at_utc']}`",
                f"- Wall-clock runtime: `{payload['wall_clock_seconds']:.1f} s`",
                "- Randomness: C2 rank-set inference is exact; the degree-matched "
                "null uses RNG seed 20260727 with 10,000 draws.",
                "",
                "## Environment",
                "",
                *[
                    f"- {name}: `{value}`"
                    for name, value in payload["environment"].items()
                ],
                "",
                "## Gates",
                "",
                f"- Frozen input hashes: {'pass' if payload['input_hash_gate'] else 'fail'}",
                f"- Scientific output equivalence: "
                f"{'pass' if comparison['all_match'] else 'fail'}",
                f"- Expected figure files: "
                f"{'pass' if payload['figure_gate']['all_present'] else 'fail'}",
                f"- C4 evidence gate: "
                f"{'pass' if payload['c4_gate']['pass'] else 'fail'}",
                f"- Method-strengthening gate: "
                f"{'pass' if payload['method_strengthening_gate']['pass'] else 'fail'}",
                f"- Normalization-sensitivity gate: "
                f"{'pass' if payload['normalization_sensitivity_gate']['pass'] else 'fail'}",
                f"- Reviewer minor-audit gate: "
                f"{'pass' if payload['reviewer_minor_audit_gate']['pass'] else 'fail'}",
                f"- Rev13 Amendments 4-5 gate: "
                f"{'pass' if payload['rev13_posthoc_gate']['pass'] else 'fail'}",
                "",
                "## Scientific output mismatches",
                "",
                *(mismatch_lines or ["- None."]),
                "",
                "## Portable command",
                "",
                "```text",
                "python -m scripts.reproduce --project-root .",
                "```",
                "",
            ]
        ),
        encoding="utf-8",
    )
    patches_path.write_text(
        "\n".join(
            [
                "# Reproduction patches",
                "",
                "No source or frozen input file was modified during this run. Required "
                "inputs were copied into the isolated run directory, and all generated "
                "artifacts were written below that directory.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return json_path, markdown_path, patches_path


def run_reproduction(
    project_root: Path,
    run_root: Path,
    input_manifest_path: Path | None = None,
) -> dict[str, Any]:
    started_at = time.perf_counter()
    project_root = project_root.resolve()
    run_root = run_root.resolve()
    if run_root == project_root:
        raise ValueError("Run directory must not be the project root")
    if input_manifest_path is None:
        manifest_path = project_root / "reproducibility" / "input_manifest.json"
    else:
        manifest_path = (
            input_manifest_path.resolve()
            if input_manifest_path.is_absolute()
            else (project_root / input_manifest_path).resolve()
        )
    expected_hashes = load_frozen_input_manifest(manifest_path)
    validated_hashes = validate_frozen_inputs(project_root, expected_hashes)
    _assert_empty_destination(run_root)
    copied_manifest_path = copy_input_manifest(manifest_path, run_root)
    copied_hashes = materialize_inputs(project_root, run_root)
    if copied_hashes != validated_hashes:
        raise ValueError("Copied input hashes differ from validated source inputs")

    run_pipeline(run_root)
    run_mechanism_enrichment(run_root)
    generate_c3_figures(run_root)
    analyze_second_review(run_root)
    generate_c4_outputs(run_root)
    run_method_strengthening(run_root)
    generate_method_strengthening_figure(run_root)
    run_normalization_sensitivity(run_root)
    run_shrinkage_missingness_audit(run_root)
    run_leave_one_seed_out(run_root)
    run_seed_excluded_scoring(run_root)

    comparison = compare_reproduction(project_root, run_root)
    figure_status = {
        relative: (run_root / relative).is_file() for relative in EXPECTED_FIGURES
    }
    c4_metrics_path = (
        run_root / "results" / "evidence_audit" / "evidence_audit_metrics.json"
    )
    c4_metrics = json.loads(c4_metrics_path.read_text(encoding="utf-8"))
    c4_gate = {
        "auc_estimable_is_false": c4_metrics.get("auc_estimable") is False,
        "strict_positive_n": c4_metrics.get("strict_positive_n"),
        "strict_negative_n": c4_metrics.get("strict_negative_n"),
        "second_reviewer_status": c4_metrics.get("second_reviewer_status"),
        "strict_binary_agreement_n": c4_metrics.get(
            "strict_binary_agreement_n"
        ),
        "strict_binary_agreement_total": c4_metrics.get(
            "strict_binary_agreement_total"
        ),
        "strict_binary_cohen_kappa": c4_metrics.get(
            "strict_binary_cohen_kappa"
        ),
    }
    c4_gate["pass"] = (
        c4_gate["auc_estimable_is_false"]
        and c4_gate["strict_positive_n"] == 2
        and c4_gate["strict_negative_n"] == 0
        and c4_gate["second_reviewer_status"]
        == "completed_and_adjudicated"
        and c4_gate["strict_binary_agreement_n"] == 19
        and c4_gate["strict_binary_agreement_total"] == 19
        and c4_gate["strict_binary_cohen_kappa"] == 1.0
    )
    method_metrics_path = (
        run_root
        / "results"
        / "method_strengthening"
        / "method_strengthening_metrics.json"
    )
    method_metrics = json.loads(method_metrics_path.read_text(encoding="utf-8"))
    method_decisions = method_metrics.get("hypothesis_decisions", {})
    method_strengthening_gate = {
        "null_draws": method_metrics.get("null_draws"),
        "H1": method_decisions.get("H1_nonredundancy", {}).get("status"),
        "BH_resolution_adequate": method_metrics.get(
            "BH_resolution_adequate"
        ),
        "H2": method_decisions.get(
            "H2_disease_context_beyond_centrality", {}
        ).get("status"),
        "H3": method_decisions.get("H3_CDK46_robustness", {}).get("status"),
    }
    method_strengthening_gate["pass"] = (
        method_strengthening_gate["null_draws"] == 10_000
        and method_strengthening_gate["BH_resolution_adequate"] is True
        and method_strengthening_gate["H1"] == "retired_descriptive_only"
        and method_strengthening_gate["H2"] == "partially_supported"
        and method_strengthening_gate["H3"] == "not_supported"
    )
    normalization_metrics = json.loads(
        (
            run_root
            / "results"
            / "normalization_sensitivity"
            / "normalization_sensitivity_metrics.json"
        ).read_text(encoding="utf-8")
    )
    normalization_variants = normalization_metrics.get("variants", {})
    normalization_gate = {
        "null_draws": normalization_metrics.get("null_draws"),
        "BH_resolution_adequate": normalization_metrics.get(
            "BH_resolution_adequate"
        ),
        "variant_names": sorted(normalization_variants),
        "column_minmax_drugs_q_lt_0_05": normalization_variants.get(
            "column_minmax", {}
        ).get("n_drugs_q_lt_0_05"),
        "uniform_ratio_gene_rho_degree": normalization_variants.get(
            "uniform_ratio_gene_rank", {}
        ).get("gene_rho_degree"),
        "symmetric_CDK46_q_across_variants": normalization_variants.get(
            "symmetric_gene_rank", {}
        ).get("CDK46_q_bh_across_variants"),
    }
    normalization_gate["pass"] = (
        normalization_gate["null_draws"] == 10_000
        and normalization_gate["BH_resolution_adequate"] is True
        and normalization_gate["variant_names"]
        == sorted(
            [
                "column_minmax",
                "column_gene_rank",
                "uniform_ratio_gene_rank",
                "symmetric_gene_rank",
            ]
        )
        and normalization_gate["column_minmax_drugs_q_lt_0_05"] == 2
        and normalization_gate["uniform_ratio_gene_rho_degree"] < 0.70
        and normalization_gate["symmetric_CDK46_q_across_variants"] < 0.05
    )
    minor_metrics = json.loads(
        (
            run_root
            / "results"
            / "reviewer_minor_audits"
            / "shrinkage_missingness_metrics.json"
        ).read_text(encoding="utf-8")
    )
    minor_audit_gate = {
        "pseudo_counts": minor_metrics.get("pseudo_counts"),
        "MIPE_missing_n": minor_metrics.get("MIPE_missing_n"),
    }
    minor_audit_gate["pass"] = (
        minor_audit_gate["pseudo_counts"] == [1.0, 3.0, 5.0, 10.0]
        and minor_audit_gate["MIPE_missing_n"] == 16
    )
    leave_one_seed_out_metrics = json.loads(
        (
            run_root
            / "results"
            / "leave_one_seed_out"
            / "leave_one_seed_out_metrics.json"
        ).read_text(encoding="utf-8")
    )
    seed_excluded_metrics = json.loads(
        (
            run_root
            / "results"
            / "seed_excluded_scoring"
            / "seed_excluded_metrics.json"
        ).read_text(encoding="utf-8")
    )
    rev13_posthoc_gate = {
        "leave_one_seed_out_run_count": leave_one_seed_out_metrics.get(
            "run_count"
        ),
        "leave_one_seed_out_minimum_rho": leave_one_seed_out_metrics.get(
            "minimum_rho"
        ),
        "leave_one_seed_out_rb1_max_abs_shift": leave_one_seed_out_metrics.get(
            "rb1_max_abs_shift"
        ),
        "seed_excluded_drug_count": seed_excluded_metrics.get("drug_count"),
        "seed_excluded_composite_spearman": seed_excluded_metrics.get(
            "composite_spearman"
        ),
        "seed_excluded_top20_intersection": seed_excluded_metrics.get(
            "composite_top20_intersection"
        ),
    }
    rev13_posthoc_gate["pass"] = (
        rev13_posthoc_gate["leave_one_seed_out_run_count"] == 180
        and abs(
            rev13_posthoc_gate["leave_one_seed_out_minimum_rho"]
            - 0.9167357359932169
        )
        < 1e-12
        and rev13_posthoc_gate["leave_one_seed_out_rb1_max_abs_shift"] == 57
        and rev13_posthoc_gate["seed_excluded_drug_count"] == 108
        and abs(
            rev13_posthoc_gate["seed_excluded_composite_spearman"]
            - 0.6459001255710277
        )
        < 1e-12
        and rev13_posthoc_gate["seed_excluded_top20_intersection"] == 7
    )
    code_hashes = {
        relative: sha256_file(project_root / relative) for relative in CODE_FILES
    }
    output_hashes = _collect_output_hashes(run_root)
    completed_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    status = (
        "verified"
        if comparison["all_match"]
        and all(figure_status.values())
        and c4_gate["pass"]
        and method_strengthening_gate["pass"]
        and normalization_gate["pass"]
        and minor_audit_gate["pass"]
        and rev13_posthoc_gate["pass"]
        else "failed"
    )
    payload: dict[str, Any] = {
        "status": status,
        "runner_version": RUNNER_VERSION,
        "completed_at_utc": completed_at,
        "wall_clock_seconds": time.perf_counter() - started_at,
        "environment": _environment(),
        "input_manifest": copied_manifest_path.relative_to(run_root).as_posix(),
        "input_manifest_sha256": sha256_file(copied_manifest_path),
        "input_hash_gate": copied_hashes == validated_hashes,
        "input_sha256": validated_hashes,
        "code_sha256": code_hashes,
        "output_sha256": output_hashes,
        "comparison": comparison,
        "figure_gate": {
            "all_present": all(figure_status.values()),
            "files": figure_status,
        },
        "c4_gate": c4_gate,
        "method_strengthening_gate": method_strengthening_gate,
        "normalization_sensitivity_gate": normalization_gate,
        "reviewer_minor_audit_gate": minor_audit_gate,
        "rev13_posthoc_gate": rev13_posthoc_gate,
    }
    records = _write_run_records(project_root, run_root, payload)
    payload["record_paths"] = [str(path) for path in records]
    if status != "verified":
        raise RuntimeError(
            f"Reproduction failed; inspect {records[1]} for gate results"
        )
    return payload


def _default_run_root(project_root: Path) -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return project_root / "repro_outputs" / f"run_{stamp}"


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        default=None,
        help="New or empty output directory. Defaults to repro_outputs/run_<UTC>.",
    )
    parser.add_argument(
        "--input-manifest",
        type=Path,
        default=None,
        help="Frozen input hash manifest; defaults to reproducibility/input_manifest.json.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    project_root = args.project_root.resolve()
    run_root = (
        resolve_run_root(project_root, args.run_dir)
        if args.run_dir is not None
        else _default_run_root(project_root)
    )
    result = run_reproduction(project_root, run_root, args.input_manifest)
    print(
        json.dumps(
            {
                "status": result["status"],
                "runner_version": result["runner_version"],
                "run_dir": str(run_root),
                "scientific_outputs_match": result["comparison"]["all_match"],
                "figures_present": result["figure_gate"]["all_present"],
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
