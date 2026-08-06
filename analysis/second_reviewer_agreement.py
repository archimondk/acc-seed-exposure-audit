"""Compare evidence reviews and integrate the blinded B02 access completion."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.stats import beta


FIELD_MAP = {
    "evidence_domain": "evidence_domain",
    "evidence_design": "evidence_design",
    "drug_specificity": "drug_specificity",
    "exposure_context": "exposure_context",
    "direction_v2": "direction",
    "strict_candidate_eligible": "strict_candidate_eligible",
    "strict_binary_label": "strict_binary_label",
}

ADOPT_SECOND_DIRECTION = {
    "B11": "regimen_effect",
    "B14": "regimen_effect",
    "B15": "class_extrapolated",
    "B19": "regimen_effect",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def normalize(value: str) -> str:
    value = value.strip()
    return "NA" if value in {"", "NA"} else value


def apply_b02_rerate(
    locked_second: list[dict[str, str]],
    rerate_path: Path,
) -> list[dict[str, str]]:
    """Return a copy of the locked review with only B02 access-completed."""
    rerate = read_csv(rerate_path)
    if len(rerate) != 1 or rerate[0].get("record_id") != "B02":
        raise ValueError("B02 rerate must contain exactly one B02 record")
    if len(locked_second) < 2 or locked_second[1].get("record_id") != "B02":
        raise ValueError("Locked review must contain B02 once at index 1")
    if list(rerate[0]) != list(locked_second[1]):
        raise ValueError("B02 rerate columns do not match the locked review")

    for identity_field in ("record_id", "source_id", "drug", "doi_or_identifier"):
        if rerate[0][identity_field] != locked_second[1][identity_field]:
            raise ValueError(f"B02 rerate changed {identity_field}")

    required = {
        "source_accessed",
        "access_level",
        *FIELD_MAP.values(),
    }
    blank = sorted(
        field for field in required if not rerate[0].get(field, "").strip()
    )
    if blank:
        raise ValueError(f"B02 rerate has blank required fields: {blank}")
    if rerate[0]["source_accessed"] != "yes":
        raise ValueError("B02 rerate must document successful source access")

    merged = [dict(row) for row in locked_second]
    merged[1] = dict(rerate[0])
    return merged


def cohen_kappa(first: list[str], second: list[str]) -> float:
    if len(first) != len(second) or not first:
        raise ValueError("Kappa inputs must be non-empty and paired")
    n = len(first)
    observed = sum(a == b for a, b in zip(first, second)) / n
    first_counts = Counter(first)
    second_counts = Counter(second)
    expected = sum(
        first_counts[label] / n * second_counts[label] / n
        for label in set(first_counts) | set(second_counts)
    )
    if expected == 1:
        return 1.0 if observed == 1 else 0.0
    return (observed - expected) / (1 - expected)


def exact_binomial_interval(
    successes: int,
    n: int,
    confidence: float = 0.95,
) -> tuple[float, float]:
    if not 0 <= successes <= n or n < 1:
        raise ValueError("Binomial interval inputs are invalid")
    alpha = 1.0 - confidence
    lower = (
        0.0
        if successes == 0
        else float(beta.ppf(alpha / 2.0, successes, n - successes + 1))
    )
    upper = (
        1.0
        if successes == n
        else float(beta.ppf(1.0 - alpha / 2.0, successes + 1, n - successes))
    )
    return lower, upper


def paired_bootstrap_kappa(
    first: list[str],
    second: list[str],
    n_resamples: int = 10_000,
    rng_seed: int = 20260728,
) -> dict[str, object]:
    if len(first) != len(second) or not first:
        raise ValueError("Bootstrap inputs must be non-empty and paired")
    rng = np.random.default_rng(rng_seed)
    values: list[float] = []
    undefined = 0
    for _ in range(n_resamples):
        indices = rng.integers(0, len(first), size=len(first))
        first_sample = [first[int(index)] for index in indices]
        second_sample = [second[int(index)] for index in indices]
        if len(set(first_sample)) < 2 or len(set(second_sample)) < 2:
            undefined += 1
            continue
        values.append(cohen_kappa(first_sample, second_sample))
    interval = (
        [float(value) for value in np.percentile(values, [2.5, 97.5])]
        if values
        else [None, None]
    )
    return {
        "ci_95": interval,
        "resamples": n_resamples,
        "defined_n": len(values),
        "undefined_n": undefined,
        "undefined_reason": (
            "a bootstrap resample contained only one marginal category"
        ),
    }


def metrics_for(
    first: list[dict[str, str]],
    second: list[dict[str, str]],
    indices: Iterable[int],
) -> dict[str, dict]:
    indices = list(indices)
    output: dict[str, dict] = {}
    for first_field, second_field in FIELD_MAP.items():
        first_values = [normalize(first[i][first_field]) for i in indices]
        second_values = [normalize(second[i][second_field]) for i in indices]
        mismatches = [
            second[i]["record_id"]
            for i, a, b in zip(indices, first_values, second_values)
            if a != b
        ]
        agreement_count = len(indices) - len(mismatches)
        agreement_ci = exact_binomial_interval(
            agreement_count,
            len(indices),
        )
        kappa_bootstrap = paired_bootstrap_kappa(
            first_values,
            second_values,
        )
        output[first_field] = {
            "second_field": second_field,
            "n": len(indices),
            "agreement_count": agreement_count,
            "agreement_rate": agreement_count / len(indices),
            "agreement_exact_ci_95": list(agreement_ci),
            "cohen_kappa": cohen_kappa(first_values, second_values),
            "kappa_paired_bootstrap_ci_95": kappa_bootstrap["ci_95"],
            "kappa_bootstrap_resamples": kappa_bootstrap["resamples"],
            "kappa_bootstrap_defined_n": kappa_bootstrap["defined_n"],
            "kappa_bootstrap_undefined_n": kappa_bootstrap["undefined_n"],
            "kappa_bootstrap_undefined_reason": kappa_bootstrap[
                "undefined_reason"
            ],
            "mismatch_record_ids": mismatches,
        }
    return output


def recommendation(record_id: str, field: str) -> str:
    if record_id == "B08":
        return (
            "retain_first_taxonomy_class_only_or_class_extrapolated;"
            "_strict_consequence_unchanged"
        )
    if record_id in {"B11", "B14", "B19"} and field == "direction_v2":
        return "adopt_second_regimen_effect_for_combination_consistency"
    if record_id == "B15" and field == "direction_v2":
        return "adopt_second_class_extrapolated_for_taxonomy_consistency"
    return "manual_adjudication"


def collect_disagreements(
    first: list[dict[str, str]],
    second: list[dict[str, str]],
) -> list[dict[str, str]]:
    disagreements: list[dict[str, str]] = []
    for first_row, second_row in zip(first, second):
        for first_field, second_field in FIELD_MAP.items():
            first_value = normalize(first_row[first_field])
            second_value = normalize(second_row[second_field])
            if first_value == second_value:
                continue
            disagreements.append(
                {
                    "record_id": second_row["record_id"],
                    "drug": second_row["drug"],
                    "field": first_field,
                    "first_value": first_value,
                    "second_value": second_value,
                    "recommendation": recommendation(
                        second_row["record_id"], first_field
                    ),
                }
            )
    return disagreements


def build_adjudicated_table(
    first: list[dict[str, str]],
    second: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Create v3 without changing either frozen source review."""
    output = [dict(row) for row in first]
    for row, second_row in zip(output, second):
        record_id = second_row["record_id"]
        if record_id == "B02":
            row["independent_second_review"] = "agree_after_access_completion"
        elif record_id == "B08":
            row["independent_second_review"] = "adjudicated_retain_first"
        elif record_id in ADOPT_SECOND_DIRECTION:
            row["direction_v2"] = ADOPT_SECOND_DIRECTION[record_id]
            row["independent_second_review"] = "adjudicated_adopt_second"
        else:
            row["independent_second_review"] = "agree"
    return output


def adjudication_rows() -> list[dict[str, str]]:
    return [
        {
            "record_id": "B02",
            "drug": "Cabozantinib",
            "issue": "source_access_incomplete",
            "decision": "resolved_by_blinded_access_completion; retain_first_labels",
            "final_direction": "positive",
            "strict_effect": "yes/positive",
            "basis": "All seven compared fields agree after the official PubMed record was reviewed.",
        },
        {
            "record_id": "B08",
            "drug": "Carboplatin",
            "issue": "specificity_and_direction_taxonomy",
            "decision": "retain_first_taxonomy",
            "final_direction": "class_extrapolated",
            "strict_effect": "no/NA unchanged",
            "basis": "The protocol maps evidence transferred from a related platinum drug or class to class extrapolation.",
        },
        {
            "record_id": "B11",
            "drug": "Lenvatinib",
            "issue": "direction_taxonomy",
            "decision": "adopt_second_taxonomy",
            "final_direction": "regimen_effect",
            "strict_effect": "no/NA unchanged",
            "basis": "The exposure was a non-decomposable combination regimen.",
        },
        {
            "record_id": "B14",
            "drug": "Temsirolimus",
            "issue": "direction_taxonomy",
            "decision": "adopt_second_taxonomy",
            "final_direction": "regimen_effect",
            "strict_effect": "no/NA unchanged",
            "basis": "The exposure was a non-decomposable combination regimen.",
        },
        {
            "record_id": "B15",
            "drug": "Sirolimus",
            "issue": "direction_taxonomy",
            "decision": "adopt_second_taxonomy",
            "final_direction": "class_extrapolated",
            "strict_effect": "no/NA unchanged",
            "basis": "The source evaluated everolimus rather than sirolimus.",
        },
        {
            "record_id": "B19",
            "drug": "Erlotinib",
            "issue": "direction_taxonomy",
            "decision": "adopt_second_taxonomy",
            "final_direction": "regimen_effect",
            "strict_effect": "no/NA unchanged",
            "basis": "The exposure was a non-decomposable combination regimen.",
        },
    ]


def build_report(
    final_metrics: dict[str, dict],
    pre_rerate_metrics: dict[str, dict],
    disagreements: list[dict[str, str]],
) -> str:
    strict_final = final_metrics["strict_candidate_eligible"]
    binary_final = final_metrics["strict_binary_label"]
    strict_pre = pre_rerate_metrics["strict_candidate_eligible"]

    lines = [
        "# Final Claude–human traceability comparison and adjudication",
        "",
        (
            "This report compares a protocol-locked, blinded classification "
            "generated by Anthropic Claude with the primary human curator's "
            "labels. It is not an independent human rereview and must not be "
            "interpreted as human inter-rater reliability. The exact Claude "
            "model/version was not recorded and could not be recovered."
        ),
        "",
        "## Main finding",
        "",
        (
            "After the blinded B02 access-completion rerating, "
            f"strict-eligibility agreement is {strict_final['agreement_count']}/"
            f"{strict_final['n']} ({strict_final['agreement_rate']:.1%}), and "
            f"strict-binary-label agreement is {binary_final['agreement_count']}/"
            f"{binary_final['n']} ({binary_final['agreement_rate']:.1%})."
        ),
        (
            "For B02 cabozantinib, the model-based classification is `clinical`, "
            "`prospective_single_arm_phase_2`, `direct`, `monotherapy`, "
            "`positive`, and strict-eligible. Every compared B02 field agrees "
            "with the frozen first review."
        ),
        (
            "Precision is limited by n = 19 and the 2/17 strict-eligibility "
            "marginal split. For 19/19 raw agreement, the exact 95% binomial "
            f"interval is [{strict_final['agreement_exact_ci_95'][0]:.3f}, "
            f"{strict_final['agreement_exact_ci_95'][1]:.3f}]. This interval "
            "describes concordance in this fixed audit set; it is not a human "
            "inter-rater reliability estimate."
        ),
        (
            f"Before source completion, strict agreement was "
            f"{strict_pre['agreement_count']}/{strict_pre['n']} "
            f"({strict_pre['agreement_rate']:.1%}); that access-related interim "
            "result is retained in the `*_pre_rerate` audit files."
        ),
        "",
        "## Field-level agreement",
        "",
        "| Field | Final agreement | Pre-rerate agreement |",
        "|---|---:|---:|",
    ]
    for field in FIELD_MAP:
        final_item = final_metrics[field]
        pre_item = pre_rerate_metrics[field]
        lines.append(
            f"| {field} | {final_item['agreement_count']}/{final_item['n']} "
            f"({final_item['agreement_rate']:.1%}) | "
            f"{pre_item['agreement_count']}/{pre_item['n']} "
            f"({pre_item['agreement_rate']:.1%}) |"
        )

    lines.extend(
        [
            "",
            "## Remaining categorical disagreements",
            "",
            "| Record | Drug | Field | First review | Second review | Adjudication |",
            "|---|---|---|---|---|---|",
        ]
    )
    for item in disagreements:
        lines.append(
            f"| {item['record_id']} | {item['drug']} | {item['field']} | "
            f"{item['first_value']} | {item['second_value']} | "
            f"{item['recommendation']} |"
        )

    lines.extend(
        [
            "",
            "## Final adjudication",
            "",
            "- **B02 Cabozantinib:** resolved by source-access completion. "
            "Retain the first-review labels; the strict result is `yes/positive`.",
            "- **B08 Carboplatin:** retain `class_only`/`class_extrapolated`; "
            "the protocol maps related-drug or class transfer to class "
            "extrapolation. Strict exclusion is unchanged.",
            "- **B11 Lenvatinib, B14 Temsirolimus and B19 Erlotinib:** adopt "
            "`regimen_effect` because exposure occurred in non-decomposable "
            "combination regimens. Strict exclusions are unchanged.",
            "- **B15 Sirolimus:** adopt `class_extrapolated` because the source "
            "tested everolimus rather than sirolimus. Strict exclusion is unchanged.",
            "- **B06 Gemcitabine and B16 Sunitinib:** no categorical change; "
            "the model notes are retained as interpretive cautions.",
            "",
            "## Boundary",
            "",
            "The frozen primary human review (`evidence_labels_v2.csv`) and locked "
            "Claude model output remain unchanged. Final taxonomy decisions are "
            "written to `evidence_labels_v3_adjudicated.csv`; downstream "
            "manuscript and figure synchronization should use v3 only after "
            "the reproducibility inputs are deliberately version-bumped.",
        ]
    )
    return "\n".join(lines) + "\n"


def analyze(project_root: Path) -> dict:
    first_path = project_root / "data/evidence/evidence_labels_v2.csv"
    locked_second_path = (
        project_root
        / "independent_review/claude_blinded_evidence_review_v1"
        / "claude_review_locked.csv"
    )
    b02_rerate_path = (
        project_root
        / "independent_review/cabozantinib_b02_rerate"
        / "claude_b02_rerated.csv"
    )
    output_dir = project_root / "results/evidence_audit"
    output_dir.mkdir(parents=True, exist_ok=True)

    first = read_csv(first_path)
    locked_second = read_csv(locked_second_path)
    second = apply_b02_rerate(locked_second, b02_rerate_path)
    if len(first) != 19 or len(locked_second) != 19:
        raise ValueError("Both frozen reviews must contain exactly 19 records")
    if [row["drug"] for row in first] != [row["drug"] for row in second]:
        raise ValueError("First- and second-review drug order does not match")
    expected_ids = [f"B{i:02d}" for i in range(1, 20)]
    if [row["record_id"] for row in second] != expected_ids:
        raise ValueError("Second-review record IDs are missing or reordered")

    all_indices = list(range(19))
    pre_rerate_metrics = metrics_for(first, locked_second, all_indices)
    final_metrics = metrics_for(first, second, all_indices)
    disagreements = collect_disagreements(first, second)

    summary = {
        "review_version": "blinded-second-review-agreement-v2-access-completed",
        "first_review": str(first_path.relative_to(project_root)),
        "locked_second_review": str(
            locked_second_path.relative_to(project_root)
        ),
        "b02_access_completion": str(
            b02_rerate_path.relative_to(project_root)
        ),
        "all_records": final_metrics,
        "locked_pre_rerate_all_records": pre_rerate_metrics,
        "evaluable_records_excluding_access_none": final_metrics,
        "evaluable_record_ids": [row["record_id"] for row in second],
        "excluded_for_access": [],
        "categorical_disagreement_count": len(disagreements),
        "disagreement_record_ids": sorted(
            {item["record_id"] for item in disagreements}
        ),
    }
    json_path = output_dir / "second_reviewer_agreement.json"
    json_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    disagreement_path = output_dir / "second_reviewer_disagreements.csv"
    with disagreement_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=(
                "record_id",
                "drug",
                "field",
                "first_value",
                "second_value",
                "recommendation",
            ),
        )
        writer.writeheader()
        writer.writerows(disagreements)

    adjudication_path = (
        output_dir / "second_reviewer_adjudication_final.csv"
    )
    adjudications = adjudication_rows()
    with adjudication_path.open(
        "w", encoding="utf-8", newline=""
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=adjudications[0].keys())
        writer.writeheader()
        writer.writerows(adjudications)

    adjudicated_path = (
        project_root / "data/evidence/evidence_labels_v3_adjudicated.csv"
    )
    adjudicated = build_adjudicated_table(first, second)
    with adjudicated_path.open(
        "w", encoding="utf-8", newline=""
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=adjudicated[0].keys())
        writer.writeheader()
        writer.writerows(adjudicated)

    report_path = output_dir / "second_reviewer_agreement.md"
    report_path.write_text(
        build_report(final_metrics, pre_rerate_metrics, disagreements),
        encoding="utf-8",
    )
    return {
        "summary": summary,
        "json": str(json_path),
        "disagreements": str(disagreement_path),
        "adjudication": str(adjudication_path),
        "adjudicated_evidence": str(adjudicated_path),
        "report": str(report_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    args = parser.parse_args()
    result = analyze(args.project_root.resolve())
    print(
        json.dumps(
            {
                "status": "ok",
                "disagreement_records": result["summary"][
                    "disagreement_record_ids"
                ],
                "excluded_for_access": result["summary"][
                    "excluded_for_access"
                ],
                "report": result["report"],
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
