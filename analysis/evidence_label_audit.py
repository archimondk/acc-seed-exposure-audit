"""Audit heterogeneous ACC evidence labels and replace the invalid clinical AUC."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle


FIGURE_VERSION = "c4-evidence-audit-v3-agreement-precision"
TARGET_WIDTH_MM = 170.0
MIN_FONT_PT = 8.0
PNG_DPI = 1000

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
GREY = "#A7A9AC"
DARK = "#202124"
LIGHT_GREY = "#E5E7EB"

CATEGORY_ORDER = (
    "Drug-specific clinical monotherapy",
    "Clinical regimen/case (confounded)",
    "Direct ACC preclinical",
    "Class extrapolation",
)
CATEGORY_COLORS = {
    "Drug-specific clinical monotherapy": GREEN,
    "Clinical regimen/case (confounded)": ORANGE,
    "Direct ACC preclinical": BLUE,
    "Class extrapolation": GREY,
}


def _read_csv(path: Path, required: Iterable[str]) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Required input is missing: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"Input has no data rows: {path}")
    missing = set(required) - set(rows[0])
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
    return rows


def _parse_yes_no(value: str, field: str, drug: str) -> bool:
    normalized = value.strip().lower()
    if normalized not in {"yes", "no"}:
        raise ValueError(f"{drug}: {field} must be yes/no, got {value!r}")
    return normalized == "yes"


def load_evidence_labels(path: Path) -> list[dict[str, Any]]:
    required = (
        "drug",
        "legacy_benchmark_member",
        "legacy_direction",
        "legacy_score",
        "evidence_domain",
        "evidence_design",
        "drug_specificity",
        "exposure_context",
        "direction_v2",
        "strict_binary_label",
        "strict_candidate_eligible",
        "strict_exclusion_reason",
        "source_title",
        "source_year",
        "doi_or_identifier",
        "source_url",
        "locator_or_support",
        "source_verification",
        "independent_second_review",
        "used_in_external_score",
        "eligible_as_independent_validation",
    )
    raw = _read_csv(path, required)
    if len({row["drug"] for row in raw}) != len(raw):
        raise ValueError("Evidence table contains duplicate drug names")

    allowed_domains = {"clinical", "preclinical", "class_extrapolation"}
    allowed_legacy = {"positive", "negative", "neutral"}
    allowed_strict = {"", "positive", "negative"}
    rows: list[dict[str, Any]] = []
    for row in raw:
        drug = row["drug"].strip()
        parsed: dict[str, Any] = dict(row)
        parsed["drug"] = drug
        for field in (
            "legacy_benchmark_member",
            "strict_candidate_eligible",
            "used_in_external_score",
            "eligible_as_independent_validation",
        ):
            parsed[field] = _parse_yes_no(row[field], field, drug)
        parsed["legacy_score"] = float(row["legacy_score"])
        parsed["source_year"] = int(row["source_year"])
        parsed["strict_binary_label"] = row["strict_binary_label"].strip()

        if row["evidence_domain"] not in allowed_domains:
            raise ValueError(f"{drug}: invalid evidence_domain")
        if row["legacy_direction"] not in allowed_legacy:
            raise ValueError(f"{drug}: invalid legacy_direction")
        if parsed["strict_binary_label"] not in allowed_strict:
            raise ValueError(f"{drug}: invalid strict_binary_label")
        if not 0 <= parsed["legacy_score"] <= 1:
            raise ValueError(f"{drug}: legacy score is outside [0,1]")
        if not row["source_url"].startswith("https://"):
            raise ValueError(f"{drug}: source_url is not HTTPS")
        if parsed["strict_candidate_eligible"] != bool(
            parsed["strict_binary_label"]
        ):
            raise ValueError(
                f"{drug}: strict eligibility and strict label are inconsistent"
            )
        if (
            parsed["used_in_external_score"]
            and parsed["eligible_as_independent_validation"]
        ):
            raise ValueError(
                f"{drug}: reused external evidence cannot be independent validation"
            )
        rows.append(parsed)
    return rows


def load_primary_scores(path: Path) -> dict[str, dict[str, Any]]:
    raw = _read_csv(path, ("rank_comp", "drug", "ADRS_comp"))
    scores: dict[str, dict[str, Any]] = {}
    for row in raw:
        drug = row["drug"].strip()
        scores[drug] = {
            "rank_comp": int(row["rank_comp"]),
            "ADRS_comp": float(row["ADRS_comp"]),
        }
    if len(scores) != 108:
        raise ValueError(f"Expected frozen 108-drug universe, got {len(scores)}")
    if sorted(row["rank_comp"] for row in scores.values()) != list(range(1, 109)):
        raise ValueError("Primary ranks must be the consecutive integers 1..108")
    return scores


def load_second_review_summary(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Second-review agreement is missing: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("excluded_for_access") != []:
        raise ValueError("Second review still contains access-excluded records")
    for field in ("strict_candidate_eligible", "strict_binary_label"):
        metric = payload.get("all_records", {}).get(field, {})
        if (
            metric.get("n") != 19
            or metric.get("agreement_count") != 19
            or metric.get("cohen_kappa") != 1.0
        ):
            raise ValueError(
                f"Final second-review agreement gate failed for {field}"
            )
    return payload


def evidence_category(row: Mapping[str, Any]) -> str:
    if row["evidence_domain"] == "class_extrapolation":
        return "Class extrapolation"
    if row["evidence_domain"] == "preclinical":
        return "Direct ACC preclinical"
    if row["exposure_context"] in {"monotherapy", "monotherapy_or_standard"}:
        return "Drug-specific clinical monotherapy"
    return "Clinical regimen/case (confounded)"


def build_audit_rows(
    evidence_rows: Sequence[Mapping[str, Any]],
    primary_scores: Mapping[str, Mapping[str, Any]],
    *,
    legacy_only: bool,
) -> list[dict[str, Any]]:
    selected = [
        row
        for row in evidence_rows
        if not legacy_only or bool(row["legacy_benchmark_member"])
    ]
    missing = sorted(
        str(row["drug"]) for row in selected if row["drug"] not in primary_scores
    )
    if missing:
        raise ValueError(f"Evidence drugs missing from primary ranking: {missing}")
    result: list[dict[str, Any]] = []
    for row in selected:
        score = primary_scores[str(row["drug"])]
        result.append(
            {
                "drug": row["drug"],
                "rank_comp": int(score["rank_comp"]),
                "ADRS_comp": float(score["ADRS_comp"]),
                "legacy_direction": row["legacy_direction"],
                "legacy_score": float(row["legacy_score"]),
                "evidence_category": evidence_category(row),
                "direction_v2": row["direction_v2"],
                "strict_binary_label": row["strict_binary_label"],
                "strict_candidate_eligible": bool(
                    row["strict_candidate_eligible"]
                ),
                "strict_exclusion_reason": row["strict_exclusion_reason"],
                "doi_or_identifier": row["doi_or_identifier"],
                "source_url": row["source_url"],
                "source_verification": row["source_verification"],
                "independent_second_review": row["independent_second_review"],
                "eligible_as_independent_validation": bool(
                    row["eligible_as_independent_validation"]
                ),
            }
        )
    return sorted(result, key=lambda item: item["rank_comp"])


def evaluate_benchmark(
    evidence_rows: Sequence[Mapping[str, Any]],
    second_review: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    legacy = [row for row in evidence_rows if row["legacy_benchmark_member"]]
    strict = [row for row in legacy if row["strict_candidate_eligible"]]
    strict_positive = [
        row for row in strict if row["strict_binary_label"] == "positive"
    ]
    strict_negative = [
        row for row in strict if row["strict_binary_label"] == "negative"
    ]
    auc_estimable = bool(strict_positive and strict_negative)
    category_counts = Counter(evidence_category(row) for row in legacy)
    review_complete = (
        second_review is not None
        and not second_review.get("excluded_for_access")
        and all(
            row["independent_second_review"] != "pending"
            for row in evidence_rows
        )
    )
    strict_review = (
        second_review.get("all_records", {}).get(
            "strict_candidate_eligible", {}
        )
        if second_review
        else {}
    )
    binary_review = (
        second_review.get("all_records", {}).get("strict_binary_label", {})
        if second_review
        else {}
    )
    return {
        "figure_version": FIGURE_VERSION,
        "legacy_n": len(legacy),
        "legacy_positive_n": sum(
            row["legacy_direction"] == "positive" for row in legacy
        ),
        "legacy_negative_n": sum(
            row["legacy_direction"] == "negative" for row in legacy
        ),
        "category_counts": {
            category: category_counts.get(category, 0)
            for category in CATEGORY_ORDER
        },
        "strict_candidate_drugs": sorted(str(row["drug"]) for row in strict),
        "strict_positive_n": len(strict_positive),
        "strict_negative_n": len(strict_negative),
        "auc_estimable": auc_estimable,
        "auc": None,
        "independence_status": (
            "not_independent: the same literature informed S_external; "
            "the audit is descriptive and not a performance validation"
        ),
        "second_reviewer_status": (
            "completed_and_adjudicated" if review_complete else "pending"
        ),
        "strict_eligibility_agreement_n": strict_review.get(
            "agreement_count"
        ),
        "strict_eligibility_agreement_total": strict_review.get("n"),
        "strict_eligibility_cohen_kappa": strict_review.get("cohen_kappa"),
        "strict_binary_agreement_n": binary_review.get("agreement_count"),
        "strict_binary_agreement_total": binary_review.get("n"),
        "strict_binary_agreement_exact_ci_95": binary_review.get(
            "agreement_exact_ci_95"
        ),
        "strict_binary_cohen_kappa": binary_review.get("cohen_kappa"),
        "strict_binary_kappa_bootstrap_defined_n": binary_review.get(
            "kappa_bootstrap_defined_n"
        ),
        "strict_binary_kappa_bootstrap_undefined_n": binary_review.get(
            "kappa_bootstrap_undefined_n"
        ),
        "decision": (
            "Retire the legacy 14-drug ROC/PR analysis. The strict drug-specific "
            "clinical subset has two positives and no unambiguous negatives, so "
            "ROC-AUC and PR-AUC are not estimable."
        ),
    }


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _set_publication_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8.5,
            "axes.titlesize": 9.5,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "legend.fontsize": 8.0,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def _minimum_visible_font(fig: plt.Figure) -> float:
    sizes = [
        float(text.get_fontsize())
        for text in fig.findobj(match=matplotlib.text.Text)
        if text.get_visible() and text.get_text().strip()
    ]
    return min(sizes) if sizes else float("nan")


def _export_figure(fig: plt.Figure, basename: Path) -> dict[str, Path]:
    basename.parent.mkdir(parents=True, exist_ok=True)
    paths = {
        "pdf": basename.with_suffix(".pdf"),
        "svg": basename.with_suffix(".svg"),
        "png": basename.with_suffix(".png"),
    }
    fig.savefig(paths["pdf"], format="pdf")
    fig.savefig(paths["svg"], format="svg")
    fig.savefig(paths["png"], format="png", dpi=PNG_DPI)
    return paths


def _inspect_png(path: Path) -> dict[str, Any]:
    from PIL import Image

    with Image.open(path) as image:
        dpi = image.info.get("dpi", (None, None))
        dpi_x = float(dpi[0]) if dpi and dpi[0] else None
        width_mm = image.size[0] / dpi_x * 25.4 if dpi_x else None
        return {
            "pixels": list(image.size),
            "dpi": [float(value) if value else None for value in dpi],
            "measured_width_mm": round(width_mm, 3) if width_mm else None,
            "dpi_check": dpi_x is not None and dpi_x >= PNG_DPI - 1,
            "width_check": width_mm is not None
            and abs(width_mm - TARGET_WIDTH_MM) <= 0.5,
            "file_size_bytes": path.stat().st_size,
        }


def plot_figure4(
    audit_rows: Sequence[Mapping[str, Any]],
    metrics: Mapping[str, Any],
    output_base: Path,
) -> tuple[dict[str, Path], dict[str, Any]]:
    _set_publication_style()
    fig, (left, right) = plt.subplots(
        1,
        2,
        figsize=(TARGET_WIDTH_MM / 25.4, 112.0 / 25.4),
        gridspec_kw={"width_ratios": [1.0, 1.28]},
        constrained_layout=True,
    )

    wrapped_categories = {
        "Drug-specific clinical monotherapy": "Drug-specific clinical\nmonotherapy",
        "Clinical regimen/case (confounded)": "Clinical regimen/case\n(confounded)",
        "Direct ACC preclinical": "Direct ACC preclinical",
        "Class extrapolation": "Class extrapolation",
    }
    counts = [int(metrics["category_counts"][category]) for category in CATEGORY_ORDER]
    row_positions = np.arange(len(CATEGORY_ORDER))[::-1]
    block_height = 0.62
    for position, count, category in zip(
        row_positions, counts, CATEGORY_ORDER
    ):
        block = Rectangle(
            (0, position - block_height / 2),
            count,
            block_height,
            facecolor=CATEGORY_COLORS[category],
            edgecolor="white",
            linewidth=0.8,
        )
        left.add_patch(block)
        left.text(
            count + 0.12,
            position,
            f"n={count}",
            ha="left",
            va="center",
            fontsize=8.3,
            fontweight="bold",
        )
    left.set_yticks(
        row_positions,
        labels=[wrapped_categories[category] for category in CATEGORY_ORDER],
    )
    left.set_xlim(0, 7.1)
    left.set_ylim(-1.55, 3.65)
    left.set_xlabel("Number of legacy labels (exact count)")
    left.set_title("a  Legacy labels are heterogeneous", loc="left", fontweight="bold")
    left.spines[["top", "right", "left"]].set_visible(False)
    left.tick_params(axis="y", length=0)
    left.grid(axis="x", color="#ECECEC", linewidth=0.6)
    left.set_axisbelow(True)
    left.text(
        0.0,
        -0.86,
        (
            "Strict binary clinical subset\n"
            f"{metrics['strict_positive_n']} positive  |  "
            f"{metrics['strict_negative_n']} negative\n"
            "ROC-AUC and PR-AUC not estimable\n"
            "Claude–human traceability: "
            f"{metrics['strict_binary_agreement_n']}/"
            f"{metrics['strict_binary_agreement_total']}\n"
            f"Exact 95% CI "
            f"{metrics['strict_binary_agreement_exact_ci_95'][0]:.3f}"
            f"–{metrics['strict_binary_agreement_exact_ci_95'][1]:.3f}\n"
            "Model-based, not human inter-rater reliability"
        ),
        ha="left",
        va="center",
        fontsize=8.4,
        fontweight="bold",
        color=DARK,
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "#F7F7F7",
            "edgecolor": "#B8B8B8",
            "linewidth": 0.8,
        },
    )

    ordered = sorted(audit_rows, key=lambda row: int(row["rank_comp"]))
    y_positions = np.arange(len(ordered))[::-1]
    for y_value, row in zip(y_positions, ordered):
        rank = int(row["rank_comp"])
        category = str(row["evidence_category"])
        is_legacy_negative = row["legacy_direction"] == "negative"
        if is_legacy_negative:
            right.scatter(
                rank,
                y_value,
                s=46,
                marker="X",
                color=CATEGORY_COLORS[category],
                linewidth=0.8,
                zorder=3,
            )
        else:
            right.scatter(
                rank,
                y_value,
                s=42,
                marker="o",
                facecolor=CATEGORY_COLORS[category],
                edgecolor="white",
                linewidth=0.7,
                zorder=3,
            )
        right.hlines(
            y_value,
            1,
            rank,
            color=LIGHT_GREY,
            linewidth=0.8,
            zorder=1,
        )
    right.axvspan(1, 20, color="#F2F7FA", zorder=0)
    right.axvline(20, color="#7A7A7A", linestyle="--", linewidth=0.8)
    right.set_yticks(
        y_positions,
        labels=[str(row["drug"]) for row in ordered],
    )
    right.set_xlim(0, 109)
    right.set_xlabel("Frozen ADRS_comp rank (1 = highest priority)")
    right.set_title("b  Evidence type on the frozen ranking", loc="left", fontweight="bold")
    right.grid(axis="x", color="#ECECEC", linewidth=0.6)
    right.spines[["top", "right", "left"]].set_visible(False)
    right.tick_params(axis="y", length=0)
    right.text(
        0.985,
        0.985,
        "Circles = old positive label\nX = old negative label",
        transform=right.transAxes,
        ha="right",
        va="top",
        fontsize=8.0,
        color="#555555",
    )
    min_font = _minimum_visible_font(fig)
    paths = _export_figure(fig, output_base)
    qa = {
        "target_width_mm": TARGET_WIDTH_MM,
        "actual_width_mm": round(fig.get_size_inches()[0] * 25.4, 3),
        "actual_height_mm": round(fig.get_size_inches()[1] * 25.4, 3),
        "minimum_font_pt": min_font,
        "font_check": min_font >= MIN_FONT_PT,
        "annotation_source": (
            "data/evidence/evidence_labels_v3_adjudicated.csv, final blinded "
            "second-review agreement, and frozen C1 ranks"
        ),
    }
    plt.close(fig)
    return paths, qa


def _write_locator_audit(
    path: Path, evidence_rows: Sequence[Mapping[str, Any]]
) -> None:
    lines = [
        "# C4 source-locator audit",
        "",
        "Judgment states: `supports`, `partial`, or `unsupported`. The judgment "
        "refers to the legacy drug-specific claim, not merely to whether the cited "
        "paper exists.",
        "",
        "| Drug | Legacy claim | Source support | Locator judgment | Benchmark consequence |",
        "|---|---|---|---|---|",
    ]
    for row in evidence_rows:
        if row["strict_candidate_eligible"]:
            judgment = "supports"
            consequence = f"Retain as strict {row['strict_binary_label']} candidate"
        elif row["evidence_domain"] == "preclinical":
            judgment = "partial"
            consequence = "Retain as preclinical context; exclude from clinical benchmark"
        elif row["evidence_domain"] == "class_extrapolation":
            judgment = "unsupported"
            consequence = "Remove drug-specific clinical label"
        else:
            judgment = "partial"
            consequence = "Retain regimen/context evidence; exclude from strict binary benchmark"
        legacy_claim = (
            f"{row['legacy_direction']} ({float(row['legacy_score']):.2f})"
            if row["legacy_benchmark_member"]
            else f"external score only ({float(row['legacy_score']):.2f})"
        )
        support = str(row["locator_or_support"]).replace("|", "\\|")
        lines.append(
            f"| {row['drug']} | {legacy_claim} | {support} | "
            f"{judgment} | {consequence} |"
        )
    lines.extend(
        [
            "",
            "All 19 records were checked against a primary paper or an ACC clinical "
            "guideline. Anthropic Claude performed a protocol-locked, blinded "
            "model-based classification of all 19 records; this was not a human "
            "rereview. After source-access completion for B02, strict "
            "eligibility and strict binary labels agreed in 19/19 records "
            "in the Claude–human traceability comparison.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def _write_data_card(path: Path, metrics: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "# C4 evidence-label data card",
                "",
                "## Intended use",
                "",
                "The table supports evidence-aware interpretation and reprioritization. "
                "It is not an independent clinical validation dataset.",
                "",
                "## Unit and scope",
                "",
                "- Unit: one named drug.",
                "- Scope: 19 drugs with curated ACC evidence; 14 appeared in the legacy binary benchmark.",
                "- Primary score join: frozen 108-drug complete-case ADRS_comp ranking.",
                "",
                "## Evidence dimensions",
                "",
                "- Domain: clinical, preclinical, or class extrapolation.",
                "- Design: guideline/standard, prospective phase II, randomized regimen, "
                "retrospective series, case report, in-vitro experiment, or no direct study.",
                "- Specificity: direct drug evidence versus class-only extrapolation.",
                "- Exposure: monotherapy/standard versus combination-regimen context.",
                "- Direction: positive, positive preclinical, mixed/limited, regimen effect, "
                "contextual case, class extrapolation, or no direct evidence.",
                "",
                "## Benchmark gate",
                "",
                f"- Legacy labels: {metrics['legacy_positive_n']} positive and "
                f"{metrics['legacy_negative_n']} negative (n={metrics['legacy_n']}).",
                f"- Strict drug-specific clinical candidates: "
                f"{metrics['strict_positive_n']} positive and "
                f"{metrics['strict_negative_n']} negative.",
                "- ROC-AUC requires both classes. Because the strict negative class is empty, "
                "ROC-AUC and PR-AUC are not estimable.",
                "",
                "## Independence and review status",
                "",
                "- The same literature informed S_external; the table cannot be presented as "
                "an independent validation set.",
                "- Source/locator review was completed by the primary reviewer workflow.",
                "- Anthropic Claude performed a protocol-locked, blinded model-based "
                "classification of all 19 records; this was not a human rereview. "
                "Strict eligibility and strict binary labels agreed in "
                f"{metrics['strict_binary_agreement_n']}/"
                f"{metrics['strict_binary_agreement_total']} records.",
                "- Six residual field disagreements across five records were "
                "adjudicated as taxonomy-only; none changed strict inclusion.",
                "- Continuous S_external scores were not independently rescored; "
                "Monte-Carlo perturbation remains a sensitivity analysis, not an "
                "independent rescoring of the continuous rubric.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def _upsert_section(path: Path, heading: str, body: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.is_file() else ""
    marker = f"## {heading}"
    if marker in existing:
        before, remainder = existing.split(marker, 1)
        next_heading = remainder.find("\n## ")
        after = remainder[next_heading:] if next_heading >= 0 else ""
        updated = before.rstrip() + "\n\n" + marker + "\n" + body.rstrip() + "\n" + after
    else:
        updated = existing.rstrip() + "\n\n" + marker + "\n" + body.rstrip() + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(updated.lstrip(), encoding="utf-8")


def _write_figure_plan_card(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "# Figure 4 plan card",
                "",
                "- figure_id: F4",
                "- purpose: expose the evidence heterogeneity that invalidates the legacy binary clinical AUC",
                "- claim: the strict drug-specific clinical subset has two positives and no unambiguous negatives; AUC is not estimable",
                "- data_required: evidence_labels_v3_adjudicated.csv, final blinded second-review agreement, plus frozen C1 ADRS_comp ranks",
                "- layout: (a) evidence-category audit; (b) evidence-type overlay on the 108-drug rank axis",
                "- source_card: new_canonical_candidate -> databases/db07-figures/resources_real.md",
                "- target_journal: mdpi",
                "- column: full",
                "- target_width_mm: 170",
                "- output_formats: PDF, SVG, PNG",
                "- minimum_font_pt: 8",
                "- caption_draft: Figure 4. Adjudicated audit of the heterogeneous evidence labels previously used as a binary clinical benchmark.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def _write_manifest(
    path: Path,
    project_root: Path,
    figure_paths: Mapping[str, Path],
    qa: Mapping[str, Any],
) -> None:
    caption = (
        "Figure 4. Adjudicated audit of the heterogeneous evidence labels previously used "
        "as a binary clinical benchmark. (a) The 14 legacy labels comprise "
        "drug-specific clinical monotherapy evidence, combination-regimen or "
        "case evidence, direct ACC preclinical evidence, and class extrapolation. "
        "Only mitotane and cabozantinib remain drug-specific clinical positive "
        "candidates, with no unambiguous drug-specific clinical negative comparator; "
        "therefore ROC-AUC and PR-AUC are not estimable. A protocol-locked, "
        "blinded Claude classification matched the primary human curator in "
        "19/19 strict eligibility and strict binary labels (exact 95% CI, "
        "0.824–1.000); this is model–human traceability, not human inter-rater "
        "reliability. "
        "residual disagreements were "
        "taxonomy-only. (b) The same evidence "
        "categories are overlaid on the frozen 108-drug ADRS_comp ranking. This "
        "figure is an evidence audit, not an independent performance validation."
    )
    body = "\n".join(
        [
            "",
            "- figure_id: `F4`",
            "- source_card: `new_canonical_candidate -> databases/db07-figures/resources_real.md`",
            f"- vector PDF: `{figure_paths['pdf'].relative_to(project_root).as_posix()}`",
            f"- vector SVG: `{figure_paths['svg'].relative_to(project_root).as_posix()}`",
            f"- bitmap PNG: `{figure_paths['png'].relative_to(project_root).as_posix()}`",
            "- section: Results evidence-label audit",
            "- target: `mdpi`, `full`, 170 mm",
            f"- checks: width={'pass' if qa['png']['width_check'] else 'fail'}; "
            f"dpi={'pass' if qa['png']['dpi_check'] else 'fail'}; "
            f"font={'pass' if qa['figure']['font_check'] else 'fail'}",
            f"- caption: {caption}",
            "",
        ]
    )
    _upsert_section(path, "F4", body)


def _write_version_history(path: Path) -> None:
    body = "\n".join(
        [
            "",
            "- Replaced the legacy 14-drug ROC/PR benchmark with an evidence audit.",
            "- Separated drug-specific clinical monotherapy evidence, regimen/case "
            "evidence, direct ACC preclinical evidence, and class extrapolation.",
            "- Strict binary clinical gate retains two positive candidates "
            "(mitotane and cabozantinib) and no unambiguous negative comparator; "
            "AUC is therefore not estimable.",
            "- The evidence table is explicitly non-independent because the same "
            "literature informed S_external.",
            "- Anthropic Claude completed a protocol-locked, blinded model-based "
            "classification of all 19 records; after B02 source-access completion, "
            "strict eligibility and strict binary labels matched the primary human "
            "curator in 19/19 records. This was not a human rereview.",
            "- Six residual field disagreements across five records were adjudicated "
            "as taxonomy-only; four direction labels were updated in v3 and no "
            "strict inclusion decision changed.",
            "",
        ]
    )
    _upsert_section(
        path,
        "V3 evidence audit and Claude model traceability review — 2026-07-28",
        body,
    )


def _write_change_report(
    path: Path,
    project_root: Path,
    metrics: Mapping[str, Any],
    outputs: Sequence[Path],
) -> None:
    category_lines = [
        f"- {category}: {metrics['category_counts'][category]}"
        for category in CATEGORY_ORDER
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "# C4 change report: evidence-label audit",
                "",
                "## Decision",
                "",
                "The legacy 14-drug ROC-AUC/PR-AUC analysis is retired. Its labels "
                "are not a coherent drug-specific clinical endpoint.",
                "",
                "## Audit result",
                "",
                *category_lines,
                f"- Strict drug-specific clinical subset: "
                f"{metrics['strict_positive_n']} positive, "
                f"{metrics['strict_negative_n']} negative.",
                "- Consequence: ROC-AUC and PR-AUC are not estimable; no replacement "
                "numeric discrimination claim is reported.",
                "",
                "## Corrected interpretations",
                "",
                "- Etoposide, doxorubicin and cisplatin inherit only regimen-level "
                "EDP-M evidence; individual effects are not identified by FIRM-ACT.",
                "- Gemcitabine and erlotinib were evaluated largely as combinations.",
                "- Palbociclib and ribociclib provide direct preclinical, not clinical, evidence.",
                "- Afatinib and osimertinib were mislabeled from gefitinib class extrapolation.",
                "- Sunitinib is reclassified as limited/mixed rather than an unambiguous negative.",
                "- Carboplatin cannot inherit cisplatin evidence.",
                "",
                "## Manuscript synchronization map",
                "",
                "- Remove the AUC 0.40, bootstrap interval, permutation P value, "
                "PR-AUC, subset AUC and leave-one-label-out claims from the abstract, "
                "Methods 2.11, Results, Discussion, Conclusion and Limitations.",
                "- Replace the old Figure 4 file and legend with the evidence-audit "
                "figure and caption in the figure manifest.",
                "- Distinguish the blinded Claude–human categorical traceability "
                "comparison from the continuous S_external scores, which were "
                "not independently rescored.",
                "- Report clinical regimen evidence, preclinical evidence and class "
                "extrapolation separately wherever individual drugs are discussed.",
                "",
                "## Integrity status",
                "",
                "- The same literature informed S_external, so this evidence set is "
                "not described as independent validation.",
                "- All 19 source records have a primary-source/guideline locator review.",
                "- The 12 unique DOIs resolve to matching title/year records in "
                "Crossref and OpenAlex; no scripted retraction signal was found. "
                "This signal check does not guarantee that every historical "
                "retraction database is covered.",
                "- The blinded Claude model classification is complete: strict "
                "eligibility and strict binary labels matched the primary human "
                "curator in 19/19 records. No human rereview was performed.",
                "- Six residual field disagreements across five records were "
                "taxonomy-only and did not change strict inclusion.",
                "",
                "## Files",
                "",
                *[
                    f"- `{output.relative_to(project_root).as_posix()}`"
                    for output in outputs
                ],
                "",
            ]
        ),
        encoding="utf-8",
    )


def generate_c4_outputs(project_root: Path) -> dict[str, Any]:
    project_root = project_root.resolve()
    evidence_path = (
        project_root
        / "data"
        / "evidence"
        / "evidence_labels_v3_adjudicated.csv"
    )
    second_review_path = (
        project_root
        / "results"
        / "evidence_audit"
        / "second_reviewer_agreement.json"
    )
    primary_path = (
        project_root
        / "results"
        / "primary_analysis"
        / "ADRS_comp_primary_108.csv"
    )
    evidence_rows = load_evidence_labels(evidence_path)
    second_review = load_second_review_summary(second_review_path)
    primary_scores = load_primary_scores(primary_path)
    audit_rows = build_audit_rows(evidence_rows, primary_scores, legacy_only=True)
    metrics = evaluate_benchmark(evidence_rows, second_review)

    results_dir = project_root / "results" / "evidence_audit"
    figure_data_path = (
        project_root
        / "figure_data"
        / "revision"
        / "Fig4_evidence_audit_primary108.csv"
    )
    metrics_path = results_dir / "evidence_audit_metrics.json"
    locator_path = results_dir / "source_locator_audit.md"
    data_card_path = results_dir / "data_card.md"
    plan_card_path = (
        project_root
        / "projects"
        / "ACC-PHARMA-NET"
        / "figures"
        / "F4_plan_card.md"
    )
    _write_csv(figure_data_path, audit_rows)
    results_dir.mkdir(parents=True, exist_ok=True)
    metrics_path.write_text(
        json.dumps(metrics, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    _write_locator_audit(locator_path, evidence_rows)
    _write_data_card(data_card_path, metrics)
    _write_figure_plan_card(plan_card_path)

    figure_paths, figure_qa = plot_figure4(
        audit_rows,
        metrics,
        project_root
        / "figures"
        / "revision"
        / "Fig4_evidence_audit_primary108",
    )
    png_qa = _inspect_png(figure_paths["png"])
    qa = {
        "figure_version": FIGURE_VERSION,
        "python": platform.python_version(),
        "matplotlib": matplotlib.__version__,
        "figure": figure_qa,
        "png": png_qa,
    }
    if not figure_qa["font_check"]:
        raise ValueError(f"Figure 4 contains text below {MIN_FONT_PT} pt")
    if not png_qa["dpi_check"]:
        raise ValueError(f"Figure 4 PNG does not meet {PNG_DPI} dpi")
    if not png_qa["width_check"]:
        raise ValueError(f"Figure 4 PNG is not {TARGET_WIDTH_MM} mm wide")
    qa_path = results_dir / "C4_figure_QA.json"
    qa_path.write_text(
        json.dumps(qa, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    manifest_path = (
        project_root
        / "projects"
        / "ACC-PHARMA-NET"
        / "figures"
        / "manifest.md"
    )
    _write_manifest(manifest_path, project_root, figure_paths, qa)
    version_history_path = (
        project_root
        / "databases"
        / "db09-projects"
        / "projects"
        / "ACC-PHARMA-NET"
        / "version_history.md"
    )
    _write_version_history(version_history_path)

    report_path = results_dir / "C4_change_report.md"
    outputs = [
        evidence_path,
        second_review_path,
        figure_data_path,
        metrics_path,
        locator_path,
        data_card_path,
        plan_card_path,
        *figure_paths.values(),
        qa_path,
        manifest_path,
        version_history_path,
    ]
    citation_verification_path = results_dir / "citation_verification.json"
    if citation_verification_path.is_file():
        outputs.append(citation_verification_path)
    _write_change_report(report_path, project_root, metrics, outputs)
    return {
        "metrics": metrics,
        "figure_paths": figure_paths,
        "figure_data": figure_data_path,
        "qa": qa,
        "qa_path": qa_path,
        "locator_audit": locator_path,
        "data_card": data_card_path,
        "plan_card": plan_card_path,
        "manifest": manifest_path,
        "version_history": version_history_path,
        "report": report_path,
    }


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    outputs = generate_c4_outputs(args.project_root)
    metrics = outputs["metrics"]
    print(
        "C4 evidence audit complete: "
        f"legacy n={metrics['legacy_n']}, "
        f"strict positive={metrics['strict_positive_n']}, "
        f"strict negative={metrics['strict_negative_n']}, "
        f"AUC estimable={metrics['auc_estimable']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
