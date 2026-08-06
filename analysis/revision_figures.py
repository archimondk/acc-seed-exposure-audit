"""Generate model-matched revision figures from the frozen C1–C2 outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import blended_transform_factory
from scipy.stats import spearmanr


FIGURE_VERSION = "c3-primary108-v1"
TARGET_JOURNAL = "mdpi"
TARGET_COLUMN = "full"
TARGET_WIDTH_MM = 170.0
MIN_FONT_PT = 8.0
PNG_DPI = 1000
WEIGHTS = tuple(index / 20 for index in range(21))
PRIMARY_WEIGHT = 0.5
SELECTED_DRUGS = (
    "Irinotecan",
    "Ixazomib",
    "Abemaciclib",
    "Olaparib",
    "Mitotane",
    "Doxorubicin",
    "Palbociclib",
    "Afatinib",
    "Ribociclib",
)
CORRELATION_VARIABLES = (
    ("C_ACC_pct", "C_ACC\npercentile"),
    ("residual_pct", "ACC-relative\nresidual percentile"),
    ("MIPE_potency_pct", "MIPE potency\npercentile"),
    ("NCI60_potency_pct", "NCI-60 potency\npercentile"),
    ("ADRS_comp", "ADRS_comp"),
)

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
GREY = "#A7A9AC"
LIGHT_BLUE = "#A6CEE3"
DARK = "#202124"


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


def load_c3_inputs(
    project_root: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    project_root = project_root.resolve()
    primary_path = (
        project_root / "results" / "primary_analysis" / "ADRS_comp_primary_108.csv"
    )
    class_path = (
        project_root
        / "results"
        / "mechanism_enrichment"
        / "mechanism_enrichment_primary108.csv"
    )
    primary_raw = _read_csv(
        primary_path,
        (
            "rank_comp",
            "drug",
            "C_ACC_pct",
            "residual_pct",
            "MIPE_potency_pct",
            "NCI60_potency_pct",
            "ADRS_comp",
        ),
    )
    class_raw = _read_csv(
        class_path,
        (
            "mechanism_class",
            "k",
            "n_universe",
            "rank_sum",
            "mean_rank",
            "p_exact",
            "q_bh",
            "member_ranks",
        ),
    )

    numeric_primary = {
        "C_ACC_pct",
        "residual_pct",
        "MIPE_potency_pct",
        "NCI60_potency_pct",
        "ADRS_comp",
    }
    primary_rows: list[dict[str, Any]] = []
    for row in primary_raw:
        parsed: dict[str, Any] = dict(row)
        parsed["rank_comp"] = int(row["rank_comp"])
        for field in numeric_primary:
            parsed[field] = float(row[field])
        primary_rows.append(parsed)
    primary_rows.sort(key=lambda row: row["rank_comp"])
    if len(primary_rows) != 108:
        raise ValueError(f"C3 requires the locked 108-drug universe, got {len(primary_rows)}")
    if [row["rank_comp"] for row in primary_rows] != list(range(1, 109)):
        raise ValueError("C1 ranks must be consecutive integers 1..108")

    class_rows: list[dict[str, Any]] = []
    for row in class_raw:
        parsed = dict(row)
        for field in ("k", "n_universe", "rank_sum"):
            parsed[field] = int(row[field])
        for field in ("mean_rank", "p_exact", "q_bh"):
            parsed[field] = float(row[field])
        parsed["member_rank_values"] = [
            int(value.strip()) for value in row["member_ranks"].split(";")
        ]
        class_rows.append(parsed)
    return primary_rows, class_rows


def _ordinal_ranks(scores: Mapping[str, float]) -> dict[str, int]:
    ordered = sorted(scores, key=lambda drug: (-scores[drug], drug))
    return {drug: index for index, drug in enumerate(ordered, start=1)}


def compute_weight_scan(
    primary_rows: Sequence[Mapping[str, Any]],
    weights: Sequence[float] = WEIGHTS,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    by_drug = {str(row["drug"]): row for row in primary_rows}
    for weight_c in weights:
        weight_residual = 1.0 - float(weight_c)
        scores = {
            drug: (
                float(weight_c) * float(row["C_ACC_pct"])
                + weight_residual * float(row["residual_pct"])
            )
            for drug, row in by_drug.items()
        }
        is_primary = math.isclose(float(weight_c), PRIMARY_WEIGHT, abs_tol=1e-12)
        if is_primary:
            # The C1 table is the ranking authority.  Reusing its locked ranks
            # prevents decimal serialization from changing the order of drugs
            # whose full-precision component sums are extremely close.
            scores = {
                drug: float(row["ADRS_comp"]) for drug, row in by_drug.items()
            }
            ranks = {
                drug: int(row["rank_comp"]) for drug, row in by_drug.items()
            }
        else:
            ranks = _ordinal_ranks(scores)
        for drug in sorted(by_drug, key=lambda item: ranks[item]):
            rows.append(
                {
                    "weight_C_ACC": float(weight_c),
                    "weight_residual": weight_residual,
                    "drug": drug,
                    "score": scores[drug],
                    "rank": ranks[drug],
                    "is_primary_weight": is_primary,
                }
            )
    return rows


def compute_primary_correlation(
    primary_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    matrix_rows: list[dict[str, Any]] = []
    for field_x, label_x in CORRELATION_VARIABLES:
        for field_y, label_y in CORRELATION_VARIABLES:
            x_values: list[float] = []
            y_values: list[float] = []
            for row in primary_rows:
                x = row.get(field_x)
                y = row.get(field_y)
                if x is not None and y is not None:
                    x_values.append(float(x))
                    y_values.append(float(y))
            rho = float(spearmanr(x_values, y_values).statistic)
            matrix_rows.append(
                {
                    "variable_x": field_x,
                    "label_x": label_x.replace("\n", " "),
                    "variable_y": field_y,
                    "label_y": label_y.replace("\n", " "),
                    "n_pairwise": len(x_values),
                    "spearman_rho": rho,
                }
            )
    return matrix_rows


def compute_exact_mean_rank_null(
    n_universe: int, class_size: int
) -> list[dict[str, Any]]:
    """Count the exact distribution of sums of k distinct ranks from 1..N."""

    maximum_sum = class_size * (2 * n_universe - class_size + 1) // 2
    counts: list[list[int]] = [
        [0] * (maximum_sum + 1) for _ in range(class_size + 1)
    ]
    counts[0][0] = 1
    for rank in range(1, n_universe + 1):
        for selected in range(min(class_size, rank), 0, -1):
            for current_sum in range(maximum_sum, rank - 1, -1):
                counts[selected][current_sum] += counts[selected - 1][
                    current_sum - rank
                ]
    total = math.comb(n_universe, class_size)
    return [
        {
            "rank_sum": rank_sum,
            "mean_rank": rank_sum / class_size,
            "combination_count": count,
            "probability": count / total,
            "n_universe": n_universe,
            "class_size": class_size,
        }
        for rank_sum, count in enumerate(counts[class_size])
        if count
    ]


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty data: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def prepare_c3_data(
    project_root: Path, output_dir: Path
) -> tuple[dict[str, Path], dict[str, Any]]:
    primary_rows, class_rows = load_c3_inputs(project_root)
    correlation_rows = compute_primary_correlation(primary_rows)
    weight_rows = compute_weight_scan(primary_rows)
    cdk = next(row for row in class_rows if row["mechanism_class"] == "CDK4/6")
    null_rows = compute_exact_mean_rank_null(cdk["n_universe"], cdk["k"])

    favorable = sum(
        row["combination_count"]
        for row in null_rows
        if row["rank_sum"] <= cdk["rank_sum"]
    )
    total = sum(row["combination_count"] for row in null_rows)
    p_from_null = favorable / total
    if not math.isclose(p_from_null, cdk["p_exact"], abs_tol=1e-15):
        raise ValueError("C3 exact null does not reproduce the locked C2 P value")

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "figure3_data": output_dir / "Fig3_component_correlation_primary108.csv",
        "figure5a_data": output_dir / "Fig5a_weight_scan_primary108.csv",
        "figure5b_data": output_dir / "Fig5b_CDK46_exact_null_primary108.csv",
        "stats": output_dir / "C3_figure_stats.json",
    }
    _write_csv(paths["figure3_data"], correlation_rows)
    _write_csv(paths["figure5a_data"], weight_rows)
    _write_csv(paths["figure5b_data"], null_rows)

    selected_scan = [
        row for row in weight_rows if row["drug"] in set(SELECTED_DRUGS)
    ]
    primary_ranks = {
        str(row["drug"]): int(row["rank_comp"]) for row in primary_rows
    }
    stats = {
        "figure_version": FIGURE_VERSION,
        "n_primary": len(primary_rows),
        "n_weight_settings": len(WEIGHTS),
        "weight_C_ACC_values": list(WEIGHTS),
        "primary_weight_C_ACC": PRIMARY_WEIGHT,
        "selected_drugs": list(SELECTED_DRUGS),
        "selected_scan_rows": len(selected_scan),
        "cdk46_members": cdk["members"].split("; "),
        "cdk46_ranks": cdk["member_rank_values"],
        "cdk46_mean_rank": cdk["mean_rank"],
        "cdk46_p_exact": cdk["p_exact"],
        "cdk46_q_bh": cdk["q_bh"],
        "cdk46_rank_sum": cdk["rank_sum"],
        "n_eligible_classes": len(class_rows),
        "exact_null_combinations": total,
        "primary_selected_ranks": {
            drug: primary_ranks[drug] for drug in SELECTED_DRUGS
        },
    }
    paths["stats"].write_text(
        json.dumps(stats, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return paths, stats


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
            "lines.linewidth": 1.2,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "savefig.bbox": None,
        }
    )


def _export_figure(fig: plt.Figure, basename: Path) -> dict[str, Path]:
    basename.parent.mkdir(parents=True, exist_ok=True)
    paths = {
        "pdf": basename.with_suffix(".pdf"),
        "svg": basename.with_suffix(".svg"),
        "png": basename.with_suffix(".png"),
    }
    fig.savefig(paths["pdf"], format="pdf", bbox_inches=None)
    fig.savefig(paths["svg"], format="svg", bbox_inches=None)
    fig.savefig(paths["png"], format="png", dpi=PNG_DPI, bbox_inches=None)
    return paths


def _minimum_visible_font(fig: plt.Figure) -> float:
    text_objects = fig.findobj(match=matplotlib.text.Text)
    sizes = [
        float(text.get_fontsize())
        for text in text_objects
        if text.get_visible() and text.get_text().strip()
    ]
    return min(sizes) if sizes else float("nan")


def _plot_figure3(
    primary_rows: Sequence[Mapping[str, Any]],
    correlation_rows: Sequence[Mapping[str, Any]],
    output_base: Path,
) -> tuple[dict[str, Path], dict[str, Any]]:
    n_variables = len(CORRELATION_VARIABLES)
    index = {
        (row["variable_y"], row["variable_x"]): float(row["spearman_rho"])
        for row in correlation_rows
    }
    matrix = np.asarray(
        [
            [index[(field_y, field_x)] for field_x, _ in CORRELATION_VARIABLES]
            for field_y, _ in CORRELATION_VARIABLES
        ]
    )
    labels = [label for _, label in CORRELATION_VARIABLES]
    cmap = LinearSegmentedColormap.from_list(
        "blue_grey_orange", [BLUE, "#F5F5F5", ORANGE]
    )

    fig, ax = plt.subplots(
        figsize=(TARGET_WIDTH_MM / 25.4, 124.0 / 25.4),
        constrained_layout=True,
    )
    image = ax.imshow(matrix, cmap=cmap, vmin=-1, vmax=1, aspect="equal")
    ax.set_xticks(range(n_variables), labels=labels, rotation=35, ha="right")
    ax.set_yticks(range(n_variables), labels=labels)
    for row_index in range(n_variables):
        for column_index in range(n_variables):
            value = matrix[row_index, column_index]
            ax.text(
                column_index,
                row_index,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=8.5,
                color="white" if abs(value) >= 0.58 else DARK,
            )
    ax.set_title(
        f"Primary-model correlation structure (Spearman, n={len(primary_rows)})",
        pad=9,
        fontweight="bold",
    )
    colorbar = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.04)
    colorbar.set_label("Spearman correlation")
    colorbar.set_ticks([-1, -0.5, 0, 0.5, 1])
    min_font = _minimum_visible_font(fig)
    paths = _export_figure(fig, output_base)
    qa = {
        "target_width_mm": TARGET_WIDTH_MM,
        "actual_width_mm": round(fig.get_size_inches()[0] * 25.4, 3),
        "actual_height_mm": round(fig.get_size_inches()[1] * 25.4, 3),
        "minimum_font_pt": min_font,
        "font_check": min_font >= MIN_FONT_PT,
        "n_primary": len(primary_rows),
        "color_scale": [-1, 1],
    }
    plt.close(fig)
    return paths, qa


def _plot_figure5(
    primary_rows: Sequence[Mapping[str, Any]],
    weight_rows: Sequence[Mapping[str, Any]],
    null_rows: Sequence[Mapping[str, Any]],
    stats: Mapping[str, Any],
    output_base: Path,
) -> tuple[dict[str, Path], dict[str, Any]]:
    selected_set = set(SELECTED_DRUGS)
    values_by_drug: dict[str, list[int]] = {drug: [] for drug in SELECTED_DRUGS}
    primary_rank_by_drug: dict[str, int] = {}
    for row in weight_rows:
        drug = str(row["drug"])
        if drug not in selected_set:
            continue
        values_by_drug[drug].append(int(row["rank"]))
        if row["is_primary_weight"]:
            primary_rank_by_drug[drug] = int(row["rank"])

    order = sorted(SELECTED_DRUGS, key=lambda drug: primary_rank_by_drug[drug])
    box_data = [values_by_drug[drug] for drug in order]
    means = np.asarray([float(row["mean_rank"]) for row in null_rows])
    probabilities = np.asarray([float(row["probability"]) for row in null_rows])
    observed = float(stats["cdk46_mean_rank"])
    tail_mask = means <= observed + 1e-12

    fig, (left, right) = plt.subplots(
        1,
        2,
        figsize=(TARGET_WIDTH_MM / 25.4, 96.0 / 25.4),
        gridspec_kw={"width_ratios": [1.48, 1.0]},
        constrained_layout=True,
    )
    box = left.boxplot(
        box_data,
        vert=False,
        tick_labels=order,
        widths=0.62,
        showfliers=False,
        patch_artist=True,
        medianprops={"color": DARK, "linewidth": 1.2},
        whiskerprops={"color": BLUE, "linewidth": 1.0},
        capprops={"color": BLUE, "linewidth": 1.0},
    )
    for patch in box["boxes"]:
        patch.set(facecolor=LIGHT_BLUE, edgecolor=BLUE, linewidth=1.0)
    positions = np.arange(1, len(order) + 1)
    left.scatter(
        [primary_rank_by_drug[drug] for drug in order],
        positions,
        marker="D",
        s=24,
        facecolor=ORANGE,
        edgecolor="white",
        linewidth=0.5,
        zorder=4,
    )
    left.axvline(20, color="#666666", linestyle="--", linewidth=0.9)
    top20_transform = blended_transform_factory(left.transData, left.transAxes)
    left.text(
        20,
        0.985,
        "Top 20",
        transform=top20_transform,
        color="#555555",
        fontsize=8.0,
        ha="center",
        va="top",
    )
    left.set_xlim(0.5, 108.5)
    left.invert_yaxis()
    left.set_xlabel("Rank across 21 two-component weights (lower = higher priority)")
    left.set_title(
        "a  Primary-score weight sensitivity",
        loc="left",
        fontweight="bold",
    )
    left.grid(axis="x", color="#E5E7EB", linewidth=0.6)
    left.spines[["top", "right"]].set_visible(False)
    left.text(
        0.985,
        0.985,
        r"Orange diamonds = primary model ($w_C=0.5$)",
        transform=left.transAxes,
        ha="right",
        va="top",
        fontsize=8.0,
        color=ORANGE,
        bbox={
            "boxstyle": "square,pad=0.15",
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.85,
        },
    )

    # Exact probability mass function, not an SD/SEM/CI uncertainty band.
    right.fill_between(
        means,
        probabilities,
        0,
        where=~tail_mask,
        step="mid",
        interpolate=True,
        color=GREY,
    )
    right.fill_between(
        means,
        probabilities,
        0,
        where=tail_mask,
        step="mid",
        interpolate=True,
        color=LIGHT_BLUE,
    )
    right.axvline(observed, color=GREEN, linewidth=1.8)
    right.set_xlabel("Mean rank of a three-drug set")
    right.set_ylabel("Exact probability")
    right.set_title(
        "b  CDK4/6 random-set assessment",
        loc="left",
        fontweight="bold",
    )
    right.text(
        0.98,
        0.045,
        (
            f"Ranks {stats['cdk46_ranks'][0]}, {stats['cdk46_ranks'][1]}, "
            f"{stats['cdk46_ranks'][2]}\n"
            f"Mean = {observed:.2f}\n"
            f"Exact P = {stats['cdk46_p_exact']:.4f}\n"
            f"BH q = {stats['cdk46_q_bh']:.4f}\n"
            "Not significant"
        ),
        transform=right.transAxes,
        ha="right",
        va="bottom",
        fontsize=8.2,
        bbox={
            "boxstyle": "round,pad=0.3",
            "facecolor": "white",
            "edgecolor": "#D1D5DB",
            "linewidth": 0.7,
        },
    )
    right.text(
        0.17,
        0.22,
        r"lower-tail null mass",
        transform=right.transAxes,
        fontsize=8.0,
        color=BLUE,
        ha="center",
    )
    right.text(
        0.73,
        0.72,
        "exact null",
        transform=right.transAxes,
        fontsize=8.0,
        color="#666666",
        ha="center",
    )
    right.spines[["top", "right"]].set_visible(False)
    min_font = _minimum_visible_font(fig)
    paths = _export_figure(fig, output_base)
    qa = {
        "target_width_mm": TARGET_WIDTH_MM,
        "actual_width_mm": round(fig.get_size_inches()[0] * 25.4, 3),
        "actual_height_mm": round(fig.get_size_inches()[1] * 25.4, 3),
        "minimum_font_pt": min_font,
        "font_check": min_font >= MIN_FONT_PT,
        "n_weight_settings": len(WEIGHTS),
        "exact_null_combinations": stats["exact_null_combinations"],
        "annotation_source": "computed C1-C2 objects",
    }
    plt.close(fig)
    return paths, qa


def _inspect_png(path: Path) -> dict[str, Any]:
    from PIL import Image

    with Image.open(path) as image:
        dpi = image.info.get("dpi", (None, None))
        dpi_x = float(dpi[0]) if dpi and dpi[0] else None
        width_mm = (
            image.size[0] / dpi_x * 25.4 if dpi_x is not None else None
        )
        return {
            "pixels": list(image.size),
            "dpi": [float(value) if value is not None else None for value in dpi],
            "measured_width_mm": round(width_mm, 3) if width_mm else None,
            "dpi_check": dpi_x is not None and dpi_x >= PNG_DPI - 1,
            "width_check": width_mm is not None
            and abs(width_mm - TARGET_WIDTH_MM) <= 0.5,
            "file_size_bytes": path.stat().st_size,
        }


def _write_manifest(
    project_root: Path,
    figure_paths: Mapping[str, Mapping[str, Path]],
    qa: Mapping[str, Any],
) -> Path:
    manifest_path = (
        project_root
        / "projects"
        / "ACC-PHARMA-NET"
        / "figures"
        / "manifest.md"
    )
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    captions = {
        "F3": (
            "Figure 3. Correlation structure of the primary model in the "
            "108-drug complete-case universe. Cells show Spearman correlations "
            "among the disease-context percentile, ACC-relative activity-residual "
            "percentile, MIPE potency percentile, NCI-60 potency percentile and "
            "ADRS_comp. Correlation describes dependence and does not establish "
            "independent predictive contribution."
        ),
        "F5": (
            "Figure 5. Primary-model weight sensitivity and exact CDK4/6 random-set "
            "assessment. (a) Rank distributions for selected drugs across 21 "
            "two-component settings, with diamonds marking the primary equal-weight "
            "model. (b) Exact null distribution of the mean rank for all 204,156 "
            "three-drug subsets of the same 108-drug universe. Abemaciclib, "
            "palbociclib and ribociclib ranked 8, 26 and 51 (mean 28.33; exact "
            "one-sided P=0.0764; BH q=0.3711 across ten eligible mechanism "
            "families), indicating a non-significant ranking trend."
        ),
    }
    source_cards = {
        "F3": "databases/db07-figures/resources_real.md::相关矩阵（数值注释热图）",
        "F5": "databases/db07-figures/resources_real.md::参数敏感性分布与随机化检验组图",
    }
    sections = {"F3": "Results 3.6", "F5": "Results 3.6 and 3.8"}
    lines = ["# ACC-PHARMA-NET figure manifest", ""]
    for figure_id in ("F3", "F5"):
        paths = figure_paths[figure_id]
        lines.extend(
            [
                f"## {figure_id}",
                "",
                f"- figure_id: `{figure_id}`",
                f"- source_card: `{source_cards[figure_id]}`",
                f"- vector PDF: `{paths['pdf'].relative_to(project_root).as_posix()}`",
                f"- vector SVG: `{paths['svg'].relative_to(project_root).as_posix()}`",
                f"- bitmap PNG: `{paths['png'].relative_to(project_root).as_posix()}`",
                f"- section: {sections[figure_id]}",
                f"- target: `{TARGET_JOURNAL}`, `{TARGET_COLUMN}`, "
                f"{TARGET_WIDTH_MM:.0f} mm",
                f"- checks: width={'pass' if qa[figure_id]['png']['width_check'] else 'fail'}; "
                f"dpi={'pass' if qa[figure_id]['png']['dpi_check'] else 'fail'}; "
                f"font={'pass' if qa[figure_id]['figure']['font_check'] else 'fail'}",
                f"- caption: {captions[figure_id]}",
                "",
            ]
        )
    manifest_path.write_text("\n".join(lines), encoding="utf-8")
    return manifest_path


def _write_version_history(project_root: Path) -> Path:
    path = (
        project_root
        / "databases"
        / "db09-projects"
        / "projects"
        / "ACC-PHARMA-NET"
        / "version_history.md"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    entry = "\n".join(
        [
            "# ACC-PHARMA-NET version history",
            "",
            "## C3 figure lock — 2026-07-27",
            "",
            "- Target style: Pharmaceutics/MDPI full width, 170 mm.",
            "- Font: Arial with Helvetica/DejaVu Sans fallback; minimum 8 pt.",
            "- Palette: Okabe–Ito blue `#0072B2`, orange `#D55E00`, "
            "green `#009E73`; neutral grey `#A7A9AC`.",
            "- Figure 3 now uses only the five variables in or directly underlying "
            "the locked 108-drug primary score.",
            "- Figure 5 now uses 21 two-component weights and the exact C2 "
            "CDK4/6 null distribution; all numerical annotations are computed.",
            "- Figure 4 is generated separately by the C4 evidence-label audit.",
            "",
        ]
    )
    path.write_text(entry, encoding="utf-8")
    return path


def _write_c3_report(
    project_root: Path,
    stats: Mapping[str, Any],
    qa: Mapping[str, Any],
    data_paths: Mapping[str, Path],
    figure_paths: Mapping[str, Mapping[str, Path]],
    manifest_path: Path,
) -> Path:
    report_dir = project_root / "results" / "figures"
    report_dir.mkdir(parents=True, exist_ok=True)
    qa_path = report_dir / "C3_figure_QA.json"
    qa_path.write_text(
        json.dumps(qa, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    report_path = report_dir / "C3_change_report.md"
    report_path.write_text(
        "\n".join(
            [
                "# C3 change report: model-matched figures",
                "",
                f"- Figure version: `{FIGURE_VERSION}`.",
                f"- Target: Pharmaceutics/MDPI full width ({TARGET_WIDTH_MM:.0f} mm).",
                "- Figure 3 was rebuilt from the C1 108-drug table; legacy "
                "S_B-neighbor/external columns and pairwise 124-drug mixtures were removed.",
                f"- Figure 5a uses {stats['n_weight_settings']} two-component weights "
                "from w_C=0 to 1 in 0.05 increments.",
                f"- Figure 5b uses all {stats['exact_null_combinations']:,} possible "
                "three-drug subsets, not a Monte-Carlo sample.",
                f"- CDK4/6 annotation is computed as ranks "
                f"{'/'.join(map(str, stats['cdk46_ranks']))}, mean "
                f"{stats['cdk46_mean_rank']:.2f}, exact P="
                f"{stats['cdk46_p_exact']:.4f}, BH q={stats['cdk46_q_bh']:.4f}.",
                "- The title and caption state “non-significant ranking trend”; "
                "the legacy positive-enrichment wording was removed.",
                "- Figure 4 was not regenerated because its evidence-label semantics "
                "are the subject of C4; retaining a temporary figure gap is preferable "
                "to presenting the mixed-evidence set as clinical validation.",
                "",
                "## Files",
                "",
                *[
                    f"- `{path.relative_to(project_root).as_posix()}`"
                    for path in [
                        *data_paths.values(),
                        *figure_paths["F3"].values(),
                        *figure_paths["F5"].values(),
                        qa_path,
                        manifest_path,
                    ]
                ],
                "",
                "## QA",
                "",
                f"- F3: width {qa['F3']['png']['measured_width_mm']} mm; "
                f"PNG dpi {qa['F3']['png']['dpi'][0]:.1f}; minimum font "
                f"{qa['F3']['figure']['minimum_font_pt']:.1f} pt.",
                f"- F5: width {qa['F5']['png']['measured_width_mm']} mm; "
                f"PNG dpi {qa['F5']['png']['dpi'][0]:.1f}; minimum font "
                f"{qa['F5']['figure']['minimum_font_pt']:.1f} pt.",
                "- Vector PDF and SVG retain editable text.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return report_path


def generate_c3_figures(project_root: Path) -> dict[str, Any]:
    project_root = project_root.resolve()
    _set_publication_style()
    data_dir = project_root / "figure_data" / "revision"
    figure_dir = project_root / "figures" / "revision"
    data_paths, stats = prepare_c3_data(project_root, data_dir)
    primary_rows, _ = load_c3_inputs(project_root)

    correlation_rows = compute_primary_correlation(primary_rows)
    weight_rows = compute_weight_scan(primary_rows)
    null_rows = compute_exact_mean_rank_null(108, 3)

    f3_paths, f3_figure_qa = _plot_figure3(
        primary_rows,
        correlation_rows,
        figure_dir / "FigS5_component_correlation_primary108",
    )
    f5_paths, f5_figure_qa = _plot_figure5(
        primary_rows,
        weight_rows,
        null_rows,
        stats,
        figure_dir / "Fig3_weight_stability_CDK46_primary108",
    )
    qa = {
        "figure_version": FIGURE_VERSION,
        "python": platform.python_version(),
        "matplotlib": matplotlib.__version__,
        "target_journal": TARGET_JOURNAL,
        "target_column": TARGET_COLUMN,
        "F3": {
            "figure": f3_figure_qa,
            "png": _inspect_png(f3_paths["png"]),
        },
        "F5": {
            "figure": f5_figure_qa,
            "png": _inspect_png(f5_paths["png"]),
        },
    }
    for figure_id in ("F3", "F5"):
        if not qa[figure_id]["figure"]["font_check"]:
            raise ValueError(f"{figure_id} contains text below {MIN_FONT_PT} pt")
        if not qa[figure_id]["png"]["dpi_check"]:
            raise ValueError(f"{figure_id} PNG does not meet {PNG_DPI} dpi")
        if not qa[figure_id]["png"]["width_check"]:
            raise ValueError(f"{figure_id} PNG is not {TARGET_WIDTH_MM} mm wide")

    figure_paths = {"F3": f3_paths, "F5": f5_paths}
    manifest_path = _write_manifest(project_root, figure_paths, qa)
    version_history_path = _write_version_history(project_root)
    report_path = _write_c3_report(
        project_root,
        stats,
        qa,
        data_paths,
        figure_paths,
        manifest_path,
    )
    return {
        "data_paths": data_paths,
        "figure_paths": figure_paths,
        "stats": stats,
        "qa": qa,
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
    outputs = generate_c3_figures(args.project_root)
    stats = outputs["stats"]
    print(
        f"C3 figures complete: n={stats['n_primary']}, "
        f"weights={stats['n_weight_settings']}, "
        f"CDK4/6 P={stats['cdk46_p_exact']:.6f}, "
        f"q={stats['cdk46_q_bh']:.6f}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
