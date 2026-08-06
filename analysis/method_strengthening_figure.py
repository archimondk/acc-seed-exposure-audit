"""Generate reviewer-facing Figure 2 from method-strengthening outputs."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analysis.revision_figures import (
    BLUE,
    DARK,
    GREEN,
    GREY,
    LIGHT_BLUE,
    MIN_FONT_PT,
    ORANGE,
    PNG_DPI,
    TARGET_COLUMN,
    TARGET_JOURNAL,
    TARGET_WIDTH_MM,
    _export_figure,
    _inspect_png,
    _minimum_visible_font,
    _set_publication_style,
)


FIGURE_VERSION = "method-strengthening-f2-v1"


def _read_csv(path: Path, required: Iterable[str]) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Required figure source is missing: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"Figure source has no rows: {path}")
    missing = set(required) - set(rows[0])
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    return rows


def load_figure_inputs(
    project_root: Path,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, Any],
]:
    root = project_root.resolve()
    result_dir = root / "results" / "method_strengthening"
    baseline_raw = _read_csv(
        result_dir / "baseline_comparison_primary108.csv",
        (
            "ranking",
            "role",
            "spearman_vs_ADRS_comp",
            "top20_jaccard",
            "cdk46_p_exact",
        ),
    )
    gene_raw = _read_csv(
        result_dir / "centrality_gene399.csv",
        (
            "gene",
            "r_ACC",
            "is_ACC_seed",
            "STRING_degree",
            "STRING_PageRank",
        ),
    )
    null_raw = _read_csv(
        result_dir / "random_seed_null_primary108.csv",
        (
            "drug",
            "observed_C_ACC_pct",
            "null_mean_C_ACC_pct",
            "null_sd_C_ACC_pct",
            "z_degree_matched",
            "empirical_p_upper",
            "q_bh_108",
        ),
    )
    metrics = json.loads(
        (result_dir / "method_strengthening_metrics.json").read_text(
            encoding="utf-8"
        )
    )
    baselines = [
        {
            **row,
            "spearman_vs_ADRS_comp": float(row["spearman_vs_ADRS_comp"]),
            "top20_jaccard": float(row["top20_jaccard"]),
            "cdk46_p_exact": float(row["cdk46_p_exact"]),
        }
        for row in baseline_raw
    ]
    genes = [
        {
            **row,
            "r_ACC": float(row["r_ACC"]),
            "is_ACC_seed": row["is_ACC_seed"].casefold() == "true",
            "STRING_degree": int(row["STRING_degree"]),
            "STRING_PageRank": float(row["STRING_PageRank"]),
        }
        for row in gene_raw
    ]
    null_rows = [
        {
            **row,
            "observed_C_ACC_pct": float(row["observed_C_ACC_pct"]),
            "null_mean_C_ACC_pct": float(row["null_mean_C_ACC_pct"]),
            "null_sd_C_ACC_pct": float(row["null_sd_C_ACC_pct"]),
            "z_degree_matched": float(row["z_degree_matched"]),
            "empirical_p_upper": float(row["empirical_p_upper"]),
            "q_bh_108": float(row["q_bh_108"]),
        }
        for row in null_raw
    ]
    if len(baselines) != 8 or len(genes) != 399 or len(null_rows) != 108:
        raise ValueError("Figure 2 requires the frozen 8/399/108 source rows")
    if metrics.get("null_draws") != 10_000:
        raise ValueError("Figure 2 requires the formal 10,000-draw null")
    return baselines, genes, null_rows, metrics


def _panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.15,
        1.08,
        label,
        transform=axis.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
        ha="left",
        color=DARK,
    )


def _plot_baselines(axis: plt.Axes, rows: Sequence[Mapping[str, Any]]) -> None:
    labels = {
        "raw_MIPE_potency": "Raw MIPE",
        "residual_alone": "Residual",
        "C_ACC_alone": "C_ACC",
        "association_count": "Association count",
        "direct_seed_overlap_fraction": "Direct seed overlap",
        "S_external": "S_external",
        "degree_matched_random_seed_mean": "Random-seed mean",
    }
    label_offsets = {
        "raw_MIPE_potency": (16, -22),
        "residual_alone": (16, 28),
        "C_ACC_alone": (16, -2),
        "association_count": (7, 8),
        "direct_seed_overlap_fraction": (-44, 20),
        "S_external": (7, 18),
        "degree_matched_random_seed_mean": (18, -14),
    }
    colors = {
        "model_component": BLUE,
        "simple_baseline": GREY,
        "exploratory_non_independent": ORANGE,
        "random_network_baseline": GREEN,
    }
    markers = {
        "model_component": "o",
        "simple_baseline": "s",
        "exploratory_non_independent": "^",
        "random_network_baseline": "D",
    }
    role_labels = {
        "simple_baseline": "Simple baseline",
        "model_component": "Model component",
        "exploratory_non_independent": "Literature (non-independent)",
        "random_network_baseline": "Random network",
    }
    plotted_roles: set[str] = set()
    for row in rows:
        if row["ranking"] == "ADRS_comp":
            continue
        role = str(row["role"])
        axis.scatter(
            row["spearman_vs_ADRS_comp"],
            row["top20_jaccard"],
            s=42,
            marker=markers[role],
            facecolor=colors[role],
            edgecolor="white",
            linewidth=0.6,
            label=role_labels[role] if role not in plotted_roles else None,
            zorder=3,
        )
        plotted_roles.add(role)
        axis.annotate(
            labels[str(row["ranking"])],
            (
                row["spearman_vs_ADRS_comp"],
                row["top20_jaccard"],
            ),
            xytext=label_offsets[str(row["ranking"])],
            textcoords="offset points",
            fontsize=8,
            color=DARK,
        )
    axis.set_xlim(-0.2, 1.02)
    axis.set_ylim(0.0, 1.02)
    axis.set_xlabel("Spearman correlation with ADRS_comp")
    axis.set_ylabel("Top-20 Jaccard overlap")
    axis.set_title("Descriptive component/baseline concordance")
    axis.legend(
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=0.92,
        fontsize=8,
        loc="upper left",
    )
    _panel_label(axis, "A")


def _plot_centrality(
    axis: plt.Axes,
    genes: Sequence[Mapping[str, Any]],
    metrics: Mapping[str, Any],
    kind: str,
) -> None:
    r_acc = np.asarray([row["r_ACC"] for row in genes], dtype=float)
    seeds = np.asarray([row["is_ACC_seed"] for row in genes], dtype=bool)
    if kind == "degree":
        x = np.log2(
            np.asarray([row["STRING_degree"] for row in genes], dtype=float) + 1
        )
        xlabel = "log2(STRING degree + 1)"
        record = metrics["centrality_correlations"]["gene_rACC_vs_degree"]
        title = "r_ACC strongly tracks STRING degree"
        panel = "B"
        extra = ""
    else:
        x = np.log10(
            np.asarray([row["STRING_PageRank"] for row in genes], dtype=float)
        )
        xlabel = "log10(STRING PageRank)"
        record = metrics["centrality_correlations"]["gene_rACC_vs_PageRank"]
        partial = metrics["centrality_correlations"][
            "gene_rACC_vs_PageRank_partial_log_degree"
        ]["rho"]
        title = "PageRank association vanishes\nafter degree control"
        panel = "C"
        extra = f"\npartial rho={partial:.3f}"
    axis.scatter(
        x[~seeds],
        r_acc[~seeds],
        s=15,
        color=LIGHT_BLUE,
        alpha=0.72,
        edgecolor="none",
        label="non-seed",
        rasterized=False,
    )
    axis.scatter(
        x[seeds],
        r_acc[seeds],
        s=25,
        facecolor=ORANGE,
        edgecolor=DARK,
        linewidth=0.35,
        label="ACC seed",
        zorder=3,
    )
    axis.set_xlabel(xlabel)
    axis.set_ylabel("r_ACC")
    axis.set_ylim(0.0, 1.03)
    axis.set_title(title)
    axis.text(
        0.04,
        0.96,
        f"Spearman rho={record['rho']:.3f}\n"
        f"95% bootstrap CI "
        f"[{record['bootstrap_ci_95'][0]:.3f}, "
        f"{record['bootstrap_ci_95'][1]:.3f}]"
        f"{extra}\nn=399",
        transform=axis.transAxes,
        va="top",
        fontsize=8,
        color=DARK,
    )
    if kind == "degree":
        axis.legend(frameon=False, fontsize=8, loc="lower right")
    _panel_label(axis, panel)


def _plot_null(
    axis: plt.Axes,
    rows: Sequence[Mapping[str, Any]],
    metrics: Mapping[str, Any],
) -> list[dict[str, Any]]:
    selected = sorted(
        rows,
        key=lambda row: (-float(row["z_degree_matched"]), str(row["drug"])),
    )[:12]
    selected = list(reversed(selected))
    y = np.arange(len(selected))
    means = np.asarray(
        [row["null_mean_C_ACC_pct"] for row in selected], dtype=float
    )
    standard_deviations = np.asarray(
        [row["null_sd_C_ACC_pct"] for row in selected], dtype=float
    )
    observed = np.asarray(
        [row["observed_C_ACC_pct"] for row in selected], dtype=float
    )
    axis.errorbar(
        means,
        y,
        xerr=standard_deviations,
        fmt="o",
        color=GREY,
        ecolor=GREY,
        markersize=4,
        linewidth=1.0,
        capsize=2,
        label="null mean ± 1 SD",
        zorder=2,
    )
    axis.scatter(
        observed,
        y,
        marker="D",
        s=30,
        color=BLUE,
        edgecolor="white",
        linewidth=0.4,
        label="observed",
        zorder=3,
    )
    axis.set_yticks(y, [str(row["drug"]) for row in selected])
    axis.set_xlim(0.0, 1.05)
    axis.set_xlabel("C_ACC percentile")
    significant_n = sum(float(row["q_bh_108"]) < 0.05 for row in rows)
    axis.set_title(
        "Degree-matched random-seed null\n"
        f"{significant_n}/108 drugs q<0.05; CDK4/6 P="
        f"{metrics['CDK46_degree_matched_empirical_p']:.4f}"
    )
    _panel_label(axis, "D")
    return [dict(row) for row in selected]


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _upsert_manifest(
    path: Path,
    project_root: Path,
    figure_paths: Mapping[str, Path],
    qa: Mapping[str, Any],
) -> None:
    caption = (
        "Figure 2. Reviewer-priority methodological strengthening of the locked "
        "108-drug ranking. (A) Descriptive agreement of components and simple "
        "baselines with ADRS_comp; this panel is not an inferential "
        "non-redundancy test. "
        "(B,C) Gene-level r_ACC was strongly associated with STRING degree and "
        "PageRank, while the PageRank partial correlation after degree control "
        "was near zero. (D) The twelve largest positive drug-level deviations "
        "from 10,000 degree-matched random-seed nulls; points show observed "
        "C_ACC percentiles and bars show null mean ±1 SD. No drug passed BH "
        "FDR for pralatrexate and pemetrexed, while the CDK4/6 group-level "
        "empirical P was not significant."
    )
    caption = (
        "Figure 2. Reviewer-priority methodological audit of the locked "
        "108-drug ranking. (A) Descriptive agreement of components and simple "
        "baselines with ADRS_comp; this is not an inferential non-redundancy "
        "test. (B,C) Gene-level r_ACC was strongly associated with STRING "
        "degree and PageRank, while the PageRank partial correlation after "
        "degree control was near zero. (D) The twelve largest positive "
        "drug-level deviations from 10,000 degree-matched random-seed nulls; "
        "points show observed C_ACC percentiles and bars show null mean ±1 SD. "
        "Pralatrexate and pemetrexed passed BH FDR, while the CDK4/6 "
        "group-level empirical P was not significant."
    )
    section = "\n".join(
        [
            "## F2",
            "",
            "- figure_id: `F2`",
            "- source_card: `new_canonical_candidate -> "
            "databases/db07-figures/resources_real.md`",
            f"- vector PDF: `{figure_paths['pdf'].relative_to(project_root).as_posix()}`",
            f"- vector SVG: `{figure_paths['svg'].relative_to(project_root).as_posix()}`",
            f"- bitmap PNG: `{figure_paths['png'].relative_to(project_root).as_posix()}`",
            "- section: Results methodological robustness",
            f"- target: `{TARGET_JOURNAL}`, `{TARGET_COLUMN}`, "
            f"{TARGET_WIDTH_MM:.0f} mm",
            f"- checks: width={'pass' if qa['png']['width_check'] else 'fail'}; "
            f"dpi={'pass' if qa['png']['dpi_check'] else 'fail'}; "
            f"font={'pass' if qa['figure']['font_check'] else 'fail'}",
            f"- caption: {caption}",
            "",
        ]
    )
    current = (
        path.read_text(encoding="utf-8")
        if path.is_file()
        else "# ACC-PHARMA-NET figure manifest\n"
    )
    start = current.find("## F2\n")
    if start >= 0:
        next_section = current.find("\n## ", start + 1)
        current = (
            current[:start]
            + section
            + (current[next_section + 1 :] if next_section >= 0 else "")
        )
    else:
        current = current.rstrip() + "\n\n" + section
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(current.rstrip() + "\n", encoding="utf-8")


def generate_method_strengthening_figure(
    project_root: Path,
) -> dict[str, Any]:
    root = project_root.resolve()
    baselines, genes, null_rows, metrics = load_figure_inputs(root)
    _set_publication_style()
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(TARGET_WIDTH_MM / 25.4, 148.0 / 25.4),
        layout="constrained",
    )
    _plot_baselines(axes[0, 0], baselines)
    _plot_centrality(axes[0, 1], genes, metrics, "degree")
    _plot_centrality(axes[1, 0], genes, metrics, "pagerank")
    selected_null = _plot_null(axes[1, 1], null_rows, metrics)
    fig.suptitle(
        "Methodological strengthening of the locked ACC drug ranking",
        fontsize=10,
        fontweight="bold",
    )
    figure_dir = root / "figures" / "revision"
    figure_dir.mkdir(parents=True, exist_ok=True)
    basename = figure_dir / "Fig2_method_strengthening_primary108"
    paths = _export_figure(fig, basename)
    minimum_font = _minimum_visible_font(fig)
    figure_qa = {
        "figure_version": FIGURE_VERSION,
        "target_width_mm": TARGET_WIDTH_MM,
        "minimum_visible_font_pt": minimum_font,
        "font_check": minimum_font >= MIN_FONT_PT,
    }
    plt.close(fig)
    png_qa = _inspect_png(paths["png"])
    qa = {
        "target_journal": TARGET_JOURNAL,
        "target_column": TARGET_COLUMN,
        "minimum_font_pt": MIN_FONT_PT,
        "png_dpi": PNG_DPI,
        "figure": figure_qa,
        "png": png_qa,
    }
    if not figure_qa["font_check"]:
        raise ValueError(f"Figure 2 contains text below {MIN_FONT_PT} pt")
    if not png_qa["dpi_check"] or not png_qa["width_check"]:
        raise ValueError("Figure 2 PNG failed DPI or physical-width checks")

    data_path = (
        root
        / "figure_data"
        / "revision"
        / "Fig2d_random_seed_top12_primary108.csv"
    )
    _write_csv(data_path, selected_null)
    qa_path = root / "results" / "figures" / "MS_F2_figure_QA.json"
    qa_path.parent.mkdir(parents=True, exist_ok=True)
    qa_path.write_text(
        json.dumps(qa, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    manifest_path = (
        root / "projects" / "ACC-PHARMA-NET" / "figures" / "manifest.md"
    )
    _upsert_manifest(manifest_path, root, paths, qa)
    return {
        "paths": paths,
        "figure_data": data_path,
        "qa": qa,
        "qa_path": qa_path,
        "manifest": manifest_path,
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
    output = generate_method_strengthening_figure(args.project_root)
    print(
        json.dumps(
            {
                "status": "ok",
                "figure_version": FIGURE_VERSION,
                "pdf": str(output["paths"]["pdf"]),
                "font_check": output["qa"]["figure"]["font_check"],
                "dpi_check": output["qa"]["png"]["dpi_check"],
                "width_check": output["qa"]["png"]["width_check"],
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
