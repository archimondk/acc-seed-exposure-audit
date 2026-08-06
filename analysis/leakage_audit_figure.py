"""Build the frozen seed–target leakage audit figure."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd


FIGURE_VERSION = "leakage-audit-figure-v1"
WIDTH_MM = 170.0
HEIGHT_MM = 78.0
PNG_DPI = 1000

BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
VERMILION = "#D55E00"
GREY = "#A7A9AC"
DARK = "#202124"
LIGHT_GREY = "#E8EAED"

VARIANT_LABELS = {
    "column_minmax": "Primary\nmin–max",
    "column_gene_rank": "Column\nrank",
    "uniform_ratio_gene_rank": "Uniform-adjusted\nrank",
    "symmetric_gene_rank": "Symmetric\nrank",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_inputs(root: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    observational = pd.read_csv(
        root / "leakage_audit" / "observational" / "Ld_108.csv"
    )
    curve = pd.read_csv(root / "leakage_audit" / "leakage_curve_108.csv")
    verdict = json.loads(
        (root / "leakage_audit" / "verdict.json").read_text(encoding="utf-8")
    )
    if len(observational) != 108:
        raise ValueError(f"Expected 108 observational drug rows, got {len(observational)}")
    expected_curve = 108 * 4
    if len(curve) != expected_curve:
        raise ValueError(f"Expected {expected_curve} curve rows, got {len(curve)}")
    if verdict["status"] != "PARTIAL_OR_NOT_SUPPORTED":
        raise ValueError(f"Unexpected frozen verdict: {verdict['status']}")
    return observational, curve, verdict


def _write_figure_data(
    root: Path,
    observational: pd.DataFrame,
    curve: pd.DataFrame,
) -> dict[str, Path]:
    out = root / "figure_data" / "revision"
    out.mkdir(parents=True, exist_ok=True)
    panel_a = out / "Fig7a_seed_overlap_primary108.csv"
    panel_b = out / "Fig7b_leakage_effects_variants.csv"
    panel_c = out / "Fig7c_intervention_delta_primary108.csv"

    observational[
        [
            "drug",
            "n_d",
            "n_seed",
            "seed_hits",
            "L_weighted",
            "obs_C_ACC_pct",
        ]
    ].to_csv(panel_a, index=False)
    abemaciclib = curve[curve["drug"].str.casefold() == "abemaciclib"].copy()
    abemaciclib[
        [
            "variant",
            "delta_z_acc_A1_minus_A2",
            "delta_z_breast_B2_minus_B1",
            "delta_z_breast_B2_lo_minus_B1",
            "delta_z_breast_B2_hi_minus_B1",
        ]
    ].to_csv(panel_b, index=False)
    primary = curve[curve["variant"] == "column_minmax"].copy()
    primary[
        [
            "drug",
            "delta_z_acc_A1_minus_A2",
            "delta_z_breast_B2_minus_B1",
        ]
    ].to_csv(panel_c, index=False)
    return {"panel_a": panel_a, "panel_b": panel_b, "panel_c": panel_c}


def _panel_a(ax: plt.Axes, observational: pd.DataFrame) -> None:
    ranked = observational.sort_values(
        ["obs_C_ACC_pct", "drug"], ascending=[False, True]
    ).copy()
    top20 = set(ranked.head(20)["drug"])
    seed_mask = observational["n_seed"] > 0

    ax.scatter(
        observational.loc[~seed_mask, "L_weighted"],
        observational.loc[~seed_mask, "obs_C_ACC_pct"],
        s=14,
        marker="o",
        facecolors="white",
        edgecolors=GREY,
        linewidths=0.7,
        label="No direct seed association",
        zorder=2,
    )
    ax.scatter(
        observational.loc[seed_mask, "L_weighted"],
        observational.loc[seed_mask, "obs_C_ACC_pct"],
        s=17,
        marker="o",
        color=BLUE,
        alpha=0.78,
        linewidths=0,
        label="≥1 seed association",
        zorder=3,
    )
    top_rows = observational[observational["drug"].isin(top20)]
    ax.scatter(
        top_rows["L_weighted"],
        top_rows["obs_C_ACC_pct"],
        s=34,
        marker="o",
        facecolors="none",
        edgecolors=ORANGE,
        linewidths=1.0,
        label="C_ACC Top 20",
        zorder=4,
    )
    ax.axhline(0.75, color=GREY, linestyle="--", linewidth=0.8, zorder=1)
    ax.text(
        0.002,
        0.765,
        "Top-quartile boundary",
        fontsize=7,
        color="#5F6368",
        va="bottom",
    )
    for drug in ("abemaciclib", "palbociclib"):
        row = observational[observational["drug"].str.casefold() == drug].iloc[0]
        ax.annotate(
            drug,
            (row["L_weighted"], row["obs_C_ACC_pct"]),
            xytext=(3, 3),
            textcoords="offset points",
            fontsize=6.5,
            color=DARK,
        )
    ax.text(
        0.03,
        0.97,
        "ρ = 0.828\nTop-20 Jaccard = 1.000\nNo-seed Top 25%: 0/62",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.8,
        bbox={"facecolor": "white", "edgecolor": LIGHT_GREY, "pad": 2.5},
    )
    ax.set_xlabel("Shrinkage-weighted seed-overlap share")
    ax.set_ylabel("Observed C_ACC percentile")
    ax.set_xlim(left=-0.0005)
    ax.set_ylim(0, 1.03)
    ax.legend(
        loc="lower left",
        fontsize=5.9,
        frameon=False,
        handletextpad=0.3,
        labelspacing=0.2,
    )


def _panel_b(ax: plt.Axes, curve: pd.DataFrame, verdict: dict) -> None:
    rows = (
        curve[curve["drug"].str.casefold() == "abemaciclib"]
        .set_index("variant")
        .loc[list(VARIANT_LABELS)]
    )
    y = list(range(len(rows)))
    acc = rows["delta_z_acc_A1_minus_A2"].to_numpy()
    breast = rows["delta_z_breast_B2_minus_B1"].to_numpy()

    for yi, left, right in zip(y, acc, breast, strict=True):
        ax.plot([left, right], [yi, yi], color=GREY, linewidth=1.0, zorder=1)
    ax.scatter(
        acc,
        y,
        s=32,
        marker="o",
        color=BLUE,
        edgecolors="white",
        linewidths=0.5,
        label="ACC: A1 − A2",
        zorder=3,
    )
    ax.scatter(
        breast,
        y,
        s=34,
        marker="s",
        color=ORANGE,
        edgecolors="white",
        linewidths=0.5,
        label="Breast: B2 − B1",
        zorder=3,
    )
    l1 = float(verdict["criteria"]["L1"]["threshold"])
    l2 = float(verdict["criteria"]["L2"]["threshold"])
    ax.axvline(l1, color=BLUE, linestyle="--", linewidth=0.9)
    ax.axvline(l2, color=ORANGE, linestyle=":", linewidth=1.1)
    ax.text(
        l1 + 0.04,
        0.01,
        "L1 = 2.0",
        transform=ax.get_xaxis_transform(),
        color=BLUE,
        fontsize=6.3,
        va="bottom",
    )
    ax.text(
        l2 - 0.04,
        0.01,
        "L2 = 1.5",
        transform=ax.get_xaxis_transform(),
        color="#9A6700",
        fontsize=6.3,
        ha="right",
        va="bottom",
    )
    ax.set_yticks(y, [VARIANT_LABELS[item] for item in rows.index])
    ax.invert_yaxis()
    ax.set_xlabel("Abemaciclib intervention effect, Δz")
    ax.set_xlim(0, 4.25)
    ax.text(
        0.97,
        0.19,
        "L4: 1/4 variants\n(required ≥3/4)",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.8,
        fontweight="bold",
        color=VERMILION,
        bbox={"facecolor": "white", "edgecolor": VERMILION, "pad": 2.5},
    )
    ax.annotate(
        "ACC",
        (acc[0], y[0]),
        xytext=(0, 7),
        textcoords="offset points",
        color=BLUE,
        fontsize=6.4,
        ha="center",
        fontweight="bold",
    )
    ax.annotate(
        "Breast",
        (breast[0], y[0]),
        xytext=(0, 7),
        textcoords="offset points",
        color="#9A6700",
        fontsize=6.4,
        ha="center",
        fontweight="bold",
    )


def _panel_c(ax: plt.Axes, curve: pd.DataFrame) -> None:
    primary = curve[curve["variant"] == "column_minmax"].copy()
    x_name = "delta_z_acc_A1_minus_A2"
    y_name = "delta_z_breast_B2_minus_B1"
    focal = {"abemaciclib", "palbociclib", "ribociclib"}
    base = ~primary["drug"].str.casefold().isin(focal)

    ax.add_patch(
        Rectangle(
            (-0.5, -0.5),
            1.0,
            1.0,
            facecolor=LIGHT_GREY,
            edgecolor=GREY,
            linewidth=0.8,
            linestyle="--",
            zorder=0,
        )
    )
    ax.scatter(
        primary.loc[base, x_name],
        primary.loc[base, y_name],
        s=14,
        marker="o",
        facecolors="white",
        edgecolors=GREY,
        linewidths=0.65,
        label="Other drugs",
        zorder=2,
    )
    styles = {
        "abemaciclib": (VERMILION, "D"),
        "palbociclib": (ORANGE, "s"),
        "ribociclib": (GREEN, "^"),
    }
    for drug, (color, marker) in styles.items():
        row = primary[primary["drug"].str.casefold() == drug].iloc[0]
        ax.scatter(
            [row[x_name]],
            [row[y_name]],
            s=42,
            marker=marker,
            color=color,
            edgecolors="white",
            linewidths=0.6,
            zorder=4,
        )
        if drug == "abemaciclib":
            offset = (-5, 5)
            horizontal_alignment = "right"
        elif drug == "palbociclib":
            offset = (4, 12)
            horizontal_alignment = "left"
        else:
            offset = (4, -12)
            horizontal_alignment = "left"
        ax.annotate(
            drug,
            (row[x_name], row[y_name]),
            xytext=offset,
            textcoords="offset points",
            fontsize=6.8,
            color=DARK,
            ha=horizontal_alignment,
        )
    ax.axhline(0, color="#5F6368", linewidth=0.6)
    ax.axvline(0, color="#5F6368", linewidth=0.6)
    ax.text(
        0.97,
        0.05,
        "NC2: 106/106 pass\nmax |Δz|: 0.319, 0.367",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.8,
        bbox={"facecolor": "white", "edgecolor": LIGHT_GREY, "pad": 2.5},
    )
    ax.set_xlabel("ACC intervention effect, Δz")
    ax.set_ylabel("Breast intervention effect, Δz")
    ax.set_xlim(-0.55, 4.25)
    ax.set_ylim(-0.55, 2.2)


def build_figure(project_root: Path) -> dict[str, object]:
    root = project_root.resolve()
    observational, curve, verdict = _load_inputs(root)
    figure_data = _write_figure_data(root, observational, curve)

    matplotlib.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 8,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.03, 1.0, 1.03]},
    )
    _panel_a(axes[0], observational)
    _panel_b(axes[1], curve, verdict)
    _panel_c(axes[2], curve)
    for label, ax in zip(("A", "B", "C"), axes, strict=True):
        ax.text(
            -0.18,
            1.06,
            label,
            transform=ax.transAxes,
            fontsize=10,
            fontweight="bold",
            va="top",
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(width=0.7, length=3)

    out_dir = root / "figures" / "revision"
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / "Fig5_seed_target_leakage_audit"
    outputs = {
        "pdf": stem.with_suffix(".pdf"),
        "svg": stem.with_suffix(".svg"),
        "png": stem.with_suffix(".png"),
    }
    fig.savefig(outputs["pdf"])
    fig.savefig(outputs["svg"])
    fig.savefig(outputs["png"], dpi=PNG_DPI)
    plt.close(fig)

    manifest = {
        "figure_version": FIGURE_VERSION,
        "target_journal": "mdpi",
        "target_column": "full",
        "target_width_mm": WIDTH_MM,
        "target_height_mm": HEIGHT_MM,
        "png_dpi": PNG_DPI,
        "inputs": {
            "leakage_audit/observational/Ld_108.csv": _sha256(
                root / "leakage_audit" / "observational" / "Ld_108.csv"
            ),
            "leakage_audit/leakage_curve_108.csv": _sha256(
                root / "leakage_audit" / "leakage_curve_108.csv"
            ),
            "leakage_audit/verdict.json": _sha256(
                root / "leakage_audit" / "verdict.json"
            ),
        },
        "figure_data": {
            str(path.relative_to(root)).replace("\\", "/"): _sha256(path)
            for path in figure_data.values()
        },
        "outputs": {
            str(path.relative_to(root)).replace("\\", "/"): _sha256(path)
            for path in outputs.values()
        },
    }
    manifest_path = out_dir / "Fig5_seed_target_leakage_audit_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return {
        "outputs": outputs,
        "manifest": manifest_path,
        "figure_data": figure_data,
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
    result = build_figure(args.project_root)
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
