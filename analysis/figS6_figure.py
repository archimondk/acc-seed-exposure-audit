"""Build Supplementary Figure S6: same-resource NCI-CCR/NCATS concordance (A-C) and descriptive ACC biomarker panel (D)."""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


FIGURE_VERSION = "figS6-crossplatform-biomarker-v1"
PNG_DPI = 600
WIDTH_IN = 11.0
HEIGHT_IN = 7.4

# Deterministic-rendering controls. SVG path ids are salted (otherwise
# matplotlib emits random UUIDs per render) and PDF/SVG creation metadata is
# pinned to a fixed timestamp so that byte-identical re-renders are possible.
SVG_HASHSALT = "figs6-crossplatform-biomarker-v1"
RENDER_TIMESTAMP = "2026-08-02T00:00:00+08:00"
FIXED_DATETIME = datetime.datetime(2026, 8, 2, 0, 0, 0)
PDF_METADATA = {
    "CreationDate": FIXED_DATETIME,
    "ModDate": FIXED_DATETIME,
    "Producer": "matplotlib",
    "Creator": "matplotlib",
}
SVG_METADATA = {"Date": FIXED_DATETIME, "Creator": "matplotlib"}

# Frozen manuscript statistics (Section 3.9 / legend S6). These values are
# annotated on the panels and must match the manuscript exactly.
STATS = {
    "CU-ACC1": {"n": 21, "rho": 0.84, "ci": (0.65, 0.92), "p": "<0.001"},
    "CU-ACC2": {"n": 20, "rho": 0.58, "ci": (0.08, 0.90), "p": "0.008"},
    "NCI-H295R": {"n": 16, "rho": 0.61, "ci": (0.14, 0.88), "p": "0.013"},
}

# Panel D genes (legend S6: steroidogenic markers, CDK4, RB1) and the six
# ACC surgical tumours.
D_GENES = ["IGF2", "NR5A1", "CYP11A1", "STAR", "SOAT1", "CDK4", "RB1"]
D_TUMORS = ["30_SC", "31_SC", "33_SC", "34_SC", "37_SC", "38_SC"]

BOOTSTRAP_SEED = 42


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_inputs(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    abc = pd.read_csv(root / "figure_data" / "Fig6abc_crossplatform_concordance.csv")
    biom = pd.read_csv(root / "figure_data" / "Fig6d_biomarker_expression.csv")

    # Structural validation against the frozen manuscript.
    counts = abc.groupby("cell_line").size().to_dict()
    for cl, meta in STATS.items():
        if counts.get(cl) != meta["n"]:
            raise ValueError(
                f"Expected {meta['n']} rows for {cl} in Fig6abc, "
                f"got {counts.get(cl)}"
            )
    missing_genes = [g for g in D_GENES if g not in set(biom["gene"])]
    if missing_genes:
        raise ValueError(f"Missing panel-D genes in Fig6d: {missing_genes}")
    missing_tumors = [t for t in D_TUMORS if t not in biom.columns]
    if missing_tumors:
        raise ValueError(f"Missing panel-D tumours in Fig6d: {missing_tumors}")
    return abc, biom


def _render(abc: pd.DataFrame, biom: pd.DataFrame, out_stem: Path) -> None:
    fig = plt.figure(figsize=(WIDTH_IN, HEIGHT_IN))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1.05],
                          hspace=0.38, wspace=0.32)

    # Panels A-C: NCI-CCR vs NCATS/MIPE activity, one ACC line each.
    for i, (cl, meta) in enumerate(STATS.items()):
        ax = fig.add_subplot(gs[0, i])
        sub = abc[abc["cell_line"] == cl].dropna()
        x = sub["NCICCR_neglogIC50"].values
        y = sub["NCATS_neglogIC50"].values
        ax.scatter(x, y, s=22, alpha=0.85, edgecolors="k", linewidths=0.3,
                   color="#1f77b4", zorder=3)
        lo = min(x.min(), y.min())
        hi = max(x.max(), y.max())
        ax.plot([lo, hi], [lo, hi], "--", color="grey", lw=0.8, zorder=1)
        ax.set_xlabel("NCI-CCR activity\n($-\\log_{10}$ GI$_{50}$)", fontsize=8.5)
        ax.set_ylabel("NCATS/MIPE activity\n(vAUC Z-score)", fontsize=8.5)
        ax.tick_params(labelsize=7.5)
        ax.set_title(f"({chr(65 + i)}) {cl}", loc="left", fontsize=10,
                     fontweight="bold")
        txt = (f"n = {meta['n']}\n"
               f"$\\rho$ = {meta['rho']:.2f}\n"
               f"95% CI {meta['ci'][0]:.2f}\u2013{meta['ci'][1]:.2f}\n"
               f"P {meta['p']}")
        ax.text(0.97, 0.05, txt, transform=ax.transAxes, fontsize=7.6,
                va="bottom", ha="right",
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.6", lw=0.6))

    # Panel D: descriptive log2(FPKM+1) heatmap across the six tumours.
    axd = fig.add_subplot(gs[1, :])
    dmat = biom.set_index("gene").reindex(D_GENES)[D_TUMORS].astype(float)
    vmax = float(dmat.values.max())
    im = axd.imshow(dmat.values, aspect="auto", cmap="YlOrRd",
                    vmin=0, vmax=vmax, interpolation="nearest")
    axd.set_xticks(range(len(D_TUMORS)))
    axd.set_xticklabels(D_TUMORS, fontsize=8.5)
    axd.set_yticks(range(len(D_GENES)))
    axd.set_yticklabels(D_GENES, fontsize=8.5)
    axd.set_title("(D) ACC surgical tumours \u2014 log$_2$(FPKM+1) expression",
                  loc="left", fontsize=10, fontweight="bold")
    for r in range(len(D_GENES)):
        for c in range(len(D_TUMORS)):
            v = dmat.values[r, c]
            axd.text(c, r, f"{v:.1f}", ha="center", va="center", fontsize=6.6,
                     color="black" if v < 0.62 * vmax else "white")
    cbar = fig.colorbar(im, ax=axd, fraction=0.035, pad=0.015)
    cbar.ax.tick_params(labelsize=7.5)
    cbar.set_label("log$_2$(FPKM+1)", fontsize=8.5)
    axd.text(0.0, -0.28,
             "Descriptive biomarker observations only; not an independent "
             "drug-response endpoint.",
             transform=axd.transAxes, fontsize=7.6, style="italic", color="0.25")

    fig.savefig(out_stem.with_suffix(".png"), dpi=PNG_DPI, bbox_inches="tight")
    fig.savefig(out_stem.with_suffix(".pdf"), bbox_inches="tight",
                metadata=PDF_METADATA)
    fig.savefig(out_stem.with_suffix(".svg"), bbox_inches="tight",
                metadata=SVG_METADATA)
    plt.close(fig)


def build_figure(root: Path) -> dict[str, Any]:
    root = root.resolve()
    # Deterministic SVG ids (matplotlib otherwise emits random UUIDs).
    plt.rcParams["svg.hashsalt"] = SVG_HASHSALT
    abc, biom = _load_inputs(root)
    out_dir = root / "figures" / "revision"
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / "FigS6_crossplatform_biomarker_ACCcmdb"
    _render(abc, biom, stem)

    manifest = {
        "figure_version": FIGURE_VERSION,
        "rendered": "frozen-input deterministic render",
        "target_journal": "mdpi",
        "statistics_annotated": STATS,
        "ci_method": "bootstrap percentile, 2000 drug resamples (verification only)",
        "bootstrap_rng_seed": BOOTSTRAP_SEED,
        "panel_d_genes": D_GENES,
        "panel_d_tumors": D_TUMORS,
        "inputs": {
            "figure_data/Fig6abc_crossplatform_concordance.csv": _sha256(
                root / "figure_data" / "Fig6abc_crossplatform_concordance.csv"
            ),
            "figure_data/Fig6d_biomarker_expression.csv": _sha256(
                root / "figure_data" / "Fig6d_biomarker_expression.csv"
            ),
        },
        "outputs": {
            f"figures/revision/FigS6_crossplatform_biomarker_ACCcmdb.{ext}": _sha256(
                root / "figures" / "revision" / f"FigS6_crossplatform_biomarker_ACCcmdb.{ext}"
            )
            for ext in ("png", "pdf", "svg")
        },
    }
    manifest_path = out_dir / "FigS6_crossplatform_biomarker_ACCcmdb_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return {
        "outputs": {
            ext: str(root / "figures" / "revision"
                     / f"FigS6_crossplatform_biomarker_ACCcmdb.{ext}")
            for ext in ("png", "pdf", "svg")
        },
        "manifest": str(manifest_path),
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
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
