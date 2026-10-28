"""Fig 1 (revision R1): study workflow with the melanoma arm and five status labels.

Static publication figure for PLOS ONE: 7.5 in wide, text 8-9 pt, RGB TIFF
(LZW) at 600 dpi plus PDF and PNG previews. Status assignments follow
manuscript Table 1.

Usage: python -I fig1_workflow_r1.py <out_dir>
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 8, "pdf.fonttype": 42,
                     "savefig.facecolor": "white"})

INK, MUTED, LINE, PANEL = "#202124", "#5F6368", "#7A858E", "#F7F8F9"
STATUS = {  # Okabe-Ito based; validated with dataviz validate_palette.js (light)
    "original": ("#0072B2", "ORIGINAL", "white"),
    "exploratory": ("#E69F00", "EXPLORATORY", INK),
    "frozen": ("#009E73", "FROZEN BEFORE RESULT", "white"),
    "posthoc": ("#8B5FBF", "RESULT-KNOWN POST HOC", "white"),
    "revision": ("#D55E00", "REVISION · 2-STAGE FREEZE", "white"),
}
W, H = 7.5, 8.75


def box(ax, x, y, w, h, edge, fill="white", lw=0.9, ls="-", r=0.06, z=1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                linewidth=lw, edgecolor=edge, facecolor=fill, linestyle=ls, zorder=z))


def arrow(ax, p0, p1, color=LINE, lw=0.9):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=7, linewidth=lw,
                                 color=color, shrinkA=0, shrinkB=0, zorder=0.5))


def section(ax, y, text):
    ax.text(0.12, y, text, ha="left", va="center", fontsize=8.5, fontweight="bold", color=MUTED, zorder=7,
            bbox=dict(boxstyle="square,pad=0.15", facecolor="white", edgecolor="none"))
    ax.plot([0.12 + 0.075 * len(text), W - 0.1], [y, y], color="#DADCE0", lw=0.6, zorder=0)


def badge(ax, x, y, status, fontsize=7.6):
    color, label, txt = STATUS[status]
    t = ax.text(x, y, label, ha="left", va="center", fontsize=fontsize, fontweight="bold", color=txt, zorder=6,
                bbox=dict(boxstyle="round,pad=0.22,rounding_size=0.35", facecolor=color, edgecolor=color))
    return t


def input_box(ax, x, y, w, h, title, lines, accent):
    box(ax, x, y, w, h, "#C9CED3", fill=PANEL)
    ax.add_patch(plt.Rectangle((x, y), 0.045, h, color=accent, zorder=2, linewidth=0))
    ax.text(x + 0.11, y + h - 0.12, title, ha="left", va="top", fontsize=8.5, fontweight="bold", color=INK)
    for i, line in enumerate(lines):
        ax.text(x + 0.11, y + h - 0.40 - 0.17 * i, line, ha="left", va="top", fontsize=8, color=INK)


def step_box(ax, x, y, w, h, title, sub, status="original", strong=False):
    color = STATUS[status][0]
    box(ax, x, y, w, h, color if strong else "#AEB6BD", lw=1.3 if strong else 0.9)
    badge(ax, x + 0.09, y + h - 0.15, status, fontsize=8)
    ax.text(x + 0.09, y + h - 0.36, title, ha="left", va="top", fontsize=8.5, fontweight="bold", color=INK,
            linespacing=1.15)
    ax.text(x + 0.09, y + 0.09, sub, ha="left", va="bottom", fontsize=8, color=MUTED, linespacing=1.15)


def audit_box(ax, x, y, w, h, status, title, lines):
    color, label, txt = STATUS[status]
    box(ax, x, y, w, h, color, lw=1.0)
    ax.add_patch(FancyBboxPatch((x, y + h - 0.27), w, 0.27, boxstyle="round,pad=0,rounding_size=0.06",
                                facecolor=color, edgecolor=color, linewidth=1.0, zorder=2))
    ax.add_patch(plt.Rectangle((x, y + h - 0.27), w, 0.12, facecolor=color, edgecolor=color, zorder=2))
    ax.text(x + 0.08, y + h - 0.135, label, ha="left", va="center", fontsize=8, fontweight="bold", color=txt,
            zorder=3)
    ax.text(x + 0.08, y + h - 0.37, title, ha="left", va="top", fontsize=8.5, fontweight="bold", color=INK,
            linespacing=1.12)
    for i, line in enumerate(lines):
        yy = y + h - 0.68 - 0.21 * i
        ax.plot(x + 0.11, yy - 0.055, marker="o", ms=2.6, color=color, zorder=3)
        ax.text(x + 0.18, yy, line, ha="left", va="top", fontsize=8, color=INK)


def draw(out: Path):
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")

    # ---- data inputs
    section(ax, 8.64, "DATA INPUTS")
    xs = [0.1 + i * (1.40 + 0.075) for i in range(5)]
    cw = 1.40
    inputs = [
        ("Association map", ["124 drugs, 399 genes", "1304 associations"], "#0072B2"),
        ("ACC module", ["50 genes → 45 seeds", "G/R/P/L/S weights"], "#009E73"),
        ("STRING v12", ["Combined score ≥ 400", "929,472 edges"], "#5F6368"),
        ("TCGA (cBioPortal)", ["ER+/HER2− breast", "BRAF-mutant melanoma"], "#8C5A3C"),
        ("Drug activity", ["MIPE 5.0, 3 ACC lines", "NCI-60 potency"], "#4E7D3A"),
    ]
    for x, (title, lines, accent) in zip(xs, inputs):
        input_box(ax, x, 7.62, cw, 0.92, title, lines, accent)

    # ---- core workflow
    section(ax, 7.36, "CORE ACC WORKFLOW")
    ch, cy = 1.02, 6.12
    core = [
        ("Network QC", "Exact edge and\nidentifier audit"),
        ("45 weighted\nseeds", "44 propagate;\nVEGFA unmapped"),
        ("Propagation\n→ r$_{ACC}$", "RWR, α = 0.40;\nmin–max, 399 genes"),
        ("Drug context\nscore C$_{ACC}$", "Shrunken gene-set\nmean (k = 3)"),
        ("Activity\nresidual", "MIPE on NCI-60\npotency"),
    ]
    for i, (x, (title, sub)) in enumerate(zip(xs, core)):
        step_box(ax, x, cy, cw, ch, title, sub)
        if i < 3:
            arrow(ax, (x + cw, cy + ch / 2), (xs[i + 1], cy + ch / 2))
    for i in (0, 1, 2, 4):
        arrow(ax, (xs[i] + cw / 2, 7.62), (xs[i] + cw / 2, cy + ch))
    # universes note and locked composite
    step_box(ax, xs[0], 4.92, xs[2] + cw - xs[0], 0.95, "Drug universes",
             "124-drug network universe (melanoma primary)\n108 complete cases for ADRS; 16 lack MIPE data")
    step_box(ax, xs[3], 4.92, xs[4] + cw - xs[3], 0.95, "ADRS$_{comp}$ (locked primary)",
             "0.50·P(C$_{ACC}$) + 0.50·P(residual)", strong=True)
    arrow(ax, (xs[3] + cw / 2, cy), (xs[3] + cw / 2, 5.87))
    arrow(ax, (xs[4] + cw / 2, cy), (xs[4] + cw / 2, 5.87))

    # ---- audits
    section(ax, 4.62, "AUDIT BRANCHES AND ANALYSIS STATUS")
    aw = (W - 0.2 - 3 * 0.1) / 4
    ax_ = [0.1 + i * (aw + 0.1) for i in range(4)]
    ah = 1.48
    row1 = [
        ("posthoc", "Topology and null audits", ["Degree, PageRank vs r$_{ACC}$", "10,000 matched-seed nulls", "3 rank-based normalizations"]),
        ("exploratory", "Seed-exposure diagnosis", ["Count / weight overlap", "Seed adjacency (O1–O4)", "Partial correlation"]),
        ("frozen", "Breast-cancer control", ["Outcome-blind TCGA seeds", "CDK4/6 recovery gate", "Gate failed (0 of 4)"]),
        ("frozen", "RB1 intervention", ["ACC −RB1; breast +RB1", "Nulls rebuilt per arm", "Rules met in 1 of 4"]),
    ]
    row2 = [
        ("posthoc", "Amendments 1–5", ["Scale-free re-analysis", "Seed-weight perturbation", "45-seed scan; seed exclusion"]),
        ("original", "Evidence boundary", ["Literature evidence types", "Cross-platform concordance", "No negative comparator"]),
        ("revision", "Second disease: melanoma", ["Stage A: cohort, seeds", "Stage B: exposure map", "H4–H6; BRAF restored"]),
        ("posthoc", "Revision diagnostics", ["Three-disease exposure", "Set-size effect", "Study-wide multiplicity"]),
    ]
    y1, y2 = 2.86, 1.16
    ax.plot([ax_[0] + aw / 2, ax_[3] + aw / 2], [4.47, 4.47], color=LINE, lw=0.9, zorder=0)
    arrow(ax, (xs[3] + cw + 0.04, 4.92), (xs[3] + cw + 0.04, 4.47))
    for x in ax_:
        arrow(ax, (x + aw / 2, 4.47), (x + aw / 2, y1 + ah))
    for x, (status, title, lines) in zip(ax_, row1):
        audit_box(ax, x, y1, aw, ah, status, title, lines)
    ax.plot([0.05, ax_[3] + aw / 2], [2.76, 2.76], color=LINE, lw=0.9, zorder=0)
    ax.plot([0.05, ax_[0] + aw / 2], [4.47, 4.47], color=LINE, lw=0.9, zorder=0)
    ax.plot([0.05, 0.05], [2.76, 4.47], color=LINE, lw=0.9, zorder=0)
    for x in ax_:
        arrow(ax, (x + aw / 2, 2.76), (x + aw / 2, y2 + ah))
    for x, (status, title, lines) in zip(ax_, row2):
        audit_box(ax, x, y2, aw, ah, status, title, lines)

    # ---- interpretation boundary and status key
    box(ax, 0.1, 0.40, W - 0.2, 0.56, "#5F6368", fill=PANEL, lw=0.9)
    ax.text(0.22, 0.68, "Interpretation boundary", ha="left", va="center", fontsize=8.5, fontweight="bold",
            color=MUTED)
    ax.text(1.85, 0.77, "Report locked and seed-excluded rankings together; rankings are audit outputs,",
            ha="left", va="center", fontsize=8.5, fontweight="bold", color=INK)
    ax.text(1.85, 0.58, "not validated efficacy predictions.", ha="left", va="center", fontsize=8.5,
            fontweight="bold", color=INK)
    x = 0.12
    ax.text(x, 0.17, "Status:", ha="left", va="center", fontsize=8, fontweight="bold", color=MUTED)
    x += 0.48
    for key, label in (("original", "Original"), ("exploratory", "Exploratory"), ("frozen", "Frozen before result"),
                       ("posthoc", "Result-known post hoc"), ("revision", "Revision, frozen in two stages")):
        ax.plot(x, 0.17, "o", ms=5, color=STATUS[key][0])
        ax.text(x + 0.09, 0.17, label, ha="left", va="center", fontsize=8, color=INK)
        x += 0.17 + 0.062 * len(label)

    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "Fig1_rgba.png", dpi=600)
    from PIL import Image
    Image.open(out / "Fig1_rgba.png").convert("RGB").save(out / "Fig1.tif", compression="tiff_lzw", dpi=(600, 600))
    (out / "Fig1_rgba.png").unlink()
    fig.savefig(out / "Fig1.pdf")
    fig.savefig(out / "Fig1_preview.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    draw(Path(sys.argv[1]))
