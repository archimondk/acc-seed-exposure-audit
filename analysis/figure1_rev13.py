"""Render the rev13 status-labelled Figure 1 for the Pharmaceutics manuscript.

The diagram is constructed entirely from vector primitives and editable text.
It contains no generative-image content and no data-derived geometry.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "revision" / "Fig1_workflow_status_audit_rev13"

COLORS = {
    "ink": "#202124",
    "muted": "#5F6368",
    "line": "#667078",
    "panel": "#F7F8F9",
    "blue": "#0072B2",
    "blue_fill": "#E7F2F8",
    "teal": "#009E73",
    "teal_fill": "#E5F4EF",
    "orange": "#D55E00",
    "orange_fill": "#FCEDE3",
    "purple": "#8B5FBF",
    "purple_fill": "#F1EAF8",
    "grey_fill": "#F0F2F3",
    "original": "#0072B2",
    "exploratory": "#E69F00",
    "frozen": "#009E73",
    "posthoc": "#8B5FBF",
}


mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
    }
)


def rounded_box(ax, x, y, w, h, text, edge, fill="white", *, fontsize=8.0,
                weight="normal", linestyle="-", status=None, zorder=3):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.006,rounding_size=0.012",
        linewidth=1.15,
        edgecolor=edge,
        facecolor=fill,
        linestyle=linestyle,
        zorder=zorder,
    )
    ax.add_patch(patch)
    text_y = y + h * (0.29 if status else 0.50)
    ax.text(
        x + w / 2, text_y, text,
        ha="center", va="center", fontsize=fontsize, color=COLORS["ink"],
        fontweight=weight, linespacing=1.14, zorder=zorder + 1,
    )
    if status:
        badge(ax, x + w - 0.009, y + h - 0.007, status, anchor="right",
              zorder=zorder + 2)
    return patch


def badge(ax, x, y, status, *, anchor="left", zorder=6):
    labels = {
        "original": "ORIG.",
        "exploratory": "EXPL.",
        "frozen": "FROZEN",
        "posthoc": "POST-HOC",
    }
    ax.text(
        x, y, labels[status],
        ha=anchor, va="top", fontsize=6.6, color="white", fontweight="bold",
        bbox={
            "boxstyle": "round,pad=0.20,rounding_size=0.7",
            "facecolor": COLORS[status],
            "edgecolor": COLORS[status],
            "linewidth": 0.6,
        },
        zorder=zorder,
    )


def arrow(ax, start, end, *, color=None, style="-|>", lw=1.35,
          connectionstyle="arc3", zorder=2):
    arr = FancyArrowPatch(
        start, end, arrowstyle=style, mutation_scale=11,
        linewidth=lw, color=color or COLORS["line"],
        connectionstyle=connectionstyle, zorder=zorder,
    )
    ax.add_patch(arr)
    return arr


def panel(ax, x, y, w, h, number, title, edge, fill):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.008,rounding_size=0.018",
        linewidth=1.45, edgecolor=edge, facecolor=fill, zorder=0,
    )
    ax.add_patch(patch)
    ax.text(x + 0.017, y + h - 0.027, number, ha="left", va="top",
            fontsize=11, fontweight="bold", color=edge)
    ax.text(x + 0.046, y + h - 0.030, title, ha="left", va="top",
            fontsize=9.0, fontweight="bold", color=edge)
    return patch


def draw():
    width_in = 170 / 25.4
    height_in = 164 / 25.4
    fig, ax = plt.subplots(figsize=(width_in, height_in))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(
        0.5, 0.975, "Study workflow and status-labelled audit structure",
        ha="center", va="top", fontsize=12.2, fontweight="bold",
        color=COLORS["ink"],
    )

    ax.text(0.020, 0.918, "Status", ha="left", va="center",
            fontsize=8.0, fontweight="bold", color=COLORS["muted"])
    key = [
        (0.105, 0.120, "original", "Original"),
        (0.240, 0.155, "exploratory", "Exploratory"),
        (0.410, 0.240, "frozen", "Frozen before result"),
        (0.665, 0.285, "posthoc", "Result-known post-hoc"),
    ]
    for x, w, status, label in key:
        p = FancyBboxPatch(
            (x, 0.900), w, 0.038,
            boxstyle="round,pad=0.002,rounding_size=0.012",
            linewidth=0.8, edgecolor=COLORS[status],
            facecolor=COLORS[status], zorder=4,
        )
        ax.add_patch(p)
        ax.text(x + w / 2, 0.919, label, ha="center", va="center",
                fontsize=7.1, fontweight="bold", color="white", zorder=5)

    y0, ph = 0.185, 0.675
    xs = [0.020, 0.252, 0.494, 0.736]
    ws = [0.220, 0.230, 0.230, 0.244]
    panel(ax, xs[0], y0, ws[0], ph, "1", "Inputs",
          COLORS["blue"], COLORS["blue_fill"])
    panel(ax, xs[1], y0, ws[1], ph, "2", "Context",
          COLORS["teal"], COLORS["teal_fill"])
    panel(ax, xs[2], y0, ws[2], ph, "3", "Ranking",
          COLORS["orange"], COLORS["orange_fill"])
    panel(ax, xs[3], y0, ws[3], ph, "4", "Audits",
          COLORS["muted"], COLORS["grey_fill"])

    # Inter-panel flow arrows sit behind the boxes.
    arrow(ax, (0.236, 0.505), (0.258, 0.505))
    arrow(ax, (0.478, 0.505), (0.500, 0.505))
    arrow(ax, (0.720, 0.505), (0.742, 0.505))

    # Panel 1: data sources (inputs are not assigned analysis-status labels).
    for y, text in [
        (0.695, "B-index network\n124 drugs\n399 genes\n1304 associations"),
        (0.565, "STRING v12\nprotein network"),
        (0.435, "MIPE 5.0\n3 ACC cell lines"),
    ]:
        rounded_box(ax, 0.043, y, 0.174, 0.096, text, COLORS["blue"],
                    fill="white", fontsize=7.2)
    rounded_box(ax, 0.043, 0.305, 0.174, 0.096,
                "NCI-60 +\nACC_CellMinerCDB", COLORS["blue"],
                fill="white", fontsize=6.2)

    # Panel 2: original construction.
    rounded_box(ax, 0.276, 0.705, 0.182, 0.087,
                "45 disease-only\nseeds", COLORS["teal"], fill="white",
                fontsize=7.5,
                status="original")
    rounded_box(ax, 0.276, 0.570, 0.182, 0.103,
                "RWR (α = 0.40)\n→ r_ACC", COLORS["teal"], fill="white",
                status="original")
    rounded_box(ax, 0.276, 0.435, 0.182, 0.103,
                "Gene-set mean\n+ shrinkage\n→ C_ACC", COLORS["teal"],
                fill="white", fontsize=7.2, status="original")
    rounded_box(ax, 0.276, 0.285, 0.182, 0.118,
                "MIPE − NCI-60\nbaseline\n→ ACC residual",
                COLORS["teal"], fill="white", fontsize=7.2,
                status="original")
    arrow(ax, (0.367, 0.704), (0.367, 0.674))
    arrow(ax, (0.367, 0.569), (0.367, 0.539))
    arrow(ax, (0.367, 0.434), (0.367, 0.404))

    # Panel 3: universes and locked score.
    rounded_box(ax, 0.518, 0.718, 0.182, 0.070,
                "124-drug\ncontext universe", COLORS["orange"], fill="white",
                fontsize=7.2,
                status="original")
    rounded_box(ax, 0.518, 0.612, 0.182, 0.080,
                "Complete-case gate\n108 primary\n16 context-only",
                COLORS["orange"], fill="white", fontsize=7.0,
                status="original")
    rounded_box(ax, 0.518, 0.470, 0.182, 0.112,
                "ADRS_comp\n0.50·P(C_ACC)\n+ 0.50·P(residual)",
                COLORS["orange"], fill="white", fontsize=6.6,
                weight="bold", status="original")
    rounded_box(ax, 0.518, 0.355, 0.182, 0.085,
                "Mechanism tests\n21-setting\nweight scan",
                COLORS["orange"], fill="white", fontsize=7.0,
                status="original")
    rounded_box(ax, 0.518, 0.235, 0.182, 0.084,
                "Supplement only\nLiterature-informed\nreprioritization",
                COLORS["purple"], fill=COLORS["purple_fill"],
                linestyle="--", fontsize=7.2, status="exploratory")
    arrow(ax, (0.609, 0.717), (0.609, 0.693))
    arrow(ax, (0.609, 0.611), (0.609, 0.583))
    arrow(ax, (0.609, 0.469), (0.609, 0.441))

    # Panel 4: each analysis block has an explicit text badge.
    audit_boxes = [
        (0.720, "Bias/null audits\nDirect seed-overlap\nbaseline", "exploratory", 7.1),
        (0.643, "Breast control gate\nRB1 intervention", "frozen", 7.1),
        (0.566, "A1: scale-free\nre-analysis", "posthoc", 7.1),
        (0.489, "A2–3: weight analyses\n+ W2–MIPE", "posthoc", 7.1),
        (0.412, "A4: 45-seed scan", "posthoc", 7.1),
        (0.335, "A5: seed-excluded\nscoring", "posthoc", 7.1),
        (0.258, "Evidence-label audit\ntraceability only", "exploratory", 7.1),
    ]
    for y, text, status, font_size in audit_boxes:
        rounded_box(ax, 0.758, y, 0.200, 0.070, text, COLORS[status],
                    fill="white", fontsize=font_size, status=status)

    # Bounded output, deliberately separated from efficacy claims.
    arrow(ax, (0.858, 0.185), (0.858, 0.153), color=COLORS["teal"])
    rounded_box(
        ax, 0.208, 0.045, 0.584, 0.088,
        "Audit and hypotheses for testing — not efficacy predictions",
        COLORS["teal"], fill=COLORS["teal_fill"], fontsize=9.4,
        weight="bold",
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT.with_suffix(".svg"), bbox_inches=None)
    fig.savefig(OUT.with_suffix(".pdf"), bbox_inches=None)
    fig.savefig(OUT.with_suffix(".png"), dpi=1000, bbox_inches=None)
    plt.close(fig)


if __name__ == "__main__":
    draw()
