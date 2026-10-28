"""Fig 3 (revision R1; replaces the original Fig 5): seed exposure, interventions and seed exclusion in three diseases.

A  Rank barcodes, 108-drug universe: every drug at its rank under each
   normalization; directly exposed drugs dark, unexposed light; top-quartile line.
B  Unexposed drugs in the top quartile against their attainable range.
C  Seed-intervention effects (delta z) for the directly manipulated drugs.
D  Melanoma: ranks of the six directly exposed drugs before and after seed exclusion.

Inputs: ACC repository (archived ACC and breast outputs) and the melanoma
stage-B outputs. Writes Fig3.tif (RGB, LZW, 600 dpi), Fig3.pdf, a PNG preview and
Fig3_source_data.csv.

Usage: python -I fig3_three_disease_r1.py <acc_repo_root> <revision_root> <out_dir>
"""
import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 8, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.labelsize": 8,
                     "ytick.labelsize": 8, "axes.labelsize": 8, "pdf.fonttype": 42, "savefig.facecolor": "white"})
INK, MUTED, GRID = "#202124", "#5F6368", "#E3E6E8"
EXPOSED, UNEXPOSED = "#202124", "#C4CAD0"
DISEASE = {"ACC": "#0072B2", "Breast": "#E69F00", "Melanoma": "#009E73"}  # validated (light)
VARIANTS = ["column_minmax", "column_gene_rank", "uniform_ratio_gene_rank", "symmetric_gene_rank"]
VLABEL = {"column_minmax": "Column min–max (primary)", "column_gene_rank": "Column rank",
          "uniform_ratio_gene_rank": "Uniform-ratio rank", "symmetric_gene_rank": "Symmetric rank"}
VSHORT = {"column_minmax": "Primary", "column_gene_rank": "Col.\nrank",
          "uniform_ratio_gene_rank": "Unif.\nrank", "symmetric_gene_rank": "Symm.\nrank"}
TOP = 27


def read(p):
    with open(p, encoding="utf-8-sig", newline="") as h:
        return list(csv.DictReader(h))


def load(root: Path, rev: Path):
    assoc = {}
    for e in read(root / "data/bindex_network/bindex_edges_1304.csv"):
        assoc.setdefault(e["drug"], set()).add(e["gene"])
    lower = {d.lower(): d for d in assoc}
    seeds = {
        "ACC": set(json.loads((root / "leakage_audit/observational/seed_set_used.json").read_text())["weights"]),
        "Breast": {r["gene"] for r in read(root / "results/positive_control/positive_control_seed_frozen.csv")
                   if r["include_primary"] == "yes"},
        "Melanoma": {r["gene"] for r in read(rev / "results/melanoma_stageA_seeds/primary_metastasis/seeds_frozen.csv")
                     if r["arm"] == "M0"},
    }
    pct = {d: {} for d in seeds}
    for r in read(root / "results/normalization_sensitivity/normalization_sensitivity_null_primary108.csv"):
        pct["ACC"].setdefault(r["variant"], {})[lower[r["drug"].lower()]] = float(r["observed_C_ACC_pct"])
    for r in read(root / "results/positive_control/positive_control_drug_ranks_primary108.csv"):
        pct["Breast"].setdefault(r["variant"], {})[r["drug"]] = float(r["observed_C_ACC_percentile"])
    for r in read(rev / "results/melanoma_stageB/run/drug_statistics.csv"):
        if r["arm"] == "M0" and r["universe"] == "108":
            pct["Melanoma"].setdefault(r["variant"], {})[r["drug"]] = float(r["percentile"])
    return assoc, seeds, pct


def ranking(p):
    order = sorted(p, key=lambda d: (-p[d], d))
    return {d: i + 1 for i, d in enumerate(order)}


def main(root: Path, rev: Path, out: Path):
    out.mkdir(parents=True, exist_ok=True)
    assoc, seeds, pct = load(root, rev)
    source = []
    fig = plt.figure(figsize=(7.5, 7.0))
    gtop = fig.add_gridspec(1, 1, left=0.205, right=0.985, top=0.965, bottom=0.555)
    gbot = fig.add_gridspec(1, 3, left=0.105, right=0.965, top=0.445, bottom=0.175, wspace=0.62)

    # ---------------------------------------------------------------- A
    axA = fig.add_subplot(gtop[0, 0])
    y, yticks, ylabels, stats = 0, [], [], {}
    names = {"ACC": "ACC", "Breast": "ER+/HER2− breast cancer", "Melanoma": "BRAF V600E/K melanoma"}
    for disease in ("ACC", "Breast", "Melanoma"):
        stats[disease] = {}
        n_exp = sum(1 for d in pct[disease]["column_minmax"] if assoc[d] & seeds[disease])
        axA.text(0.6, y + 0.62, f"{names[disease]}: {n_exp} of 108 drugs directly exposed", ha="left",
                 va="center", fontsize=8, fontweight="bold", color=DISEASE[disease], zorder=5,
                 bbox=dict(boxstyle="square,pad=0.1", facecolor="white", edgecolor="none"))
        for v in VARIANTS:
            rk = ranking(pct[disease][v])
            exp = {d for d in rk if assoc[d] & seeds[disease]}
            xs_e = [rk[d] for d in exp]
            xs_u = [rk[d] for d in rk if d not in exp]
            axA.vlines(xs_u, y - 0.22, y + 0.22, color=UNEXPOSED, lw=0.9)
            axA.vlines(xs_e, y - 0.36, y + 0.36, color=EXPOSED, lw=1.1)
            n_un_top = sum(1 for d in rk if d not in exp and rk[d] <= TOP)
            stats[disease][v] = (n_un_top, max(TOP - len(exp), 0), min(TOP, len(rk) - len(exp)), len(exp))
            for d in rk:
                source.append({"panel": "A", "disease": disease, "variant": v, "drug": d, "rank_108": rk[d],
                               "directly_exposed": d in exp})
            yticks.append(y)
            ylabels.append(VLABEL[v])
            y -= 1
        y -= 0.75
    axA.axvline(TOP + 0.5, color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
    axA.text(TOP + 1.2, y + 0.82, "Top quartile (rank ≤ 27)", ha="left", va="center", fontsize=8, color=MUTED)
    axA.set_yticks(yticks, ylabels)
    axA.set_xlim(0.3, 108.7)
    axA.set_ylim(y + 0.55, 1.0)
    axA.set_xticks([1, 27, 54, 81, 108])
    axA.set_xlabel("Rank in the 108-drug universe (1 = highest context score)")
    for s in ("top", "right", "left"):
        axA.spines[s].set_visible(False)
    axA.tick_params(axis="y", length=0)
    axA.vlines([], 0, 0)
    axA.plot([], [], color=EXPOSED, lw=1.2, label="Directly exposed drug")
    axA.plot([], [], color=UNEXPOSED, lw=1.2, label="Unexposed drug")
    axA.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, frameon=False, fontsize=8,
               handlelength=1.2, borderaxespad=0.1)
    axA.text(-0.255, 1.0, "A", transform=axA.transAxes, fontsize=11, fontweight="bold", va="bottom")

    # ---------------------------------------------------------------- B
    axB = fig.add_subplot(gbot[0, 0])
    rows = [("ACC", 2), ("Breast", 1), ("Melanoma", 0)]
    for disease, yy in rows:
        lo, hi = stats[disease]["column_minmax"][1], stats[disease]["column_minmax"][2]
        axB.plot([lo, hi], [yy, yy], color="#DADDE1", lw=7, solid_capstyle="round", zorder=1)
        rb = [stats[disease][v][0] for v in VARIANTS[1:]]
        axB.scatter(rb, [yy] * 3, s=30, facecolor="white", edgecolor=DISEASE[disease], lw=1.2, zorder=3,
                    marker="o")
        prim = stats[disease]["column_minmax"][0]
        axB.scatter([prim], [yy], s=42, color=DISEASE[disease], edgecolor="white", lw=1.0, zorder=4, marker="o")
        axB.text(prim, yy + 0.30, f"{prim}", ha="center", va="bottom", fontsize=8, color=INK)
        source.append({"panel": "B", "disease": disease, "variant": "column_minmax",
                       "unexposed_in_top27": prim, "attainable_min": lo, "attainable_max": hi})
        for v in VARIANTS[1:]:
            source.append({"panel": "B", "disease": disease, "variant": v, "unexposed_in_top27": stats[disease][v][0],
                           "attainable_min": lo, "attainable_max": hi})
    axB.set_yticks([2, 1, 0], ["ACC", "Breast", "Melanoma"])
    axB.set_ylim(-0.6, 2.75)
    axB.set_xlim(-1.5, 28.5)
    axB.set_xticks([0, 9, 18, 27])
    axB.set_xlabel("Unexposed drugs in top 27")
    axB.scatter([], [], s=42, color=INK, label="Min–max (primary)")
    axB.scatter([], [], s=30, facecolor="white", edgecolor=INK, lw=1.2, label="Rank-based (×3)")
    axB.plot([], [], color="#DADDE1", lw=6, label="Attainable range")
    axB.legend(loc="upper left", bbox_to_anchor=(-0.35, -0.24), ncol=1, frameon=False, fontsize=8,
               handletextpad=0.5, borderaxespad=0)
    for s in ("top", "right"):
        axB.spines[s].set_visible(False)
    axB.text(-0.42, 1.02, "B", transform=axB.transAxes, fontsize=11, fontweight="bold", va="bottom")

    # ---------------------------------------------------------------- C
    axC = fig.add_subplot(gbot[0, 1])
    curve = {(r["variant"], r["drug"]): r for r in read(root / "leakage_audit/leakage_curve_108.csv")}
    h6 = {(r["variant"], r["drug"]): float(r["delta_z"]) for r in read(rev / "results/melanoma_stageB/diagnostics/H6_delta_z_124.csv")}
    series = [
        ("Abemaciclib, ACC −RB1", "ACC", "o", [float(curve[(v, "Abemaciclib")]["delta_z_acc_A1_minus_A2"]) for v in VARIANTS]),
        ("Abemaciclib, breast +RB1", "Breast", "s", [float(curve[(v, "Abemaciclib")]["delta_z_breast_B2_minus_B1"]) for v in VARIANTS]),
        ("Dacarbazine, melanoma +BRAF", "Melanoma", "^", [h6[(v, "Dacarbazine")] for v in VARIANTS]),
        ("Tioguanine, melanoma +BRAF", "Melanoma", "D", [h6[(v, "Tioguanine")] for v in VARIANTS]),
    ]
    axC.axhspan(-0.5, 0.5, color="#EEF0F2", zorder=0)
    axC.axhline(2.0, color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
    axC.text(3.45, 2.08, "Δz = 2", ha="right", va="bottom", fontsize=8, color=MUTED)
    axC.text(-0.45, -0.45, "±0.5 band", ha="left", va="bottom", fontsize=8, color=MUTED)
    offs = [-0.21, -0.07, 0.07, 0.21]
    for (label, disease, marker, vals), off in zip(series, offs):
        xs = [i + off for i in range(4)]
        filled = marker != "D"
        axC.scatter(xs, vals, marker=marker, s=26, zorder=3, color=DISEASE[disease] if filled else "white",
                    edgecolor=DISEASE[disease], lw=1.1, label=label)
        for v, val in zip(VARIANTS, vals):
            source.append({"panel": "C", "series": label, "variant": v, "delta_z": val})
    axC.set_xticks(range(4), [VSHORT[v] for v in VARIANTS])
    axC.set_xlim(-0.5, 3.5)
    axC.set_ylim(-1.0, 4.4)
    axC.set_ylabel("Intervention effect, Δz")
    axC.legend(loc="upper left", bbox_to_anchor=(-0.30, -0.24), frameon=False, fontsize=8, handletextpad=0.3,
               borderaxespad=0)
    for s in ("top", "right"):
        axC.spines[s].set_visible(False)
    axC.text(-0.36, 1.02, "C", transform=axC.transAxes, fontsize=11, fontweight="bold", va="bottom")

    # ---------------------------------------------------------------- D
    axD = fig.add_subplot(gbot[0, 2])
    dec = [r for r in read(rev / "results/melanoma_stageB/diagnostics/seed_nonseed_decomposition.csv")
           if r["arm"] == "M0" and r["universe"] == "124" and r["variant"] == "column_minmax" and r["seed_hits"]]
    dec.sort(key=lambda r: int(r["rank_full"]))
    ends = sorted((int(r["rank_nonseed"]), r["drug"]) for r in dec)
    label_y, last = {}, -1e9
    for b, name in ends:  # keep labels >= 8 ranks apart, then centre the stack on the data
        pos = max(b, last + 8)
        label_y[name], last = pos, pos
    shift = (sum(b for b, _ in ends) - sum(label_y.values())) / len(ends)
    label_y = {k: v + shift for k, v in label_y.items()}
    for r in dec:
        a, b = int(r["rank_full"]), int(r["rank_nonseed"])
        axD.plot([0, 1], [a, b], color=DISEASE["Melanoma"], lw=1.0, marker="o", ms=3.5, zorder=3)
        ly = label_y[r["drug"]]
        axD.plot([1.03, 1.17], [b, ly], color=MUTED, lw=0.5, zorder=2)
        axD.text(1.2, ly, r["drug"], ha="left", va="center", fontsize=8, color=INK)
        source.append({"panel": "D", "drug": r["drug"], "rank_locked_124": a, "rank_seed_excluded_124": b,
                       "seed_hits": r["seed_hits"]})
    axD.axhline(31.5, color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
    axD.text(1.05, 29.5, "Top quartile", ha="left", va="bottom", fontsize=8, color=MUTED)
    axD.set_xticks([0, 1], ["Locked", "Seed-\nexcluded"])
    axD.set_xlim(-0.2, 2.45)
    axD.set_ylim(124.5, -3)
    axD.set_yticks([1, 31, 62, 93, 124])
    axD.set_ylabel("Rank among 124 drugs")
    for s in ("top", "right"):
        axD.spines[s].set_visible(False)
    axD.text(-0.42, 1.02, "D", transform=axD.transAxes, fontsize=11, fontweight="bold", va="bottom")

    fig.savefig(out / "Fig3_rgba.png", dpi=600)
    Image.open(out / "Fig3_rgba.png").convert("RGB").save(out / "Fig3.tif", compression="tiff_lzw", dpi=(600, 600))
    (out / "Fig3_rgba.png").unlink()
    fig.savefig(out / "Fig3.pdf")
    fig.savefig(out / "Fig3_preview.png", dpi=150)
    plt.close(fig)
    keys = sorted({k for r in source for k in r}, key=lambda k: ["panel", "disease", "series", "variant", "drug"].index(k)
                  if k in ["panel", "disease", "series", "variant", "drug"] else 99)
    with open(out / "Fig3_source_data.csv", "w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=keys, lineterminator="\n")
        w.writeheader()
        w.writerows(source)
    print(json.dumps({d: {v: stats[d][v] for v in VARIANTS} for d in stats}, indent=1))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
