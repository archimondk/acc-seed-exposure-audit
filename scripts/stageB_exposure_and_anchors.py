"""Stage B, steps 1-2 for protocol melanoma_R1_v0.3 (no propagation).

Step 1: exposure map of the stage-A melanoma seeds (primary M0/M1) on the
frozen 399-gene drug association network: per-seed drug lists, per-drug
exposure counts, L_count / L_weight, coverage c in the 108 and 124 universes,
status of the six BRAF/MEK drugs and the intervention set X.

Step 2: recheck of the ACC and breast anchors for H5 from archived outputs
(108 universe, primary column_minmax variant), for both L_count and L_weight.

Reads only: bindex association files, the archived ACC/breast seed tables and
archived ACC/breast C percentiles. It does not read STRING or compute any
melanoma drug score.

Usage:
    python -I stageB_exposure_and_anchors.py <acc_repo_root> <stageA_results_dir> <out_dir>
"""
import csv
import json
import sys
from fractions import Fraction
from pathlib import Path

from scipy import stats

K = 3
SIX = ["Dabrafenib", "Vemurafenib", "Encorafenib", "Trametinib", "Cobimetinib", "Binimetinib"]
TOP25 = {108: 27, 124: 31}


def read(path):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.12g}" if isinstance(v, float) else v) for k, v in r.items()})


def load_network(root):
    net = root / "data" / "bindex_network"
    assoc = {}
    for e in read(net / "bindex_edges_1304.csv"):
        assoc.setdefault(e["drug"], set()).add(e["gene"])
    genes = {r["gene"] for r in read(net / "bindex_genes_399.csv")}
    drugs = [r["drug"] for r in read(net / "bindex_drugs_124.csv")]
    u108 = {r["drug"] for r in read(root / "results" / "positive_control" / "positive_control_drug_ranks_primary108.csv")}
    assert len(drugs) == 124 and len(u108) == 108 and len(genes) == 399
    return assoc, genes, drugs, u108


def overlap(assoc, drug, weights):
    hits = sorted(assoc[drug] & set(weights))
    n = len(assoc[drug])
    return hits, Fraction(len(hits), n + K), sum(weights[g] for g in hits) / (n + K)


def top_set(pct, universe_n):
    order = sorted(pct, key=lambda d: (-pct[d], d))
    return set(order[:TOP25[universe_n]])


def anchor(name, assoc, weights, pct):
    drugs = sorted(pct)
    n = len(drugs)
    lc, lw, exposed = {}, {}, set()
    for d in drugs:
        hits, c, w = overlap(assoc, d, weights)
        lc[d], lw[d] = c, w
        if hits:
            exposed.add(d)
    unexposed = [d for d in drugs if d not in exposed]
    top = top_set(pct, n)
    q75 = {d for d in drugs if pct[d] >= 0.75}
    rho_c = stats.spearmanr([float(lc[d]) for d in drugs], [pct[d] for d in drugs]).statistic
    rho_w = stats.spearmanr([lw[d] for d in drugs], [pct[d] for d in drugs]).statistic
    u_top = Fraction(len(set(unexposed) & top), len(unexposed))
    return {
        "disease": name, "universe": n,
        "exposed": len(exposed), "unexposed": len(unexposed),
        "c": f"{len(exposed)}/{n}", "c_float": len(exposed) / n,
        "rho_L_count": float(rho_c), "rho_L_weight": float(rho_w),
        "unexposed_in_top_rank_cutoff": len(set(unexposed) & top),
        "unexposed_in_pct_ge_0.75": len(set(unexposed) & q75),
        "u": f"{u_top.numerator}/{u_top.denominator}", "u_float": float(u_top),
        "max_unexposed_pct": max(pct[d] for d in unexposed) if unexposed else float("nan"),
    }


def main(root, stage_a, out):
    assoc, genes, drugs, u108 = load_network(root)
    lower = {d.lower(): d for d in drugs}

    # ---------------- step 1: melanoma exposure map (primary M0/M1)
    seed_rows = read(stage_a / "primary_metastasis" / "seeds_frozen.csv")
    arms = {}
    for r in seed_rows:
        arms.setdefault(r["arm"], {})[r["gene"]] = float(r["normalized_weight"])
    seed_map = []
    for arm, w in arms.items():
        for g in sorted(w):
            linked = sorted(d for d in drugs if g in assoc[d])
            seed_map.append({"arm": arm, "seed": g, "normalized_weight": w[g], "in_399": g in genes,
                             "n_drugs_124": len(linked), "n_drugs_108": sum(d in u108 for d in linked),
                             "drugs": ";".join(linked)})
    restored = sorted(set(arms["M1"]) - set(arms["M0"]))
    x_set = sorted(d for d in drugs if assoc[d] & set(restored))
    drug_rows = []
    for d in drugs:
        row = {"drug": d, "n_assoc": len(assoc[d]), "in_108": d in u108, "six_drug": d in SIX, "in_X": d in x_set}
        for arm in ("M0", "M1"):
            hits, c, w = overlap(assoc, d, arms[arm])
            row[f"{arm}_seed_hits"] = ";".join(hits)
            row[f"{arm}_n_hits"] = len(hits)
            row[f"{arm}_L_count"] = float(c)
            row[f"{arm}_L_weight"] = w
        drug_rows.append(row)
    coverage = {}
    for arm in ("M0", "M1"):
        for n, universe in ((108, [d for d in drugs if d in u108]), (124, drugs)):
            exposed = [d for d in universe if any(r["drug"] == d and r[f"{arm}_n_hits"] for r in drug_rows)]
            coverage[f"{arm}_{n}"] = {"exposed": len(exposed), "n": n, "c": f"{Fraction(len(exposed), n)}",
                                      "c_float": len(exposed) / n}
    six_status = {d: {"M0_hits": next(r["M0_seed_hits"] for r in drug_rows if r["drug"] == d),
                      "M1_hits": next(r["M1_seed_hits"] for r in drug_rows if r["drug"] == d)} for d in SIX}
    write(out / "exposure_map_M0_M1_seeds.csv", seed_map)
    write(out / "exposure_map_M0_M1_drugs.csv", drug_rows)

    # ---------------- step 2: anchors (108, column_minmax)
    acc_w = json.loads((root / "leakage_audit" / "observational" / "seed_set_used.json").read_text())["weights"]
    total = sum(acc_w.values())
    acc_w = {g: v / total for g, v in acc_w.items()}
    acc_pct = {lower[r["drug"].lower()]: float(r["C_ACC_pct"])
               for r in read(root / "results" / "primary_analysis" / "ADRS_comp_primary_108.csv")}
    br_raw = {r["gene"]: float(r["raw_weight"]) for r in
              read(root / "results" / "positive_control" / "positive_control_seed_frozen.csv") if r["include_primary"] == "yes"}
    total = sum(br_raw.values())
    br_w = {g: v / total for g, v in br_raw.items()}
    br_pct = {r["drug"]: float(r["observed_C_ACC_percentile"]) for r in
              read(root / "results" / "positive_control" / "positive_control_drug_ranks_primary108.csv")
              if r["variant"] == "column_minmax"}
    anchors = [anchor("ACC", assoc, acc_w, acc_pct), anchor("Breast", assoc, br_w, br_pct)]
    archived = {r["baseline"]: float(r["spearman"]) for r in
                read(root / "leakage_audit" / "observational" / "O1_replacement_baseline.csv") if r["target"] == "obs_C_ACC_pct"}
    write(out / "anchor_recheck_108_column_minmax.csv", anchors)

    summary = {"melanoma_coverage": coverage, "restored_genes_M1": restored, "X": x_set,
               "six_drug_exposure": six_status, "anchors": anchors,
               "archived_O1_acc_rho": archived,
               "seed_file": str(stage_a / "primary_metastasis" / "seeds_frozen.csv")}
    (out / "stageB_exposure_anchor_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
