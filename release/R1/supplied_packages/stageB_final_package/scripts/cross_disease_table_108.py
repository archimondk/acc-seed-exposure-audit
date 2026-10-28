"""Cross-disease seed-exposure table (108-drug universe, four normalizations).

ACC and breast values are recomputed from archived per-drug percentiles in
the ACC repository; melanoma values from the stage-B audit outputs. No
propagation is run. This is a result-known descriptive table (task 5).

Usage:
    python -I cross_disease_table_108.py <acc_repo_root> <stageB_run_dir> <out_csv>
"""
import csv
import json
import sys
from fractions import Fraction
from pathlib import Path

from scipy import stats

K, TOP = 3, 27
VARIANTS = ("column_minmax", "column_gene_rank", "uniform_ratio_gene_rank", "symmetric_gene_rank")


def read(p):
    with open(p, encoding="utf-8-sig", newline="") as h:
        return list(csv.DictReader(h))


def metrics(disease, variant, pct, assoc, seeds):
    names = sorted(pct)
    order = sorted(names, key=lambda d: (-pct[d], d))
    rank = {d: i + 1 for i, d in enumerate(order)}
    exposed = {d for d in names if assoc[d] & seeds}
    unexposed = [d for d in names if d not in exposed]
    top = set(order[:TOP])
    lc = [len(assoc[d] & seeds) / (len(assoc[d]) + K) for d in names]
    p = [pct[d] for d in names]
    n = [len(assoc[d]) for d in names]
    e = len(exposed)
    return {"disease": disease, "variant": variant, "exposed": e, "unexposed": len(unexposed),
            "c": f"{e}/{len(names)}",
            "unexposed_in_top27": len(top & set(unexposed)),
            "arith_min": max(TOP - e, 0), "arith_max": min(TOP, len(unexposed)),
            "at_arith_min": len(top & set(unexposed)) == max(TOP - e, 0),
            "first_unexposed_rank": min(rank[d] for d in unexposed),
            "exposed_in_top27": len(top & exposed),
            "rho_Lcount_C": float(stats.spearmanr(lc, p).statistic),
            "rho_nassoc_C": float(stats.spearmanr(n, p).statistic)}


def main(root, run_dir, out):
    assoc = {}
    for e in read(root / "data/bindex_network/bindex_edges_1304.csv"):
        assoc.setdefault(e["drug"], set()).add(e["gene"])
    lower = {d.lower(): d for d in assoc}
    acc_seeds = set(json.loads((root / "leakage_audit/observational/seed_set_used.json").read_text())["weights"])
    br_seeds = {r["gene"] for r in read(root / "results/positive_control/positive_control_seed_frozen.csv")
                if r["include_primary"] == "yes"}
    acc = {}
    for r in read(root / "results/normalization_sensitivity/normalization_sensitivity_null_primary108.csv"):
        acc.setdefault(r["variant"], {})[lower[r["drug"].lower()]] = float(r["observed_C_ACC_pct"])
    br = {}
    for r in read(root / "results/positive_control/positive_control_drug_ranks_primary108.csv"):
        br.setdefault(r["variant"], {})[r["drug"]] = float(r["observed_C_ACC_percentile"])
    mel_seeds = {"ARID2", "CDKN2A", "PTEN", "TP53"}
    mel = {}
    for r in read(run_dir / "drug_statistics.csv"):
        if r["arm"] == "M0" and r["universe"] == "108":
            mel.setdefault(r["variant"], {})[r["drug"]] = float(r["percentile"])
    rows = []
    for v in VARIANTS:
        rows.append(metrics("ACC", v, acc[v], assoc, acc_seeds))
        rows.append(metrics("Breast", v, br[v], assoc, br_seeds))
        rows.append(metrics("Melanoma", v, mel[v], assoc, mel_seeds))
    with open(out, "w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{x:.6f}" if isinstance(x, float) else x) for k, x in r.items()})
    for r in rows:
        print(r)


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
