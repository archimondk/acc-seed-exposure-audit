"""R1-12: inventory of inferential tests and study-wide BH sensitivity.

Family H (hypothesis-level): every summary-level inferential P value reported
in the original manuscript or added in the revision. Values come from archived
result files where they exist; the few that exist only in the manuscript text
are entered verbatim and marked source='manuscript'.
Family D (drug-level): every per-drug degree-matched empirical P value
(ACC 108 x 4 variants, breast 108 x 4, melanoma 124 x 4 x 2 arms), pooled.

Usage: python -I test_inventory_global_bh.py <acc_repo_root> <melanoma_stageB_run_dir> <out_dir>
"""
import csv
import json
import sys
from pathlib import Path


def bh(p):
    m = len(p)
    order = sorted(range(m), key=lambda i: (p[i], i))
    q, running = [1.0] * m, 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, p[i] * m / rank)
        q[i] = min(1.0, running)
    return q


def read(path):
    with open(path, encoding="utf-8-sig", newline="") as h:
        return list(csv.DictReader(h))


root, run, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
out.mkdir(parents=True, exist_ok=True)
H = []


def add(test, p, status, family_in_paper, source):
    H.append({"test": test, "p": float(p), "status": status, "original_family": family_in_paper, "source": source})


for r in read(root / "results/mechanism_enrichment/mechanism_enrichment_primary108.csv"):
    add(f"ACC mechanism-class rank sum: {r['mechanism_class']}", r["p_exact"], "original", "10 mechanism classes (BH)",
        "mechanism_enrichment_primary108.csv")
for r in read(root / "results/normalization_sensitivity/normalization_sensitivity_summary.csv"):
    add(f"ACC CDK4/6 group vs degree-matched null: {r['variant']}", r["CDK46_degree_matched_empirical_p"],
        "original (primary) / post hoc (rank variants)", "4 variants (BH)", "normalization_sensitivity_summary.csv")
for r in read(root / "results/positive_control/positive_control_variant_summary.csv"):
    add(f"Breast CDK4/6 group vs null: {r['variant']}", r["CDK46_empirical_p_upper"], "frozen-before-result",
        "4 variants (BH)", "positive_control_variant_summary.csv")
for r in read(root / "leakage_audit/amendment1_scale_free_effect.csv"):
    add(f"Amendment 1 abemaciclib vs unexposed: {r['disease']} {r['variant']}", r["empirical_P_vs_unexposed"],
        "result-known post hoc", "none", "amendment1_scale_free_effect.csv")
decisions = json.loads((run / "decisions.json").read_text(encoding="utf-8"))
for v, g in decisions["H4"]["124"]["group_by_variant"].items():
    add(f"Melanoma BRAF/MEK group vs null (124): {v}", g["p"], "frozen-before-result (revision)", "4 variants (BH)",
        "melanoma decisions.json")
for test, p, status in [
    ("Broad cytotoxic rank shift, Mann-Whitney", 0.034, "original"),
    ("Broad cytotoxics in Top 20 after adjustment, hypergeometric", 0.17, "original"),
    ("Broad cytotoxics in Top 20 before adjustment, hypergeometric", 0.07, "original"),
    ("S_B-neighbor LOMCO Spearman", 0.266, "original"),
    ("MIPE missingness by mechanism group, Fisher", 0.0237, "original"),
    ("Partial Spearman C_ACC vs matched-null z given L_count", 0.667, "exploratory"),
    ("L_d(weighted) vs C_ACC Spearman", 2.37e-28, "exploratory"),
    ("Cross-platform concordance CU-ACC1", 0.001, "original (reported as P < 0.001; bound used)"),
    ("Cross-platform concordance CU-ACC2", 0.008, "original"),
    ("Cross-platform concordance NCI-H295R", 0.013, "original"),
]:
    add(test, p, status, "none", "manuscript")
qh = bh([h["p"] for h in H])
for h, q in zip(H, qh):
    h["q_study_wide_bh"] = q
    h["survives_q_lt_0_05"] = q < 0.05
H.sort(key=lambda h: h["p"])
with open(out / "test_inventory_hypothesis_level.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(H[0]), lineterminator="\n")
    w.writeheader()
    w.writerows(H)

D = []
for r in read(root / "results/normalization_sensitivity/normalization_sensitivity_null_primary108.csv"):
    D.append(("ACC", r["variant"], r["drug"], float(r["empirical_p_upper"])))
for r in read(root / "results/positive_control/positive_control_drug_ranks_primary108.csv"):
    D.append(("Breast", r["variant"], r["drug"], float(r["empirical_p_upper"])))
for r in read(run / "drug_statistics.csv"):
    if r["universe"] == "124":
        D.append((f"Melanoma {r['arm']}", r["variant"], r["drug"], float(r["p_upper"])))
qd = bh([d[3] for d in D])
rows = [{"disease": a, "variant": v, "drug": d, "p": p, "q_pooled_bh": q} for (a, v, d, p), q in zip(D, qd)]
rows.sort(key=lambda r: r["p"])
with open(out / "test_inventory_drug_level_pooled.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
summary = {"hypothesis_level_tests": len(H), "surviving_hypothesis_level": [h["test"] for h in H if h["survives_q_lt_0_05"]],
           "drug_level_tests": len(rows), "min_pooled_drug_q": rows[0]["q_pooled_bh"],
           "drug_level_q_lt_0_05": [(r["disease"], r["variant"], r["drug"]) for r in rows if r["q_pooled_bh"] < 0.05],
           "smallest_drug_level": [(r["disease"], r["variant"], r["drug"], r["p"], round(r["q_pooled_bh"], 4)) for r in rows[:8]]}
(out / "test_inventory_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
for h in H[:14]:
    print(f"{h['p']:.3g}\t{h['q_study_wide_bh']:.3g}\t{h['test']}")
