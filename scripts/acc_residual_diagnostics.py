"""R1-9: residual diagnostics for the ACC-relative activity residual (108 drugs).

Recomputed from the archived ADRS_comp_primary_108.csv (average-rank
percentiles), independent of the legacy diagnostic script.
Usage: python -I acc_residual_diagnostics.py <acc_repo_root> <out_json>
"""
import csv, json, sys
from pathlib import Path
import numpy as np
from scipy import stats

root, out = Path(sys.argv[1]), Path(sys.argv[2])
rows = list(csv.DictReader(open(root / "results/primary_analysis/ADRS_comp_primary_108.csv", encoding="utf-8-sig")))
drug = [r["drug"] for r in rows]
y = np.array([float(r["MIPE_potency_pct"]) for r in rows])
x = np.array([float(r["NCI60_potency_pct"]) for r in rows])
c = np.array([float(r["C_ACC_pct"]) for r in rows])
archived_resid = np.array([float(r["residual"]) for r in rows])

def pct(v):
    return (stats.rankdata(v, method="average") - 1) / (len(v) - 1)

def top20(score):
    return set(np.array(drug)[np.lexsort((np.array(drug), -score))][:20])

ols = stats.linregress(x, y)
res_ols = y - (ols.intercept + ols.slope * x)
X = np.column_stack([np.ones_like(x), x])
h = np.einsum("ij,jk,ik->i", X, np.linalg.inv(X.T @ X), X)
hc3 = np.linalg.inv(X.T @ X) @ (X.T * (res_ols / (1 - h)) ** 2) @ X @ np.linalg.inv(X.T @ X)
ts = stats.theilslopes(y, x)
res_ts = y - (ts.intercept + ts.slope * x)
rng = np.random.default_rng(20261019)
boot = []
for _ in range(5000):
    i = rng.integers(0, len(x), len(x))
    boot.append(stats.linregress(x[i], y[i]).slope)
adrs_ols = 0.5 * c + 0.5 * pct(res_ols)
adrs_ts = 0.5 * c + 0.5 * pct(res_ts)
# post hoc sensitivity: drop the six most pan-cancer-potent drugs (top ~5% of NCI-60 potency)
keep = x < np.sort(x)[-6]
ols_k = stats.linregress(x[keep], y[keep])
res_k = y[keep] - (ols_k.intercept + ols_k.slope * x[keep])
summary = {
    "n": len(x),
    "ols_slope": ols.slope, "ols_intercept": ols.intercept, "ols_r": ols.rvalue,
    "ols_slope_hc3_se": float(np.sqrt(hc3[1, 1])),
    "ols_slope_bootstrap95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
    "max_abs_diff_vs_archived_residual": float(np.max(np.abs(res_ols - archived_resid))),
    "shapiro_W": float(stats.shapiro(res_ols).statistic), "shapiro_P": float(stats.shapiro(res_ols).pvalue),
    "residual_skew": float(stats.skew(res_ols)), "residual_excess_kurtosis": float(stats.kurtosis(res_ols)),
    "theil_sen_slope": ts.slope, "theil_sen_95": [ts.low_slope, ts.high_slope],
    "spearman_resid_ols_vs_ts": float(stats.spearmanr(res_ols, res_ts).statistic),
    "spearman_adrs_ols_vs_ts": float(stats.spearmanr(adrs_ols, adrs_ts).statistic),
    "top20_overlap_adrs_ols_vs_ts": len(top20(adrs_ols) & top20(adrs_ts)),
    "posthoc_drop6_most_potent": {"dropped": sorted(np.array(drug)[~keep].tolist()), "slope": ols_k.slope,
                                   "spearman_resid_rank_vs_full": float(stats.spearmanr(res_k, res_ols[keep]).statistic)},
}
out.write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
