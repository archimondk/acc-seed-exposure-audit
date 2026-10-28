"""Frozen decision rules of protocol melanoma_R1_v0.3 (pure functions).

H4 (PROJECT_PLAN 5.2), H5 (4.2), H6 (4.2) and supporting statistics.
No I/O and no dependence on the propagation code, so every branch can be
tested with constructed inputs.
"""
from __future__ import annotations

import math
from fractions import Fraction
from typing import Mapping, Sequence

import numpy as np
from scipy import stats

PRIMARY = "column_minmax"
VARIANTS = ("column_minmax", "column_gene_rank", "uniform_ratio_gene_rank", "symmetric_gene_rank")


# ------------------------------------------------------------ statistics
def benjamini_hochberg(p_values: Sequence[float]) -> list[float]:
    m = len(p_values)
    order = sorted(range(m), key=lambda i: (p_values[i], i))
    q = [1.0] * m
    running = 1.0
    for rank in range(m, 0, -1):
        index = order[rank - 1]
        running = min(running, p_values[index] * m / rank)
        q[index] = min(1.0, running)
    return q


def empirical_upper(observed: float, null: np.ndarray) -> tuple[int, float]:
    """Return (exceedance count b, plus-one P = (b + 1) / (B + 1)), float64."""
    null = np.asarray(null, dtype=np.float64)
    if null.ndim != 1 or null.size == 0 or not np.isfinite(null).all():
        raise ValueError("null must be a non-empty finite vector")
    b = int(np.count_nonzero(null >= np.float64(observed)))
    return b, (b + 1.0) / (null.size + 1.0)


def clopper_pearson(b: int, n: int, level: float = 0.95) -> tuple[float, float]:
    alpha = 1.0 - level
    lower = 0.0 if b == 0 else float(stats.beta.ppf(alpha / 2, b, n - b + 1))
    upper = 1.0 if b == n else float(stats.beta.ppf(1 - alpha / 2, b + 1, n - b))
    return lower, upper


def z_score(observed: float, null: np.ndarray) -> float:
    """z with ddof = 1; NA (nan) when the null has zero variance.

    Zero variance is tested exactly (all null values identical) rather than
    via std() > 0, because the floating-point std of a constant vector can be
    ~1e-17 instead of 0. No epsilon is used.
    """
    null = np.asarray(null, dtype=np.float64)
    if null.size < 2 or np.ptp(null) == 0:
        return math.nan
    sd = float(null.std(ddof=1))
    return (float(observed) - float(null.mean())) / sd


def top_k(names: Sequence[str], scores: Sequence[float], k: int) -> list[str]:
    """Descending score, ties broken by drug name (lexicographic)."""
    return [n for n, _ in sorted(zip(names, scores, strict=True), key=lambda t: (-t[1], t[0]))[:k]]


def ordinal_ranks(names: Sequence[str], scores: Sequence[float]) -> dict[str, int]:
    order = top_k(names, scores, len(names))
    return {name: i + 1 for i, name in enumerate(order)}


def spearman(x: Sequence[float], y: Sequence[float]) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return math.nan
    return float(stats.spearmanr(x, y).statistic)


# -------------------------------------------------------------------- H4
def evaluate_h4(group: Mapping[str, Mapping[str, float]], top_flags: Mapping[str, bool],
                braf_i: Sequence[str], mek_i: Sequence[str], min_top: int = 3,
                p_threshold: float = 0.05, q_threshold: float = 0.05, min_concordant: int = 3) -> dict:
    """group[variant] = {'p': group P, 'observed': mean pct, 'null_mean': mean of null group means}.

    top_flags: drug -> in Top25% under the primary variant of the evaluated universe.
    """
    p_list = [group[v]["p"] for v in VARIANTS]
    q = dict(zip(VARIANTS, benjamini_hochberg(p_list), strict=True))
    n_top = sum(bool(top_flags[d]) for d in list(braf_i) + list(mek_i))
    has_braf = any(top_flags[d] for d in braf_i)
    has_mek = any(top_flags[d] for d in mek_i)
    concordant = sum(group[v]["observed"] > group[v]["null_mean"] for v in VARIANTS)
    criteria = {
        "primary_group_p_lt_0_05": group[PRIMARY]["p"] < p_threshold,
        "primary_group_q_lt_0_05": q[PRIMARY] < q_threshold,
        "at_least_3_top25_with_one_BRAFi_and_one_MEKi": n_top >= min_top and has_braf and has_mek,
        "direction_concordant_at_least_3_variants": concordant >= min_concordant,
    }
    passed = sum(criteria.values())
    status = "pass" if passed == 4 else ("partial_recovery" if passed >= 2 else "fail")
    return {"status": status, "criteria": criteria, "q_by_variant": q, "n_six_in_top25": n_top,
            "BRAFi_in_top25": has_braf, "MEKi_in_top25": has_mek, "concordant_variants": concordant}


# -------------------------------------------------------------------- H5
def h5_feasible_bounds(l_values: Sequence[float], top_cutoff: int) -> dict:
    """Range of rho(L, C) and u attainable by any untied ordering of C given L."""
    lv = sorted((float(x) for x in l_values), reverse=True)
    n = len(lv)
    exposed = sum(x > 0 for x in lv)
    unexposed = n - exposed
    ranks_desc = list(range(n, 0, -1))
    rho_max = spearman(lv, ranks_desc)
    rho_min = spearman(lv, ranks_desc[::-1])
    if unexposed == 0:
        u_min = u_max = None
    else:
        u_min = Fraction(max(top_cutoff - exposed, 0), unexposed)
        u_max = Fraction(min(top_cutoff, unexposed), unexposed)
    return {"rho_min": rho_min, "rho_max": rho_max, "u_min": u_min, "u_max": u_max}


def _in_interval(value, low, high) -> bool:
    return (low is None or value >= low) and (high is None or value <= high)


def evaluate_h5(c_mel: Fraction, rho_mel: float, u_mel: Fraction | None,
                anchors: Mapping[str, Mapping[str, object]], bounds: Mapping[str, object] | None = None,
                tie_band: Fraction = Fraction(1, 50)) -> dict:
    """anchors: {'ACC': {'c': Fraction, 'rho': float, 'u': Fraction}, 'Breast': {...}}."""
    acc, br = anchors["ACC"], anchors["Breast"]
    if c_mel in (0, 1):
        return {"verdict": "NOT_EVALUABLE", "reason": "c equals 0 or 1"}
    near_acc = abs(c_mel - acc["c"]) < tie_band
    near_br = abs(c_mel - br["c"]) < tie_band
    if near_acc and near_br:
        return {"verdict": "NOT_EVALUABLE", "reason": "c within 0.02 of both anchors"}
    if near_acc:      # compare with breast only, one-sided
        branch, rho_iv, u_iv = "near_ACC_one_sided", (br["rho"], None), (None, br["u"])
    elif near_br:     # compare with ACC only, one-sided
        branch, rho_iv, u_iv = "near_Breast_one_sided", (None, acc["rho"]), (acc["u"], None)
    elif br["c"] < c_mel < acc["c"]:
        branch, rho_iv, u_iv = "between_anchors", (br["rho"], acc["rho"]), (acc["u"], br["u"])
    elif c_mel > acc["c"]:
        branch, rho_iv, u_iv = "above_ACC", (acc["rho"], None), (None, acc["u"])
    else:
        branch, rho_iv, u_iv = "below_Breast", (None, br["rho"]), (br["u"], None)
    rho_status = "NA" if (rho_mel is None or math.isnan(rho_mel)) else (
        "consistent" if _in_interval(rho_mel, *rho_iv) else "inconsistent")
    u_status = "NA" if u_mel is None else ("consistent" if _in_interval(u_mel, *u_iv) else "inconsistent")
    statuses = [rho_status, u_status]
    if statuses.count("NA") == 2:
        verdict = "NOT_EVALUABLE"
    elif "NA" in statuses:
        verdict = "INCOMPLETE"
    elif statuses.count("consistent") == 2:
        verdict = "SUPPORTED"
    elif statuses.count("consistent") == 1:
        verdict = "PARTIAL"
    else:
        verdict = "NOT_SUPPORTED"
    annotation = ""
    if bounds is not None:
        rho_forced = _in_interval(bounds["rho_min"], *rho_iv) and _in_interval(bounds["rho_max"], *rho_iv)
        u_forced = (bounds["u_min"] is not None and _in_interval(bounds["u_min"], *u_iv)
                    and _in_interval(bounds["u_max"], *u_iv))
        if rho_forced and u_forced:
            annotation = "ARITHMETICALLY_DETERMINED"
        elif rho_forced or u_forced:
            annotation = "PARTLY_ARITHMETICALLY_DETERMINED:" + ("rho" if rho_forced else "u")
    return {"verdict": verdict, "branch": branch, "rho_interval": rho_iv, "u_interval": u_iv,
            "rho_status": rho_status, "u_status": u_status, "annotation": annotation}


# -------------------------------------------------------------------- H6
def evaluate_h6(delta_z: Mapping[str, Mapping[str, float]], exposed_set: Sequence[str],
                six: Sequence[str], i1_min: float = 2.0, nc_max: float = 0.5,
                nc_na_fraction: float = 0.05, min_variants: int = 3) -> dict:
    """delta_z[variant][drug] = z_M1 - z_M0 (nan allowed)."""
    if not exposed_set:
        return {"verdict": "NOT_EVALUABLE", "reason": "M1 identical to M0"}
    x = list(exposed_set)

    def isna(v):
        return v is None or (isinstance(v, float) and math.isnan(v))

    def i1(variant):
        vals = [delta_z[variant][d] for d in x]
        if any(isna(v) for v in vals):
            return "NA"
        return bool(max(vals) >= i1_min and all(v > 0 for v in vals))

    i1_by_variant = {v: i1(v) for v in VARIANTS}
    i2_count = sum(i1_by_variant[v] is True for v in VARIANTS)
    i2 = i2_count >= min_variants

    def nc_a(variant):
        vals = [delta_z[variant][d] for d in six]
        if any(isna(v) for v in vals):
            return "NA"
        return bool(max(abs(v) for v in vals) <= nc_max)

    def nc_b(variant):
        others = [d for d in delta_z[variant] if d not in x]
        na = [d for d in others if isna(delta_z[variant][d])]
        vals = [abs(delta_z[variant][d]) for d in others if not isna(delta_z[variant][d])]
        max_abs = max(vals) if vals else math.nan
        if not vals or len(na) > nc_na_fraction * len(others):
            status = "NA"
        else:
            status = bool(max_abs <= nc_max)
        return {"status": status, "max_abs": max_abs, "n_na": len(na), "na_drugs": sorted(na)}

    nca = {v: nc_a(v) for v in VARIANTS}
    ncb = {v: nc_b(v) for v in VARIANTS}
    primary_x = [delta_z[PRIMARY][d] for d in x if not isna(delta_z[PRIMARY][d])]
    f_m1 = "NOT_EVALUABLE" if not primary_x else bool(all(v <= 0 for v in primary_x))
    i1p = i1_by_variant[PRIMARY]
    if i1p == "NA":
        verdict = "NOT_EVALUABLE"
    elif i1p is False:
        verdict = "NOT_SUPPORTED"
    elif nca[PRIMARY] is True and ncb[PRIMARY]["status"] is True and i2:
        verdict = "SUPPORTED"
    else:
        verdict = "PARTIAL"
    return {"verdict": verdict, "I1_by_variant": i1_by_variant, "I2_count": i2_count, "I2": i2,
            "I1_NA_variants": sum(i1_by_variant[v] == "NA" for v in VARIANTS),
            "NC_M1a_by_variant": nca, "NC_ROBUST": all(nca[v] is True for v in VARIANTS),
            "NC_M1b_by_variant": ncb, "F_M1": f_m1}


def scale_free_rank(delta_z: Mapping[str, float], exposed_set: Sequence[str]) -> dict[str, float]:
    """Per X drug: (count of unexposed drugs with delta_z >= x + 1) / (n_unexposed + 1)."""
    others = [v for d, v in delta_z.items() if d not in exposed_set and not math.isnan(v)]
    out = {}
    for d in exposed_set:
        v = delta_z[d]
        out[d] = math.nan if math.isnan(v) else (sum(o >= v for o in others) + 1) / (len(others) + 1)
    return out
