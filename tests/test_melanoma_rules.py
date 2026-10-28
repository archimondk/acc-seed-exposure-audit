"""Gold-standard tests for the frozen H4/H5/H6 rules (melanoma_R1_v0.3)."""
import math
import random
from fractions import Fraction

import numpy as np
import pytest

from analysis import melanoma_rules as r
from analysis.mechanism_enrichment import benjamini_hochberg as repo_bh

V = r.VARIANTS
SIX = ["Dabrafenib", "Vemurafenib", "Encorafenib", "Trametinib", "Cobimetinib", "Binimetinib"]
BRAFI, MEKI = SIX[:3], SIX[3:]
ANCH = {"ACC": {"c": Fraction(46, 108), "rho": 0.8220085221393144, "u": Fraction(0, 1)},
        "Breast": {"c": Fraction(14, 108), "rho": 0.5772974445641333, "u": Fraction(13, 94)}}


# ------------------------------------------------------------ statistics
def test_bh_matches_repository_implementation():
    rng = random.Random(1)
    for _ in range(50):
        p = [rng.random() for _ in range(rng.randint(1, 40))]
        assert r.benjamini_hochberg(p) == pytest.approx(repo_bh(p))


def test_empirical_upper_counts_ties_and_plus_one():
    b, p = r.empirical_upper(0.5, np.array([0.1, 0.5, 0.7, 0.2]))
    assert b == 2 and p == pytest.approx(3 / 5)
    b, p = r.empirical_upper(0.99, np.zeros(10000))
    assert b == 0 and p == pytest.approx(1 / 10001)


def test_empirical_upper_float64_not_rounded():
    obs = 0.5 + 1e-9  # float32 would merge it with 0.5
    b, _ = r.empirical_upper(obs, np.array([0.5], dtype=np.float64))
    assert b == 0


def test_clopper_pearson_known_values():
    lo, hi = r.clopper_pearson(0, 10000)
    assert lo == 0 and hi == pytest.approx(1 - 0.025 ** (1 / 10000), rel=1e-6)
    lo, hi = r.clopper_pearson(10, 10)
    assert hi == 1 and lo == pytest.approx(0.025 ** (1 / 10), rel=1e-6)


def test_z_score_zero_sd_is_nan():
    assert math.isnan(r.z_score(0.3, np.full(100, 0.2)))
    assert r.z_score(1.0, np.array([0.0, 2.0])) == pytest.approx(0.0)


def test_top_k_breaks_ties_by_name():
    assert r.top_k(["b", "a", "c"], [1.0, 1.0, 0.5], 2) == ["a", "b"]
    assert r.ordinal_ranks(["b", "a"], [1.0, 1.0]) == {"a": 1, "b": 2}


# -------------------------------------------------------------------- H4
def group(p_primary=0.01, others=(0.02, 0.03, 0.04), above=(True, True, True, True)):
    ps = [p_primary, *others]
    return {v: {"p": p, "observed": 0.6 if a else 0.4, "null_mean": 0.5} for v, p, a in zip(V, ps, above)}


def test_h4_pass():
    flags = {d: d in ("Dabrafenib", "Trametinib", "Cobimetinib") for d in SIX}
    out = r.evaluate_h4(group(), flags, BRAFI, MEKI)
    assert out["status"] == "pass"


def test_h4_q_uses_primary_variant_after_bh():
    # primary P = 0.03 with others 0.5 -> q_primary = 0.03*4/1 = 0.12 -> fails q only
    flags = {d: d in ("Dabrafenib", "Trametinib", "Cobimetinib") for d in SIX}
    out = r.evaluate_h4(group(0.03, (0.5, 0.6, 0.7)), flags, BRAFI, MEKI)
    assert out["criteria"]["primary_group_p_lt_0_05"] is True
    assert out["criteria"]["primary_group_q_lt_0_05"] is False
    assert out["q_by_variant"][r.PRIMARY] == pytest.approx(0.12)
    assert out["status"] == "partial_recovery"


def test_h4_requires_both_classes():
    flags = {d: d in MEKI for d in SIX}  # three MEK inhibitors only
    out = r.evaluate_h4(group(), flags, BRAFI, MEKI)
    assert out["criteria"]["at_least_3_top25_with_one_BRAFi_and_one_MEKi"] is False


def test_h4_fail_label():
    flags = {d: False for d in SIX}
    out = r.evaluate_h4(group(0.5, (0.6, 0.7, 0.8), (False,) * 4), flags, BRAFI, MEKI)
    assert out["status"] == "fail"


# -------------------------------------------------------------------- H5
@pytest.mark.parametrize("c,rho,u,branch,verdict", [
    (Fraction(30, 108), 0.70, Fraction(5, 100), "between_anchors", "SUPPORTED"),
    (Fraction(30, 108), 0.90, Fraction(5, 100), "between_anchors", "PARTIAL"),
    (Fraction(30, 108), 0.90, Fraction(20, 100), "between_anchors", "NOT_SUPPORTED"),
    (Fraction(60, 108), 0.85, Fraction(0, 1), "above_ACC", "SUPPORTED"),
    (Fraction(5, 108), 0.30, Fraction(22, 103), "below_Breast", "SUPPORTED"),
    (Fraction(5, 108), 0.60, Fraction(22, 103), "below_Breast", "PARTIAL"),
    (Fraction(45, 108), 0.70, Fraction(1, 10), "near_ACC_one_sided", "SUPPORTED"),
    (Fraction(15, 108), 0.70, Fraction(1, 10), "near_Breast_one_sided", "SUPPORTED"),
])
def test_h5_branches(c, rho, u, branch, verdict):
    out = r.evaluate_h5(c, rho, u, ANCH)
    assert out["branch"] == branch and out["verdict"] == verdict


def test_h5_endpoint_equality_is_consistent():
    out = r.evaluate_h5(Fraction(5, 108), ANCH["Breast"]["rho"], Fraction(13, 94), ANCH)
    assert out["verdict"] == "SUPPORTED"


def test_h5_missing_and_degenerate():
    assert r.evaluate_h5(Fraction(5, 108), math.nan, Fraction(22, 103), ANCH)["verdict"] == "INCOMPLETE"
    assert r.evaluate_h5(Fraction(5, 108), math.nan, None, ANCH)["verdict"] == "NOT_EVALUABLE"
    assert r.evaluate_h5(Fraction(0, 108), 0.1, Fraction(1, 4), ANCH)["verdict"] == "NOT_EVALUABLE"
    assert r.evaluate_h5(Fraction(1, 1), 0.1, None, ANCH)["verdict"] == "NOT_EVALUABLE"


def test_h5_tie_band_uses_exact_fraction():
    # |c - c_breast| exactly 0.02 is NOT within the band (strict <)
    c = Fraction(14, 108) + Fraction(1, 50)
    out = r.evaluate_h5(c, 0.6, Fraction(1, 10), ANCH)
    assert out["branch"] == "between_anchors"


def test_h5_arithmetic_bound_detected_for_melanoma_structure():
    lvals = [0.2, 0.1, 0.1, 0.05, 0.03] + [0.0] * 103
    bounds = r.h5_feasible_bounds(lvals, 27)
    assert bounds["u_min"] == Fraction(22, 103) and bounds["u_max"] == Fraction(27, 103)
    assert 0.35 < bounds["rho_max"] < 0.38
    out = r.evaluate_h5(Fraction(5, 108), 0.2, Fraction(24, 103), ANCH, bounds)
    assert out["verdict"] == "SUPPORTED" and out["annotation"] == "ARITHMETICALLY_DETERMINED"


def test_h5_no_annotation_when_not_forced():
    lvals = [0.1] * 46 + [0.0] * 62
    bounds = r.h5_feasible_bounds(lvals, 27)
    out = r.evaluate_h5(Fraction(46, 108), 0.82, Fraction(0, 1), ANCH, bounds)
    assert out["annotation"] == ""


# -------------------------------------------------------------------- H6
def dz_table(x_vals, six_val=0.1, other_val=0.2, n_other=20, na=None):
    na = na or {}
    t = {}
    for v in V:
        row = {"Dacarbazine": x_vals[v][0], "Tioguanine": x_vals[v][1]}
        row.update({d: six_val for d in SIX})
        row.update({f"D{i}": other_val for i in range(n_other)})
        for d in na.get(v, []):
            row[d] = math.nan
        t[v] = row
    return t


X = ["Dacarbazine", "Tioguanine"]
GOOD = {v: (2.5, 0.4) for v in V}


def test_h6_supported():
    out = r.evaluate_h6(dz_table(GOOD), X, SIX)
    assert out["verdict"] == "SUPPORTED" and out["I2_count"] == 4 and out["NC_ROBUST"] is True


def test_h6_partial_when_nc_fails_in_primary():
    out = r.evaluate_h6(dz_table(GOOD, six_val=0.8), X, SIX)
    assert out["verdict"] == "PARTIAL" and out["NC_M1a_by_variant"][r.PRIMARY] is False


def test_h6_nc_only_primary_counts_for_verdict():
    t = dz_table(GOOD)
    for v in V[1:]:
        t[v]["D0"] = 3.0  # NC-M1b fails outside primary
    out = r.evaluate_h6(t, X, SIX)
    assert out["verdict"] == "SUPPORTED"
    assert all(out["NC_M1b_by_variant"][v]["status"] is False for v in V[1:])


def test_h6_partial_when_i2_fails():
    vals = dict(GOOD)
    for v in V[1:3]:
        vals[v] = (0.5, 0.4)
    out = r.evaluate_h6(dz_table(vals), X, SIX)
    assert out["I2_count"] == 2 and out["verdict"] == "PARTIAL"


def test_h6_not_supported_and_falsification():
    vals = {v: (-0.3, -0.1) for v in V}
    out = r.evaluate_h6(dz_table(vals), X, SIX)
    assert out["verdict"] == "NOT_SUPPORTED" and out["F_M1"] is True


def test_h6_requires_all_x_positive():
    vals = {v: (2.5, -0.1) for v in V}
    assert r.evaluate_h6(dz_table(vals), X, SIX)["verdict"] == "NOT_SUPPORTED"


def test_h6_na_rules():
    vals = {v: (2.5, 0.4) for v in V}
    vals[r.PRIMARY] = (math.nan, 0.4)
    assert r.evaluate_h6(dz_table(vals), X, SIX)["verdict"] == "NOT_EVALUABLE"
    vals = dict(GOOD)
    vals[V[1]] = (math.nan, 0.4)
    vals[V[2]] = (0.1, 0.4)
    out = r.evaluate_h6(dz_table(vals), X, SIX)
    assert out["I1_NA_variants"] == 1 and out["I2_count"] == 2 and out["verdict"] == "PARTIAL"
    out = r.evaluate_h6(dz_table(GOOD, na={r.PRIMARY: ["Trametinib"]}), X, SIX)
    assert out["NC_M1a_by_variant"][r.PRIMARY] == "NA" and out["verdict"] == "PARTIAL"
    out = r.evaluate_h6(dz_table(GOOD, n_other=20, na={r.PRIMARY: ["D0", "D1"]}), X, SIX)
    # 2 NA of 26 non-X drugs > 5% -> NA
    assert out["NC_M1b_by_variant"][r.PRIMARY]["status"] == "NA" and out["verdict"] == "PARTIAL"
    out = r.evaluate_h6(dz_table(GOOD, n_other=60, na={r.PRIMARY: ["D0"]}), X, SIX)
    assert out["NC_M1b_by_variant"][r.PRIMARY]["status"] is True and out["NC_M1b_by_variant"][r.PRIMARY]["n_na"] == 1


def test_h6_empty_x():
    assert r.evaluate_h6(dz_table(GOOD), [], SIX)["verdict"] == "NOT_EVALUABLE"


def test_scale_free_rank():
    out = r.scale_free_rank({"X1": 1.0, "a": 0.5, "b": 1.5, "c": math.nan}, ["X1"])
    assert out["X1"] == pytest.approx((1 + 1) / (2 + 1))
