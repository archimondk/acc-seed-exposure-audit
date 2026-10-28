"""Reference tests for the S4-S6 branch implementation (IMPLEMENTATION_NOTE_S4_S6.md)."""
import csv
import math

import numpy as np
import pytest

from analysis import melanoma_seeds_s4_s6 as s


def _design(n, rng):
    x1 = rng.integers(0, 2, n).astype(float)
    x2 = rng.normal(size=n)
    cat = rng.integers(0, 3, n)
    X = np.column_stack([np.ones(n), x1, x2, (cat == 1).astype(float), (cat == 2).astype(float)])
    y = 0.7 * x1 + 0.3 * x2 + rng.standard_t(4, n) * (1 + x1)
    return X, y


def test_hc3_equals_sum_of_leave_one_out_shifts():
    """V_HC3 = sum_i (b - b_(-i))(b - b_(-i))^T, an identity checked by actual refits."""
    rng = np.random.default_rng(1)
    X, y = _design(60, rng)
    fit = s.hc3_ols(X, y)
    assert fit["status"] == "ok"
    b = fit["beta"]
    acc = np.zeros((X.shape[1], X.shape[1]))
    for i in range(len(y)):
        keep = np.arange(len(y)) != i
        bi = np.linalg.lstsq(X[keep], y[keep], rcond=None)[0]
        acc += np.outer(b - bi, b - bi)
    np.testing.assert_allclose(fit["cov"], acc, rtol=1e-10, atol=1e-14)
    assert fit["se_hc3"] == pytest.approx(math.sqrt(acc[1, 1]), rel=1e-10)


def test_hc3_p_values_use_normal_reference_and_t_descriptive():
    rng = np.random.default_rng(2)
    X, y = _design(80, rng)
    fit = s.hc3_ols(X, y)
    z = fit["coef"] / fit["se_hc3"]
    from scipy import stats
    assert fit["p_normal"] == pytest.approx(2 * stats.norm.sf(abs(z)), rel=1e-12)
    assert fit["p_t"] == pytest.approx(2 * stats.t.sf(abs(z), 80 - 5), rel=1e-12)
    assert fit["df"] == 75


def test_hc3_matches_statsmodels_when_available():
    sm = pytest.importorskip("statsmodels.api")
    rng = np.random.default_rng(3)
    X, y = _design(70, rng)
    ours = s.hc3_ols(X, y)
    ref = sm.OLS(y, X).fit(cov_type="HC3")
    assert ours["se_hc3"] == pytest.approx(ref.bse[1], rel=1e-10)
    assert ours["p_normal"] == pytest.approx(ref.pvalues[1], rel=1e-8)


def test_hc3_not_estimable_cases():
    rng = np.random.default_rng(4)
    X, y = _design(60, rng)
    assert s.hc3_ols(np.column_stack([X, X[:, 1]]), y)["status"] == "rank_deficient"
    assert s.hc3_ols(X[:34], y[:34])["status"] == "df_below_30"
    Xl = np.column_stack([X, np.eye(60)[:, 0]])  # dummy observed in one case only -> leverage 1
    assert s.hc3_ols(Xl, y)["status"] == "leverage_one"


def test_site_proxy_and_conflicts():
    assert s.site_proxy("Extremities|Regional Lymph Node") == "Regional Lymph Node"
    assert s.site_proxy("Distant Metastasis|Trunk") == "Distant Metastasis"
    assert s.site_proxy("Regional Lymph Node|Distant Metastasis") == "other_or_ambiguous"
    assert s.site_proxy("Other") == "other_or_ambiguous"
    assert s.site_proxy("NA|Trunk") == "other_or_ambiguous"
    rows = [{"SAMPLE_ID": "a", "SAMPLE_TYPE": "Metastasis", "TUMOR_TISSUE_SITE": "Extremities|Primary"},
            {"SAMPLE_ID": "b", "SAMPLE_TYPE": "Primary", "TUMOR_TISSUE_SITE": "Primary|Trunk"},
            {"SAMPLE_ID": "c", "SAMPLE_TYPE": "Metastasis", "TUMOR_TISSUE_SITE": "Distant Metastasis"},
            {"SAMPLE_ID": "d", "SAMPLE_TYPE": "Metastasis", "TUMOR_TISSUE_SITE": "Primary"}]
    assert s.site_conflicts(rows) == {"a", "d"}


def test_center_map_threshold_and_missing():
    pats = [f"p{i}" for i in range(25)]
    center = {p: ("A" if i < 10 else "B" if i < 19 else "") for i, p in enumerate(pats)}
    m = s.center_map(pats, center)
    assert {m[p] for p in pats[:10]} == {"A"}
    assert {m[p] for p in pats[10:]} == {"other_or_unknown"}


def test_dummies_reference_is_first_level_and_absent_levels_dropped():
    m, names = s.dummies(["b", "a", "c", "a"], "x")
    assert names == ["x=b", "x=c"]
    assert m.tolist() == [[1, 0], [0, 0], [0, 1], [0, 0]]
    m2, names2 = s.dummies(["z", "z"], "x")
    assert m2.shape == (2, 0) and names2 == []


def _toy_rppa(path, samples, signal_rows):
    rng = np.random.default_rng(5)
    with path.open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t")
        w.writerow(["Composite.Element.REF"] + samples)
        for i in range(198):
            if i == 0:
                name = "GENEA|GENEA_pS1"
            elif i == 1:
                name = "GENEB GENEC|multi"
            else:
                name = f"G{i}|G{i}_ab"
            vals = rng.normal(size=len(samples))
            if i in signal_rows:
                vals = vals + np.asarray([2.0 if sm.startswith("E") else 0.0 for sm in samples])
            w.writerow([name] + [f"{v:.6f}" for v in vals])


def test_adjusted_rppa_selection_and_placeholders(tmp_path):
    ek = [f"E{i}" for i in range(70)]
    ref = [f"R{i}" for i in range(70)]
    samples = ek + ref
    _toy_rppa(tmp_path / "rppa.txt", samples, signal_rows={0, 1})
    cohort = {f"P{s_}": {"sample_id": s_, "group": "EK" if s_.startswith("E") else "REF"} for s_ in samples}
    rng = np.random.default_rng(6)
    cov = {p: {"age": float(rng.integers(30, 80)), "sex": "Male" if i % 2 else "Female",
               "site": ["Regional Lymph Node", "Distant Metastasis", "other_or_ambiguous"][i % 3],
               "center": ["A", "B", "other_or_unknown"][(i // 3) % 3]} for i, p in enumerate(sorted(cohort))}
    rows, selected, meta = s.adjusted_rppa(tmp_path / "rppa.txt", cohort, cov)
    assert meta == {"patients_complete_age_sex": 140, "EK": 70, "REF": 70}
    by = {r["feature"]: r for r in rows}
    assert by["GENEB GENEC|multi"]["not_tested_reason"] == "multi_gene_annotation"
    assert by["GENEB GENEC|multi"]["p_hc3_normal"] == 1.0
    assert set(selected) == {"GENEA"}
    assert selected["GENEA"]["L"] == pytest.approx(min(abs(by["GENEA|GENEA_pS1"]["adjusted_group_coefficient"]) / 2, 1))
    assert len(rows) == 198


def test_adjusted_rppa_drops_patients_without_covariates(tmp_path):
    ek = [f"E{i}" for i in range(70)]
    ref = [f"R{i}" for i in range(70)]
    _toy_rppa(tmp_path / "rppa.txt", ek + ref, signal_rows=set())
    cohort = {f"P{x}": {"sample_id": x, "group": "EK" if x.startswith("E") else "REF"} for x in ek + ref}
    cov = {p: {"age": 50.0, "sex": "Male", "site": "Regional Lymph Node", "center": "A"}
           for p in sorted(cohort)[:100]}
    rows, _, meta = s.adjusted_rppa(tmp_path / "rppa.txt", cohort, cov)
    assert meta["patients_complete_age_sex"] == 100
    # one reference group shrinks below 50 -> nothing testable
    assert all(r["not_tested_reason"] for r in rows)


def test_adjusted_rppa_collinear_covariates_are_not_estimable(tmp_path):
    """Frozen rule: only unobserved levels / all-zero dummies are dropped; remaining
    rank deficiency is recorded as not estimable with P = 1, never repaired."""
    ek = [f"E{i}" for i in range(70)]
    ref = [f"R{i}" for i in range(70)]
    _toy_rppa(tmp_path / "rppa.txt", ek + ref, signal_rows={0})
    cohort = {f"P{x}": {"sample_id": x, "group": "EK" if x.startswith("E") else "REF"} for x in ek + ref}
    cov = {p: {"age": 40.0 + i % 7, "sex": "Male" if i % 2 else "Female",
               "site": ["s1", "s2", "s3"][i % 3], "center": ["c1", "c2", "c3"][i % 3]}
           for i, p in enumerate(sorted(cohort))}
    rows, selected, _ = s.adjusted_rppa(tmp_path / "rppa.txt", cohort, cov)
    tested = [r for r in rows if r["not_tested_reason"] != "multi_gene_annotation"]
    assert all(r["not_tested_reason"] == "not_estimable_rank_deficient" for r in tested)
    assert all(r["p_hc3_normal"] == 1.0 for r in rows) and not selected
