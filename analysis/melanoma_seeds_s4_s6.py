"""Stage-A sensitivity branches S4, S5 and S6 for protocol melanoma_R1_v0.3.

Rules: PROJECT_PLAN.md sections 1.3, 2.2 and 7 (frozen). Implementation details not
written in the frozen text are fixed in 04_melanoma_protocol/IMPLEMENTATION_NOTE_S4_S6.md,
written before this script was first run; its SHA-256 is recorded in the run manifest.

  S4  exclude metastatic samples whose TUMOR_TISSUE_SITE contains 'Primary' (gate: exactly
      two such samples among the primary-cohort selections); all other rules unchanged.
  S5  RPPA selection by per-feature OLS  z ~ group + AGE/10 + SEX + site proxy + center
      with HC3 standard errors (normal reference), BH over 198, q < 0.01 and |coef| >= 0.5;
      mutation component unchanged from the primary cohort.
  S6  primary tumours (expected 36 vs 28): description only; no P values, no BH, no seeds.

Reads no drug data, no STRING graph and no propagation output.

Usage:
    python -I -m analysis.melanoma_seeds_s4_s6 --data <raw_dir> --protocol <protocol_dir> --out <out_dir>
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import sys
from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import scipy
from scipy import stats

from analysis import melanoma_seeds as ms

ANALYSIS_VERSION = "melanoma-seeds-S4-S6-v1"
NOTE = "IMPLEMENTATION_NOTE_S4_S6.md"
SITE_LABELS = ("Distant Metastasis", "Regional Cutaneous or Subcutaneous Tissue", "Regional Lymph Node")
CENTER_MIN = 10
MIN_DF = 30
S6_EXPECTED = (36, 28)


# ---------------------------------------------------------------- HC3 OLS
def hc3_ols(X: np.ndarray, y: np.ndarray, coef_index: int = 1) -> dict[str, object]:
    """OLS with HC3 covariance. Returns coefficient, HC3 SE, normal and t P values.

    Status 'ok' or a not-estimable reason ('rank_deficient', 'df_below_30', 'leverage_one').
    """
    n, p = X.shape
    if np.linalg.matrix_rank(X) < p:
        return {"status": "rank_deficient"}
    df = n - p
    if df < MIN_DF:
        return {"status": "df_below_30", "df": df}
    xtx_inv = np.linalg.inv(X.T @ X)
    beta = xtx_inv @ X.T @ y
    resid = y - X @ beta
    h = np.einsum("ij,jk,ik->i", X, xtx_inv, X)
    if np.any(1.0 - h <= 1e-12):
        return {"status": "leverage_one", "df": df}
    omega = (resid / (1.0 - h)) ** 2
    cov = xtx_inv @ (X.T * omega) @ X @ xtx_inv
    se = math.sqrt(cov[coef_index, coef_index])
    b = float(beta[coef_index])
    if se == 0:
        return {"status": "zero_se", "df": df}
    t = b / se
    return {"status": "ok", "coef": b, "se_hc3": se, "df": df, "n": n, "k": p,
            "p_normal": float(2 * stats.norm.sf(abs(t))), "p_t": float(2 * stats.t.sf(abs(t), df)),
            "beta": beta, "cov": cov}


def dummies(values: Sequence[str], prefix: str) -> tuple[np.ndarray, list[str]]:
    """Treatment coding with the lexicographically first observed level as reference;
    levels absent from `values` produce no column."""
    levels = sorted(set(values))
    cols = [lv for lv in levels[1:]]
    if not cols:
        return np.zeros((len(values), 0)), []
    m = np.asarray([[1.0 if v == lv else 0.0 for lv in cols] for v in values])
    return m, [f"{prefix}={lv}" for lv in cols]


def site_proxy(tumor_tissue_site: str) -> str:
    tokens = {t.strip() for t in tumor_tissue_site.split("|")}
    hits = [lab for lab in SITE_LABELS if lab in tokens]
    return hits[0] if len(hits) == 1 else "other_or_ambiguous"


def center_map(patients: Sequence[str], center_of: Mapping[str, str]) -> dict[str, str]:
    counts = Counter(center_of.get(p, "") for p in patients)
    keep = {c for c, n in counts.items() if c and c not in ms.MISSING_TOKENS and n >= CENTER_MIN}
    return {p: (center_of.get(p, "") if center_of.get(p, "") in keep else "other_or_unknown") for p in patients}


# ------------------------------------------------------------------ S4
def site_conflicts(sample_rows: Sequence[Mapping[str, str]]) -> set[str]:
    out = set()
    for r in sample_rows:
        tokens = {t.strip() for t in r.get("TUMOR_TISSUE_SITE", "").split("|")}
        if r["SAMPLE_TYPE"].strip() == "Metastasis" and "Primary" in tokens:
            out.add(r["SAMPLE_ID"].strip())
    return out


# ------------------------------------------------------------------ S5
def adjusted_rppa(rppa_path: Path, cohort: Mapping[str, Mapping[str, str]], covariates: Mapping[str, Mapping[str, object]]
                  ) -> tuple[list[dict[str, object]], dict[str, dict[str, object]], dict[str, object]]:
    """cohort: {patient: {'sample_id', 'group'}}; covariates: {patient: {'age','sex','site','center'}}
    for patients with complete age and sex only."""
    with rppa_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        raw = [row for row in reader if row]
    if len(raw) != ms.RPPA_FAMILY_N:
        raise ValueError(f"expected {ms.RPPA_FAMILY_N} RPPA features, found {len(raw)}")
    column = {name: i for i, name in enumerate(header)}
    pats = sorted(p for p in cohort if p in covariates)
    n_case = sum(cohort[p]["group"] == "EK" for p in pats)
    n_ref = sum(cohort[p]["group"] == "REF" for p in pats)
    rows, pvals = [], []
    for record in raw:
        feature = record[0]
        genes = tuple(g.strip().upper() for g in feature.split("|", 1)[0].split() if g.strip())
        use, yv = [], []
        for p in pats:
            v = ms.parse_value(record[column[cohort[p]["sample_id"]]])
            if v is not None:
                use.append(p)
                yv.append(v)
        a = sum(cohort[p]["group"] == "EK" for p in use)
        b = len(use) - a
        reason, fit = "", None
        if not genes:
            reason = "no_gene_annotation"
        elif len(genes) > 1:
            reason = "multi_gene_annotation"
        elif a < ms.RPPA_MIN_N or b < ms.RPPA_MIN_N:
            reason = "fewer_than_50_finite_in_a_group"
        elif a < ms.RPPA_MIN_FRACTION * n_case or b < ms.RPPA_MIN_FRACTION * n_ref:
            reason = "finite_fraction_below_0.5_in_a_group"
        elif np.ptp(yv) == 0:
            reason = "constant_feature"
        else:
            grp = np.asarray([1.0 if cohort[p]["group"] == "EK" else 0.0 for p in use])
            age = np.asarray([float(covariates[p]["age"]) / 10.0 for p in use])
            parts = [np.ones(len(use)), grp, age]
            names = ["intercept", "group_EK", "age_per_10y"]
            for key in ("sex", "site", "center"):
                m, nm = dummies([str(covariates[p][key]) for p in use], key)
                for j in range(m.shape[1]):
                    if m[:, j].any():
                        parts.append(m[:, j])
                        names.append(nm[j])
            X = np.column_stack(parts)
            fit = hc3_ols(X, np.asarray(yv, float))
            if fit["status"] != "ok":
                reason = f"not_estimable_{fit['status']}"
        if reason:
            p_sel, coef, se, pt, df = 1.0, math.nan, math.nan, math.nan, (fit or {}).get("df", "")
        else:
            p_sel, coef, se, pt, df = fit["p_normal"], fit["coef"], fit["se_hc3"], fit["p_t"], fit["df"]
        rows.append({"feature": feature, "genes": ";".join(genes), "case_n": a, "ref_n": b,
                     "adjusted_group_coefficient": coef, "se_hc3": se, "residual_df": df,
                     "p_hc3_normal": p_sel, "p_hc3_t_descriptive": pt, "not_tested_reason": reason})
        pvals.append(p_sel)
    q = ms.benjamini_hochberg(pvals)
    selected: dict[str, dict[str, object]] = {}
    for row, qq in zip(rows, q, strict=True):
        row["q_bh_198"] = qq
        ok = (not row["not_tested_reason"]) and qq < ms.RPPA_Q and abs(row["adjusted_group_coefficient"]) >= ms.RPPA_MIN_EFFECT
        row["selected"] = "yes" if ok else "no"
        if not ok:
            continue
        effect = abs(row["adjusted_group_coefficient"])
        lineage = min(effect / ms.RPPA_EFFECT_SCALE, 1.0)
        gene = row["genes"]
        cur = selected.get(gene)
        if cur is None or lineage > cur["L"] or (lineage == cur["L"] and row["feature"] < cur["feature"]):
            selected[gene] = {"feature": row["feature"], "effect": effect, "q": qq,
                              "direction": "up" if row["adjusted_group_coefficient"] > 0 else "down", "L": lineage}
    meta = {"patients_complete_age_sex": len(pats), "EK": n_case, "REF": n_ref}
    return rows, selected, meta


# ------------------------------------------------------------------ S6
def describe_primary(rppa_path: Path, case_samples: Sequence[str], ref_samples: Sequence[str],
                     metastatic: Mapping[str, Mapping[str, object]]) -> list[dict[str, object]]:
    with rppa_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        raw = [row for row in reader if row]
    column = {name: i for i, name in enumerate(header)}
    out = []
    for record in raw:
        feature = record[0]
        a = [v for v in (ms.parse_value(record[column[s]]) for s in case_samples) if v is not None]
        b = [v for v in (ms.parse_value(record[column[s]]) for s in ref_samples) if v is not None]
        diff = float(np.median(a) - np.median(b)) if a and b else math.nan
        m = metastatic.get(feature, {})
        md = float(m.get("median_difference", "nan") or "nan")
        same = "NA" if math.isnan(diff) or math.isnan(md) or diff == 0 or md == 0 else ("yes" if (diff > 0) == (md > 0) else "no")
        out.append({"feature": feature, "case_n": len(a), "ref_n": len(b), "median_difference_primary": diff,
                    "meets_primary_testing_threshold": "yes" if min(len(a), len(b)) >= ms.RPPA_MIN_N else "no",
                    "median_difference_metastatic": md, "metastatic_selected": m.get("selected", ""),
                    "same_direction_as_metastatic": same})
    return out


# ------------------------------------------------------------------ main
def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--protocol", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    frozen = ms.check_stage_a_freeze(args.protocol)
    amendments = ms.check_amendments(args.protocol, ms.REQUIRED_AMENDMENTS)
    note_path = args.protocol / NOTE
    if not note_path.exists():
        raise RuntimeError("implementation note missing")
    out = args.out
    if out.exists() and any(out.iterdir()):
        raise RuntimeError(f"refusing to overwrite non-empty {out}")
    data = args.data
    with (data / "data_rppa_zscores.txt").open(encoding="utf-8-sig") as handle:
        rppa_header = handle.readline().rstrip("\n").split("\t")
    samples = ms.read_cbio_tsv(data / "data_clinical_sample.txt")
    patients = ms.read_cbio_tsv(data / "data_clinical_patient.txt")
    inputs = {"samples": samples, "patients": patients,
              "sequenced": ms.read_case_list(data / "cases_sequenced.txt"),
              "rppa_samples": set(rppa_header[1:]),
              "braf": ms.braf_records_by_sample(data / "data_mutations.txt")}
    manifest: dict[str, object] = {
        "analysis_version": ANALYSIS_VERSION, "protocol_id": ms.PROTOCOL_ID,
        "frozen_files_verified": frozen, "amendments_verified": amendments,
        "implementation_note_sha256": ms.sha256(note_path),
        "inputs_sha256": {f: ms.sha256(data / f) for f in
                          ["data_clinical_patient.txt", "data_clinical_sample.txt", "data_mutations.txt",
                           "data_rppa_zscores.txt", "cases_sequenced.txt"]},
        "code_sha256": {Path(__file__).name: ms.sha256(Path(__file__)), Path(ms.__file__).name: ms.sha256(Path(ms.__file__))},
        "python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__,
        "platform": platform.platform()}
    primary = ms.run_branch("primary_metastasis", data, inputs, "Metastasis")
    pc = primary["counts"]
    if (pc["EK_patients"], pc["REF_patients"]) != (106, 125):
        raise RuntimeError(f"primary cohort not reproduced: {pc}")
    out.mkdir(parents=True, exist_ok=True)

    # ---- S4
    selected_primary = {v["sample_id"] for v in primary["chosen"].values()}
    conflicts = site_conflicts(samples)
    in_primary = sorted(conflicts & selected_primary)
    site_of = {r["SAMPLE_ID"].strip(): r.get("TUMOR_TISSUE_SITE", "") for r in samples}
    manifest["S4_conflict_samples_in_primary"] = [{"sample_id": s, "tumor_tissue_site": site_of[s]} for s in in_primary]
    manifest["S4_conflict_samples_all"] = sorted(conflicts)
    if len(in_primary) != 2:
        manifest["status"] = "STOPPED_S4_CONFLICT_COUNT"
        (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
        print(json.dumps(manifest["S4_conflict_samples_in_primary"], indent=2))
        return 2
    s4_inputs = dict(inputs, samples=[r for r in samples if r["SAMPLE_ID"].strip() not in conflicts])
    s4 = ms.run_branch("S4_site_conflict_excluded", data, s4_inputs, "Metastasis")
    d4 = out / s4["name"]
    ms.write_csv(d4 / "cohort_selected_patients.csv", [{"patient_id": p, **v} for p, v in sorted(s4["chosen"].items())])
    ms.write_csv(d4 / "mutation_statistics.csv", s4["mutation"])
    ms.write_csv(d4 / "rppa_feature_tests.csv", s4["rppa"])
    ms.write_csv(d4 / "seed_candidates.csv", s4["seeds"])
    w4 = ms.arm_weights(s4["seeds"], "M0") + ms.arm_weights(s4["seeds"], "M1")
    if w4:
        ms.write_csv(d4 / "seeds_frozen.csv", w4)

    # ---- S5
    pinfo = {r["PATIENT_ID"].strip(): r for r in patients}
    sinfo = {r["SAMPLE_ID"].strip(): r for r in samples}
    cov_all, dropped = {}, []
    for p, v in primary["chosen"].items():
        age = ms.parse_value(pinfo.get(p, {}).get("AGE", ""))
        sex = pinfo.get(p, {}).get("SEX", "").strip()
        if age is None or sex in ms.MISSING_TOKENS:
            dropped.append(p)
            continue
        s = sinfo[v["sample_id"]]
        cov_all[p] = {"age": age, "sex": sex, "site": site_proxy(s.get("TUMOR_TISSUE_SITE", "")),
                      "center_raw": s.get("TISSUE_SOURCE_SITE", "").strip()}
    cmap = center_map(sorted(cov_all), {p: c["center_raw"] for p, c in cov_all.items()})
    for p in cov_all:
        cov_all[p]["center"] = cmap[p]
    rows5, sel5, meta5 = adjusted_rppa(data / "data_rppa_zscores.txt", primary["chosen"], cov_all)
    seeds5 = ms.build_seed_rows(primary["mutation"], sel5)
    for r in seeds5:
        r["rppa_abs_adjusted_coefficient"] = r.pop("rppa_abs_median_difference")
    d5 = out / "S5_hc3_adjusted"
    ms.write_csv(d5 / "rppa_feature_tests_adjusted.csv", rows5)
    ms.write_csv(d5 / "covariates.csv", [{"patient_id": p, "group": primary["chosen"][p]["group"],
                                          **{k: c[k] for k in ("age", "sex", "site", "center", "center_raw")}}
                                         for p, c in sorted(cov_all.items())])
    ms.write_csv(d5 / "seed_candidates.csv", seeds5)
    w5 = ms.arm_weights(seeds5, "M0") + ms.arm_weights(seeds5, "M1")
    if w5:
        ms.write_csv(d5 / "seeds_frozen.csv", w5)
    s5_counts = {**meta5, "dropped_missing_age_or_sex": sorted(dropped),
                 "center_levels": dict(Counter(c["center"] for c in cov_all.values())),
                 "site_levels": dict(Counter(c["site"] for c in cov_all.values())),
                 "rppa_features_tested": sum(not r["not_tested_reason"] for r in rows5),
                 "not_tested_reasons": dict(Counter(r["not_tested_reason"] for r in rows5 if r["not_tested_reason"])),
                 "rppa_features_selected": sum(r["selected"] == "yes" for r in rows5),
                 "selected_t_reference_descriptive": sum(
                     1 for r, qq in zip(rows5, ms.benjamini_hochberg(
                         [1.0 if r["not_tested_reason"] else r["p_hc3_t_descriptive"] for r in rows5]))
                     if not r["not_tested_reason"] and qq < ms.RPPA_Q and abs(r["adjusted_group_coefficient"]) >= ms.RPPA_MIN_EFFECT),
                 "seeds_M0": sum(r["include_M0"] == "yes" for r in seeds5),
                 "seeds_M1": sum(r["include_M1"] == "yes" for r in seeds5),
                 "selected_genes": sorted(sel5)}

    # ---- S6
    s6 = ms.build_cohort(samples, patients, inputs["sequenced"], inputs["rppa_samples"], inputs["braf"], "Primary")[0]
    case6 = sorted(v["sample_id"] for v in s6.values() if v["group"] == "EK")
    ref6 = sorted(v["sample_id"] for v in s6.values() if v["group"] == "REF")
    manifest["S6_counts"] = {"EK": len(case6), "REF": len(ref6), "expected": list(S6_EXPECTED)}
    if (len(case6), len(ref6)) != S6_EXPECTED:
        manifest["status"] = "STOPPED_S6_COUNT_MISMATCH"
        (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
        return 3
    d6 = out / "S6_primary_descriptive"
    met = {r["feature"]: r for r in primary["rppa"]}
    desc = describe_primary(data / "data_rppa_zscores.txt", case6, ref6, met)
    ms.write_csv(d6 / "rppa_description.csv", desc)
    ek6 = {v["sample_id"]: p for p, v in s6.items() if v["group"] == "EK"}
    mut6 = [{k: r[k] for k in ("gene", "mutated_patients", "denominator", "fraction")}
            for r in ms.mutation_statistics(data / "data_mutations.txt", ek6)]
    ms.write_csv(d6 / "mutation_description.csv", mut6)
    comparable = [r for r in desc if r["same_direction_as_metastatic"] in ("yes", "no")]
    manifest["S6_summary"] = {"features": len(desc),
                              "features_meeting_testing_threshold": sum(r["meets_primary_testing_threshold"] == "yes" for r in desc),
                              "direction_comparable": len(comparable),
                              "same_direction": sum(r["same_direction_as_metastatic"] == "yes" for r in comparable),
                              "note": "descriptive only; no P values, no BH, no seeds"}
    manifest["branch_counts"] = {"S4_site_conflict_excluded": s4["counts"], "S5_hc3_adjusted": s5_counts}
    manifest["status"] = "COMPLETED"
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    print(json.dumps({k: manifest[k] for k in ("S4_conflict_samples_in_primary", "branch_counts", "S6_counts", "S6_summary")},
                     indent=2, ensure_ascii=False, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
