#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Amendment 1 re-analysis, v2 -- auto-discovers the arm result files instead of
assuming a fixed directory layout.

Pass 1  inventory every CSV under leakage_audit/ that looks like a drug-rank
        table, and try to map each to an arm (A1/A2/B1/B2/B2_lo/B2_hi)
Pass 2  if the mapping succeeds, compute the scale-free intervention effect
        within each normalization variant

WHY
The frozen L1/L2 rules used an ABSOLUTE delta-z threshold calibrated on the
primary column-min-max scale. Rank-percentile gene scaling compresses the
r_ACC dynamic range and therefore compresses every effect magnitude, so the
same absolute cut is not comparable across variants. This computes, within
each variant, how extreme the abemaciclib shift is relative to that variant's
own unexposed-drug shift distribution.

Post-hoc. Report as such under protocol amendment 1.

Run:  python variant_scale_free_effect_v2.py
"""

import glob
import argparse
import json
import os
import re

import numpy as np
import pandas as pd

_root_parser = argparse.ArgumentParser(add_help=False)
_root_parser.add_argument("--project-root", default=os.getcwd())
_root_args, _unknown_args = _root_parser.parse_known_args()
DATA_DIR = os.path.abspath(_root_args.project_root)
AUDIT_DIR = os.path.join(DATA_DIR, "leakage_audit")
OUT_DIR = AUDIT_DIR

VARIANTS = ["column_minmax", "column_gene_rank",
            "uniform_ratio_gene_rank", "symmetric_gene_rank"]
PRIMARY_DRUG = "abemaciclib"
NEG_CONTROL = "ribociclib"
MANIPULATED = "RB1"

# checked in this order; first match wins
ARM_PATTERNS = [
    ("B2_lo", [r"b2[_\-]?lo"]),
    ("B2_hi", [r"b2[_\-]?hi"]),
    ("A1", [r"\ba1\b", r"a1[_\-]", r"acc[_\-]ref"]),
    ("A2", [r"\ba2\b", r"a2[_\-]", r"minus[_\-]?rb1", r"without[_\-]?rb1"]),
    ("B1", [r"\bb1\b", r"b1[_\-]", r"breast[_\-]ref"]),
    ("B2", [r"\bb2\b", r"b2[_\-]", r"plus[_\-]?rb1", r"with[_\-]?rb1"]),
]

pd.set_option("display.width", 230)
pd.set_option("display.max_columns", 40)
pd.set_option("display.max_rows", 200)


def read_any(p, **kw):
    for enc in ("utf-8-sig", "utf-8", "gbk", "latin-1"):
        try:
            return pd.read_csv(p, encoding=enc, **kw)
        except UnicodeDecodeError:
            continue
    return pd.read_csv(p, encoding="latin-1", **kw)


def find_cols(df):
    """Return (drug_col, z_col, variant_col_or_None, arm_col_or_None)."""
    low = {str(c).strip().lower(): c for c in df.columns}
    dcol = low.get("drug") or low.get("drug_name") or low.get("compound")
    zcol = next((c for k, c in low.items()
                 if k.startswith("z_") or k == "z"), None)
    vcol = low.get("variant") or low.get("normalization") or low.get("norm")
    acol = low.get("arm") or low.get("arm_id")
    return dcol, zcol, vcol, acol


def arm_of(path):
    p = path.lower().replace("\\", "/")
    for arm, pats in ARM_PATTERNS:
        if any(re.search(pat, p) for pat in pats):
            return arm
    return None


# ================================================================ PASS 1
print("\n" + "=" * 80)
print("PASS 1  inventory of drug-rank tables under leakage_audit/")
print("=" * 80)

if not os.path.isdir(AUDIT_DIR):
    raise SystemExit(f"  no leakage_audit directory at {AUDIT_DIR}")

frames = {}          # arm -> dataframe with _drug/_z/variant
combined = []        # tables that carry their own arm column

for p in sorted(glob.glob(os.path.join(AUDIT_DIR, "**", "*.csv"), recursive=True)):
    try:
        head = read_any(p, nrows=3)
    except Exception:
        continue
    dcol, zcol, vcol, acol = find_cols(head)
    if not (dcol and zcol):
        continue

    df = read_any(p)
    dcol, zcol, vcol, acol = find_cols(df)
    rel = os.path.relpath(p, DATA_DIR)
    arm = arm_of(p)
    n_v = df[vcol].nunique() if vcol else 1
    print(f"\n  {rel}")
    print(f"      rows={len(df)}  drug='{dcol}'  z='{zcol}'  "
          f"variant='{vcol}'  arm_col='{acol}'  n_variants={n_v}")
    print(f"      inferred arm from path: {arm or '(none)'}")

    df = df.copy()
    df["_drug"] = df[dcol].astype(str).str.strip().str.lower()
    df["_z"] = pd.to_numeric(df[zcol], errors="coerce")
    df["_variant"] = df[vcol].astype(str) if vcol else VARIANTS[0]

    if acol:
        df["_arm"] = df[acol].astype(str)
        combined.append(df)
        print("      -> table carries its own arm column; will be split")
    elif arm and arm not in frames:
        frames[arm] = df
        print(f"      -> assigned to arm {arm}")

# split any combined table
for df in combined:
    for a, sub in df.groupby("_arm"):
        key = arm_of(a) or a
        if key not in frames:
            frames[key] = sub
            print(f"  split: arm column value '{a}' -> {key} ({len(sub)} rows)")

print("\n" + "-" * 80)
print(f"  arms resolved: {sorted(frames)}")
need = {"A1", "A2", "B1", "B2"}
missing = need - set(frames)
if missing:
    print(f"  MISSING: {sorted(missing)}")
    print("\n  Fix one of:")
    print("    * rename the arm folders/files so the path contains A1/A2/B1/B2, or")
    print("    * add an 'arm' column to a combined drug-rank table, or")
    print("    * tell me the actual paths and I'll hard-map them.")
    raise SystemExit(0)

# ================================================================ PASS 2
print("\n" + "=" * 80)
print("PASS 2  scale-free intervention effect")
print("=" * 80)

# exposure: which drugs are directly associated with the manipulated seed
ep = [p for p in glob.glob(os.path.join(DATA_DIR, "**", "bindex_edges_1304.csv"),
                           recursive=True)
      if "repro_outputs" not in p.lower() and "tmp" not in p.lower()]
e = read_any(sorted(ep, key=len)[0])
dc = next((c for c in e.columns if str(c).strip().lower() in
           ("drug", "drug_name", "compound")), e.columns[0])
gc = next(c for c in e.columns
          if e[c].astype(str).str.fullmatch(r"[A-Z0-9\-]{2,15}").mean() > 0.8)
e = pd.DataFrame({"drug": e[dc].astype(str).str.strip().str.lower(),
                  "gene": e[gc].astype(str).str.strip().str.upper()}).drop_duplicates()
exposed_all = set(e.loc[e["gene"] == MANIPULATED, "drug"])

universe = set(frames["A1"]["_drug"])
exposed = exposed_all & universe
print(f"\n  {MANIPULATED}-associated drugs in the network : {sorted(exposed_all)}")
print(f"  of which inside the locked universe        : {sorted(exposed)}")
outside = exposed_all - universe
if outside:
    print(f"  NOT in the locked universe                 : {sorted(outside)}")
print(f"  universe size                              : {len(universe)}")

rows, curves = [], []
PAIRS = [("ACC", "A1", "A2"), ("Breast", "B2", "B1")]

for disease, hi, lo in PAIRS:
    H, L = frames[hi], frames[lo]
    vs = sorted(set(H["_variant"]) & set(L["_variant"]),
                key=lambda x: VARIANTS.index(x) if x in VARIANTS else 99)
    for v in vs:
        h = H[H["_variant"] == v].set_index("_drug")["_z"]
        l = L[L["_variant"] == v].set_index("_drug")["_z"]
        d = (h - l).dropna()
        if d.empty:
            continue
        un = d[~d.index.isin(exposed)]
        dz = float(d.get(PRIMARY_DRUG, np.nan))
        sd = float(un.std(ddof=1))
        rows.append({
            "disease": disease, "variant": v,
            "dz_abemaciclib": dz,
            "dz_palbociclib": float(d.get("palbociclib", np.nan)),
            "dz_trilaciclib": float(d.get("trilaciclib", np.nan)),
            "dz_ribociclib": float(d.get(NEG_CONTROL, np.nan)),
            "n_exposed": int(d.index.isin(exposed).sum()),
            "n_unexposed": len(un),
            "unexposed_SD": sd,
            "unexposed_max_abs": float(un.abs().max()),
            "standardized_effect": dz / sd if sd > 0 else np.nan,
            "abs_rank_of_abemaciclib": int((d.abs() >= abs(dz)).sum()),
            "empirical_P_vs_unexposed": (int((un >= dz).sum()) + 1) / (len(un) + 1),
        })
        for drug, val in d.items():
            curves.append({"disease": disease, "variant": v, "drug": drug,
                           "delta_z": float(val),
                           "exposed": bool(drug in exposed)})

res = pd.DataFrame(rows)
if res.empty:
    raise SystemExit("  no overlapping variants between arms")

print("\n" + "=" * 80)
print("Result")
print("=" * 80)
print(res[["disease", "variant", "dz_abemaciclib", "unexposed_SD",
           "unexposed_max_abs", "standardized_effect",
           "abs_rank_of_abemaciclib", "empirical_P_vs_unexposed"]]
      .round(4).to_string(index=False))

print("\n  Interpretation per row:")
for _, r in res.iterrows():
    tag = ("largest shift of all drugs"
           if r["abs_rank_of_abemaciclib"] == 1
           else f"rank {int(r['abs_rank_of_abemaciclib'])} by |delta-z|")
    print(f"    {r['disease']:<7} {r['variant']:<24} dz={r['dz_abemaciclib']:+.3f}"
          f"   {r['standardized_effect']:>7.1f} x unexposed SD   {tag}")

print("\n  Directional concordance (abemaciclib):")
for dis in res["disease"].unique():
    s = res[res["disease"] == dis]
    print(f"    {dis:<7} positive in {int((s['dz_abemaciclib'] > 0).sum())}/{len(s)} variants")

print("\n  Other RB1-exposed drugs and the negative control:")
print(res[["disease", "variant", "dz_palbociclib", "dz_trilaciclib",
           "dz_ribociclib"]].round(4).to_string(index=False))

os.makedirs(OUT_DIR, exist_ok=True)
res.to_csv(os.path.join(OUT_DIR, "amendment1_scale_free_effect.csv"),
           index=False, encoding="utf-8-sig")
pd.DataFrame(curves).to_csv(os.path.join(OUT_DIR, "amendment1_delta_z_all_drugs.csv"),
                            index=False, encoding="utf-8-sig")
with open(os.path.join(OUT_DIR, "amendment1_summary.json"), "w", encoding="utf-8") as fh:
    json.dump({
        "note": ("Post-hoc. Frozen L1/L2 used an absolute delta-z cut calibrated "
                 "on the primary variant; not scale-comparable across variants."),
        "RB1_exposed_in_universe": sorted(exposed),
        "RB1_exposed_outside_universe": sorted(outside),
        "abemaciclib_largest_shift_in_all_variants": bool(
            (res["abs_rank_of_abemaciclib"] == 1).all()),
        "direction_concordant": bool((res["dz_abemaciclib"] > 0).all()),
        "standardized_effect_range": [float(res["standardized_effect"].min()),
                                      float(res["standardized_effect"].max())],
    }, fh, indent=2, ensure_ascii=False)

print(f"\n  outputs -> {OUT_DIR}\n")
