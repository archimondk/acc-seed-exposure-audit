"""Deterministic propagation of the S4 and S5 seed branches (no nulls).

Same outputs and conventions as the seed_branches.csv block of analysis.melanoma_audit:
124-drug universe, four normalizations, each branch arm compared with the primary arm of
the same name. Before scoring, the engine reproduces the stage-B primary M0/M1 percentiles
(drug_statistics.csv, 124 drugs) and stops if any differs by more than 5e-11.

Usage (from the repository root):
  python -m analysis.melanoma_branches_s4_s5 --project <revision_dir> \
      --branches <revision_dir>/results/melanoma_S4_S6 --string-dir <dir> --out <new_dir>
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
from pathlib import Path

import numpy as np
import scipy

from analysis import melanoma_audit as ma
from analysis import melanoma_rules as rules

TOL = 5e-11
BRANCHES = ("S4_site_conflict_excluded", "S5_hc3_adjusted")


def load_arms(path: Path) -> dict[str, dict[str, float]]:
    rows = ma.read_csv(path)
    if any(r["in_STRING_v12"] != "yes" for r in rows):
        raise RuntimeError(f"unmapped seed in {path}; protocol forbids synthetic isolated seed nodes")
    arms: dict[str, dict[str, float]] = {}
    for r in rows:
        arms.setdefault(r["arm"], {})[r["gene"]] = float(r["normalized_weight"])
    for arm, w in arms.items():
        if len(w) < 2:
            raise RuntimeError(f"{path}: arm {arm} has fewer than 2 mapped seeds")
    return arms


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", type=Path, required=True)
    ap.add_argument("--branches", type=Path, required=True)
    ap.add_argument("--string-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    if args.out.exists() and any(args.out.iterdir()):
        raise RuntimeError(f"refusing to overwrite non-empty {args.out}")
    freezes = ma.check_freezes(args.project)
    config = json.loads((args.project / "04_melanoma_protocol" / "protocol_config.json").read_text(encoding="utf-8"))
    assoc, _, stage_a = ma.load_inputs(args.project)
    primary = stage_a["primary_metastasis"]
    branches = {b: load_arms(args.branches / b / "seeds_frozen_string_checked.csv") for b in BRANCHES}
    six = list(config["six_drugs"])
    alpha, k = float(config["alpha_restart"]), float(config["pseudo_count_k"])
    cut = int(config["top25_cutoff_rank"]["124"])
    seeds = set(primary["M0"]) | set(primary["M1"]) | {g for arms in branches.values() for w in arms.values() for g in w}
    info = args.string_dir / "9606.protein.info.v12.0.txt.gz"
    links = args.string_dir / "9606.protein.links.v12.0.txt.gz"
    engine = ma.Engine(info, links, assoc, seeds, int(config["string_threshold"]), alpha)
    idx = np.arange(len(engine.drugs))
    names = list(engine.drugs)

    # regression against the stage-B run
    archived = {(r["arm"], r["variant"], r["drug"]): float(r["percentile"])
                for r in ma.read_csv(args.project / "results" / "melanoma_stageB" / "run" / "drug_statistics.csv")
                if r["universe"] == "124"}
    base, worst = {}, 0.0
    for arm in ("M0", "M1"):
        gv = engine.gene_variants([primary[arm]])
        base[arm] = {v: ma.pct(engine.c_raw(gv[v], k), idx)[:, 0] for v in ma.VARIANTS}
        for v in ma.VARIANTS:
            for n, x in zip(names, base[arm][v]):
                worst = max(worst, abs(x - archived[(arm, v, n)]))
    if worst > TOL:
        raise RuntimeError(f"stage-B primary percentiles not reproduced (max diff {worst})")

    rows = []
    for bname, arms in branches.items():
        for arm, weights in sorted(arms.items()):
            gv = engine.gene_variants([weights])
            lc, _, hits = ma.exposure_vectors(engine, weights)
            exposed = np.asarray([bool(h) for h in hits])
            unexposed = [n for n, e in zip(names, exposed) if not e]
            for v in ma.VARIANTS:
                other = ma.pct(engine.c_raw(gv[v], k), idx)[:, 0]
                top = set(rules.top_k(names, other, cut))
                rows.append({"branch": bname, "arm": arm, "variant": v, "n_seeds": len(weights),
                             "seeds": ";".join(f"{g}:{weights[g]:.6f}" for g in sorted(weights)),
                             "exposed": int(exposed.sum()), "unexposed_in_top25": len(set(unexposed) & top),
                             "arith_min": max(cut - int(exposed.sum()), 0),
                             "rho_Lcount_C": rules.spearman(lc, other),
                             **ma.compare_rankings(names, base[arm][v], other, six)})
    args.out.mkdir(parents=True, exist_ok=True)
    ma.write_csv(args.out / "seed_branches_S4_S5.csv", rows)
    manifest = {"freezes_verified": freezes, "regression_max_abs_diff": worst, "tolerance": TOL,
                "branch_seed_files_sha256": {b: ma.sha256(args.branches / b / "seeds_frozen_string_checked.csv") for b in BRANCHES},
                "code_sha256": {p.name: ma.sha256(p) for p in (Path(__file__), Path(ma.__file__), Path(rules.__file__))},
                "string_sha256": {info.name: ma.sha256(info), links.name: ma.sha256(links)},
                "python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__,
                "platform": platform.platform()}
    (args.out / "run_manifest.json").write_text(json.dumps(ma.jsonable(manifest), indent=2), encoding="utf-8")
    print(json.dumps({"regression_max_abs_diff": worst, "rows": len(rows)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
