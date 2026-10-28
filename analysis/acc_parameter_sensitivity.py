"""R1-5 / R1-6: ACC single-factor sensitivity on the current primary pipeline.

Fixed: the 45 disease-only ACC seeds and weights (seed_set_used.json), the
1304-edge association network, four normalizations. Varied one at a time:
  restart probability alpha in {0.2, 0.3, 0.5, 0.7}   (STRING >= 400, k = 3)
  STRING combined-score threshold in {150, 700, 900}  (alpha = 0.4, k = 3)
  shrinkage pseudo-count k in {1, 5, 10}              (alpha = 0.4, STRING >= 400)
Baseline alpha = 0.4, STRING >= 400, k = 3 must reproduce the archived
column_minmax C_ACC percentiles (108 drugs) before anything is reported.
Deterministic; no null models. Result-known post hoc sensitivity (revision).

Usage (repository root):
  python -m analysis.acc_parameter_sensitivity --string-dir <dir> --out <dir>
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from analysis import melanoma_rules as rules
from analysis.melanoma_audit import REPO_ROOT, Engine, compare_rankings, pct, read_csv, write_csv

VARIANTS = rules.VARIANTS
TOL = 5e-11  # archived values carry 12 significant digits
FOCAL = ("Abemaciclib", "Palbociclib", "Ribociclib", "Pralatrexate", "Pemetrexed")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--string-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    assoc: dict[str, set[str]] = {}
    for e in read_csv(REPO_ROOT / "data/bindex_network/bindex_edges_1304.csv"):
        assoc.setdefault(e["drug"], set()).add(e["gene"])
    seeds = json.loads((REPO_ROOT / "leakage_audit/observational/seed_set_used.json").read_text())["weights"]
    lower = {d.lower(): d for d in assoc}
    archived = {lower[r["drug"].lower()]: float(r["C_ACC_pct"])
                for r in read_csv(REPO_ROOT / "results/primary_analysis/ADRS_comp_primary_108.csv")}
    info = args.string_dir / "9606.protein.info.v12.0.txt.gz"
    links = args.string_dir / "9606.protein.links.v12.0.txt.gz"

    base_engine = Engine(info, links, assoc, set(seeds), 400, 0.4)
    u108 = np.asarray([base_engine.drug_index[d] for d in sorted(archived)])
    names108 = sorted(archived)
    universes = {"108": (u108, names108), "124": (np.arange(len(base_engine.drugs)), list(base_engine.drugs))}
    base_gv = base_engine.gene_variants([seeds])
    base = {(u, v): pct(base_engine.c_raw(base_gv[v]), idx)[:, 0] for u, (idx, _) in universes.items() for v in VARIANTS}
    check = max(abs(a - archived[n]) for a, n in zip(base[("108", "column_minmax")], names108))
    if check > TOL:
        raise SystemExit(f"baseline does not reproduce archived C_ACC percentiles (max diff {check})")

    seed_set = set(seeds)
    lcount = {d: len(assoc[d] & seed_set) / (len(assoc[d]) + 3) for d in assoc}

    def row(factor, value, engine, gv, k):
        out = []
        for u, (idx, names) in universes.items():
            top_n = 27 if u == "108" else 31
            for v in VARIANTS:
                other = pct(engine.c_raw(gv[v], k), idx)[:, 0]
                r = {"factor": factor, "value": value, "universe": u, "variant": v,
                     **compare_rankings(names, base[(u, v)], other, [d for d in FOCAL if d in names])}
                top = set(rules.top_k(names, other, top_n))
                exposed = {d for d in names if assoc[d] & seed_set}
                r["unexposed_in_top25"] = len(top - exposed)
                r["arith_min"] = max(top_n - len(exposed), 0)
                r["rho_Lcount_C"] = rules.spearman([lcount[d] for d in names], other)
                out.append(r)
        return out

    rows = row("baseline", "alpha0.4_t400_k3", base_engine, base_gv, 3.0)
    for a in (0.2, 0.3, 0.5, 0.7):
        rows += row("alpha", a, base_engine, base_engine.gene_variants([seeds], alpha=a), 3.0)
    for kk in (1.0, 5.0, 10.0):
        rows += row("k", kk, base_engine, base_gv, kk)
    graph_sizes = {400: (len(base_engine.graph.node_names), base_engine.graph.source_edge_rows)}
    for t in (150, 700, 900):
        eng = Engine(info, links, assoc, set(seeds), t, 0.4)
        graph_sizes[t] = (len(eng.graph.node_names), eng.graph.source_edge_rows)
        rows += row("string_threshold", t, eng, eng.gene_variants([seeds]), 3.0)
    write_csv(args.out / "acc_parameter_sensitivity.csv", rows)
    (args.out / "acc_parameter_sensitivity_manifest.json").write_text(json.dumps(
        {"baseline_max_abs_diff_vs_archived": check, "graph_nodes_edge_rows": graph_sizes,
         "seeds": len(seeds), "note": "VEGFA is not a STRING v12 preferred name and enters as an isolated node"},
        indent=2), encoding="utf-8")
    print(json.dumps({"rows": len(rows), "baseline_check": check, "graphs": graph_sizes}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
