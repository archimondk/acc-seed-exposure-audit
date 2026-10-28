"""Post-freeze ACC/breast 124-drug anchors for the melanoma H5 audit.

This supplement preserves the Stage-B outputs. It verifies the archived
108-drug percentiles before evaluating the frozen H5 rule on 124 drugs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy

from analysis import melanoma_rules as rules
from analysis.melanoma_audit import Engine, check_freezes, exposure_vectors, pct, read_csv
from analysis.method_strengthening import load_disease_seed_weights

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "manuscript_revision" / "plos_one_R1_20261018" / "stageB_runtime_20261019"
STAGE_B_RUN = ROOT / "results" / "melanoma_stageB_20261019" / "run"
TOLERANCE = 5e-11


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def csv_write(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write an empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: repr(value) if isinstance(value, float) else value
                             for key, value in row.items()})


def fraction(text: str) -> Fraction:
    return Fraction(text)


def score_disease(disease: str, engine: Engine, weights: dict[str, float],
                  u108: list[str], cutoffs: dict[str, int]) -> tuple[list[dict], list[dict], dict]:
    gene_variants = engine.gene_variants([weights])
    lc, lw, hits = exposure_vectors(engine, weights)
    raw_scores = {v: engine.c_raw(gene_variants[v], 3.0) for v in rules.VARIANTS}
    universes = {"108": u108, "124": list(engine.drugs)}
    score_rows, comparison_rows, exact = [], [], {}
    archived_names = set(u108)

    for universe, drug_names in universes.items():
        idx = np.asarray([engine.drug_index[name] for name in drug_names], dtype=int)
        exposed = np.asarray([bool(hits[i]) for i in idx])
        n_unexposed = int((~exposed).sum())
        if n_unexposed == 0:
            raise RuntimeError(f"{disease} {universe}: H5 unexposed denominator is zero")
        exact[universe] = {}
        for variant in rules.VARIANTS:
            percentiles = pct(raw_scores[variant], idx)[:, 0]
            ranks = rules.ordinal_ranks(drug_names, percentiles)
            top = set(rules.top_k(drug_names, percentiles, cutoffs[universe]))
            unexposed_names = {name for name, is_exposed in zip(drug_names, exposed) if not is_exposed}
            top_unexposed = len(top & unexposed_names)
            c = Fraction(int(exposed.sum()), len(drug_names))
            u = Fraction(top_unexposed, n_unexposed)
            rho_count = rules.spearman(lc[idx], percentiles)
            rho_weight = rules.spearman(lw[idx], percentiles)
            max_unexposed = float(np.max(percentiles[~exposed]))
            exact[universe][variant] = {"c": c, "u": u, "rho": rho_count,
                                         "rho_weight": rho_weight}
            comparison_rows.append({
                "disease": disease, "universe": universe, "variant": variant,
                "n_drugs": len(drug_names), "exposed": int(exposed.sum()),
                "unexposed": n_unexposed, "coverage_c": str(c),
                "unexposed_in_top25pct": top_unexposed, "u": str(u),
                "max_unexposed_pct": max_unexposed,
                "rho_Lcount_C": rho_count, "rho_Lweight_C": rho_weight,
                "graph_nodes": len(engine.graph.node_names),
            })
            for row_index, drug in enumerate(drug_names):
                graph_index = int(idx[row_index])
                score_rows.append({
                    "disease": disease, "universe": universe,
                    "variant": variant, "drug": drug,
                    "in_archived_108": drug in archived_names, "rank": ranks[drug],
                    "C_raw": float(raw_scores[variant][graph_index, 0]),
                    "percentile": float(percentiles[row_index]),
                    "seed_hits": ";".join(hits[graph_index]),
                    "L_count": float(lc[graph_index]),
                    "L_weight": float(lw[graph_index]),
                })
    return score_rows, comparison_rows, exact


def validate_archived_108(acc_scores: list[dict], breast_scores: list[dict],
                          acc_exact: dict, breast_exact: dict, project: Path) -> dict:
    archived_acc = read_csv(ROOT / "results" / "primary_analysis" / "ADRS_comp_primary_108.csv")
    archived_breast = read_csv(ROOT / "results" / "positive_control" /
                               "positive_control_drug_ranks_primary108.csv")
    breast_108 = {(r["variant"], r["drug"]): r for r in breast_scores
                  if r["universe"] == "108"}
    acc_108 = {r["drug"]: r for r in acc_scores
               if r["universe"] == "108" and r["variant"] == rules.PRIMARY}
    if len(archived_acc) != 108 or len(archived_breast) != 432:
        raise RuntimeError("unexpected archived 108-drug row counts")
    acc_diffs = [abs(float(acc_108[r["drug"]]["percentile"]) - float(r["C_ACC_pct"]))
                 for r in archived_acc]
    breast_diffs = [abs(float(breast_108[(r["variant"], r["drug"])]["percentile"]) -
                        float(r["observed_C_ACC_percentile"])) for r in archived_breast]
    acc_max, breast_max = max(acc_diffs), max(breast_diffs)
    if acc_max > TOLERANCE or breast_max > TOLERANCE:
        raise RuntimeError(f"108-drug regression failed: ACC={acc_max}, Breast={breast_max}")

    saved = json.loads((project / "results" / "melanoma_stageB_prep" /
                        "stageB_exposure_anchor_summary.json").read_text(encoding="utf-8"))
    for row in saved["anchors"]:
        exact = acc_exact if row["disease"] == "ACC" else breast_exact
        metric = exact["108"][rules.PRIMARY]
        if metric["c"] != fraction(row["c"]) or metric["u"] != fraction(row["u"]):
            raise RuntimeError(f"{row['disease']} 108-drug count anchor mismatch")
        if abs(metric["rho"] - float(row["rho_L_count"])) > TOLERANCE:
            raise RuntimeError(f"{row['disease']} 108-drug rho anchor mismatch")
    return {"status": "PASS", "tolerance": TOLERANCE,
            "ACC_rows": 108, "ACC_max_abs_diff": acc_max,
            "Breast_rows": 432, "Breast_max_abs_diff": breast_max,
            "stageB_108_anchor_match": True}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    if out.exists() and any(out.iterdir()):
        raise RuntimeError(f"output directory is not empty: {out}")

    frozen = check_freezes(PROJECT)
    input_manifest = json.loads((ROOT / "reproducibility" / "input_manifest.json")
                                .read_text(encoding="utf-8"))["sha256"]
    required = ["9606.protein.info.v12.0.txt.gz", "9606.protein.links.v12.0.txt.gz",
                "data/bindex_network/bindex_edges_1304.csv",
                "data/ACC_P0.5C_gene_weights_v1.csv"]
    for name in required:
        if sha256(ROOT / name) != input_manifest[name]:
            raise RuntimeError(f"frozen input hash mismatch: {name}")

    association: dict[str, set[str]] = {}
    for row in read_csv(ROOT / required[2]):
        association.setdefault(row["drug"], set()).add(row["gene"])
    archived_breast = read_csv(ROOT / "results" / "positive_control" /
                               "positive_control_drug_ranks_primary108.csv")
    u108 = sorted({row["drug"] for row in archived_breast})
    if len(association) != 124 or len(u108) != 108 or not set(u108) <= set(association):
        raise RuntimeError("expected 124 B-index drugs and archived 108-drug subset")
    excluded16 = sorted(set(association) - set(u108))
    if len(excluded16) != 16:
        raise RuntimeError("expected exactly 16 previously excluded drugs")

    breast_weights = {row["gene"]: float(row["raw_weight"]) for row in
                      read_csv(ROOT / "results" / "positive_control" /
                               "positive_control_seed_frozen.csv")
                      if row["include_primary"] == "yes"}
    acc_weights = load_disease_seed_weights(ROOT / required[3])
    string_info, string_links = ROOT / required[0], ROOT / required[1]
    cutoffs = {"108": 27, "124": 31}
    all_scores, comparisons, exact = {}, [], {}
    for disease, weights, expected_nodes in (("ACC", acc_weights, 19701),
                                             ("Breast", breast_weights, 19700)):
        engine = Engine(string_info, string_links, association, set(weights), 400, 0.4)
        if len(engine.graph.node_names) != expected_nodes:
            raise RuntimeError(f"{disease} graph node count changed")
        scores, comp, disease_exact = score_disease(disease, engine, weights, u108, cutoffs)
        all_scores[disease] = scores
        comparisons.extend(comp)
        exact[disease] = disease_exact

    # The 108-drug checks include every archived breast normalization and ACC
    # primary normalization. No 124-drug H5 decision is made before they pass.
    regression = validate_archived_108(
        all_scores["ACC"], all_scores["Breast"], exact["ACC"], exact["Breast"], PROJECT)

    black_manifest_path = STAGE_B_RUN / "run_manifest.json"
    black_manifest = json.loads(black_manifest_path.read_text(encoding="utf-8"))
    for name, digest in black_manifest["string_sha256"].items():
        if sha256(ROOT / name) != digest:
            raise RuntimeError(f"Stage-B STRING hash mismatch: {name}")
    for name, digest in black_manifest["code_sha256"].items():
        if sha256(ROOT / "analysis" / name) != digest:
            raise RuntimeError(f"Stage-B analysis code changed: {name}")

    black_metric_rows = [r for r in read_csv(STAGE_B_RUN / "audit_metrics.csv")
                         if r["arm"] == "M0"]
    if len(black_metric_rows) != 8:
        raise RuntimeError("expected eight M0 melanoma metric rows")
    for row in black_metric_rows:
        comparisons.append({
            "disease": "Melanoma", "universe": row["universe"],
            "variant": row["variant"], "n_drugs": int(row["universe"]),
            "exposed": int(row["exposed"]), "unexposed": int(row["unexposed"]),
            "coverage_c": row["coverage_c"],
            "unexposed_in_top25pct": int(row["unexposed_in_top25"]),
            "u": row["u"], "max_unexposed_pct": float(row["max_unexposed_pct"]),
            "rho_Lcount_C": float(row["rho_Lcount_C"]),
            "rho_Lweight_C": float(row["rho_Lweight_C"]),
            "graph_nodes": int(black_manifest["graph"]["nodes"]),
        })
    comparisons.sort(key=lambda r: (int(r["universe"]),
                                    ("ACC", "Breast", "Melanoma").index(r["disease"]),
                                    rules.VARIANTS.index(r["variant"])))

    black_decisions_path = STAGE_B_RUN / "decisions.json"
    black_decisions = json.loads(black_decisions_path.read_text(encoding="utf-8"))
    black124 = black_decisions["H5"]["124"]
    if black124["verdict"] != "PENDING_124_ANCHORS":
        raise RuntimeError("Stage-B original H5 124 status unexpectedly changed")
    bounds = black124["feasible_bounds"]
    parsed_bounds = {"rho_min": float(bounds["rho_min"]),
                     "rho_max": float(bounds["rho_max"]),
                     "u_min": fraction(bounds["u_min"]),
                     "u_max": fraction(bounds["u_max"])}
    anchors = {disease: {"c": exact[disease]["124"][rules.PRIMARY]["c"],
                         "rho": exact[disease]["124"][rules.PRIMARY]["rho"],
                         "u": exact[disease]["124"][rules.PRIMARY]["u"]}
               for disease in ("ACC", "Breast")}
    h5_124 = rules.evaluate_h5(fraction(black124["c"]),
                               float(black124["rho_Lcount_C"]),
                               fraction(black124["u"]), anchors, parsed_bounds)
    if h5_124["verdict"] == "NOT_EVALUABLE":
        raise RuntimeError("124-drug H5 unexpectedly not evaluable")

    out.mkdir(parents=True, exist_ok=True)
    csv_write(out / "ACC_drug_scores_124.csv",
              [r for r in all_scores["ACC"] if r["universe"] == "124"])
    csv_write(out / "Breast_drug_scores_124.csv",
              [r for r in all_scores["Breast"] if r["universe"] == "124"])
    csv_write(out / "cross_disease_comparison_108_124.csv", comparisons)
    excluded_rows = [r for r in all_scores["Breast"]
                     if r["universe"] == "124" and not r["in_archived_108"]]
    if len(excluded_rows) != 16 * len(rules.VARIANTS):
        raise RuntimeError("expected four scores for each of 16 added breast drugs")
    csv_write(out / "Breast_new_16_drug_scores.csv", excluded_rows)

    def serializable(obj):
        if isinstance(obj, Fraction):
            return str(obj)
        if isinstance(obj, dict):
            return {key: serializable(value) for key, value in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [serializable(value) for value in obj]
        return obj

    decision = {
        "H5_124": serializable(h5_124),
        "primary_normalization": rules.PRIMARY,
        "anchors_124": serializable(anchors),
        "melanoma_124": {"c": black124["c"],
                         "rho_Lcount_C": black124["rho_Lcount_C"],
                         "u": black124["u"],
                         "feasible_bounds": black124["feasible_bounds"]},
        "original_H5_108": black_decisions["H5"]["108"],
        "interpretation_note": "The arithmetic-determination annotation must accompany any SUPPORTED verdict.",
    }
    (out / "H5_124_supplement.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    output_names = ("ACC_drug_scores_124.csv", "Breast_drug_scores_124.csv",
                    "Breast_new_16_drug_scores.csv", "cross_disease_comparison_108_124.csv",
                    "H5_124_supplement.json")
    manifest = {
        "status": "PASS", "analysis": "post-Stage-B cross-disease 124-drug H5 anchor supplement",
        "frozen_records_verified": {name: len(records) for name, records in frozen.items()},
        "archived_108_regression": regression,
        "n_drugs_108": 108, "n_drugs_124": 124,
        "added_drugs": excluded16,
        "graph_nodes": {"ACC": 19701, "Breast": 19700,
                        "Melanoma": int(black_manifest["graph"]["nodes"])},
        "cutoffs": cutoffs, "variants": list(rules.VARIANTS),
        "python": sys.version.split()[0], "numpy": np.__version__,
        "scipy": scipy.__version__, "platform": platform.platform(),
        "input_sha256": {name: sha256(ROOT / name) for name in required},
        "other_input_sha256": {
            "positive_control_seed_frozen.csv": sha256(ROOT / "results" / "positive_control" /
                                                       "positive_control_seed_frozen.csv"),
            "positive_control_drug_ranks_primary108.csv": sha256(ROOT / "results" /
                                                                  "positive_control" /
                                                                  "positive_control_drug_ranks_primary108.csv"),
            "ADRS_comp_primary_108.csv": sha256(ROOT / "results" / "primary_analysis" /
                                                "ADRS_comp_primary_108.csv"),
            "stageB_decisions.json": sha256(black_decisions_path),
            "stageB_audit_metrics.csv": sha256(STAGE_B_RUN / "audit_metrics.csv"),
        },
        "code_sha256": {name: sha256(ROOT / "analysis" / name)
                        for name in ("melanoma_cross_disease_124.py", "melanoma_audit.py",
                                     "melanoma_rules.py", "method_strengthening.py")},
        "output_sha256": {name: sha256(out / name) for name in output_names},
    }
    (out / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "archived_108": regression,
                      "n_added_breast_drugs": len(excluded16),
                      "anchors_124": serializable(anchors),
                      "H5_124": serializable(h5_124)}, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
