"""Stage-B network audit for protocol melanoma_R1_v0.3.

Propagation, C_MEL scoring, degree-matched nulls, audit metrics (PROJECT_PLAN
section 4.1), H4/H5/H6 (sections 4.2, 5.2), M0 leave-one-seed-out, S1 (108
universe), S10 (alpha / STRING threshold / k) and the stage-A seed branches.

Place this file and melanoma_rules.py in the analysis/ package of the ACC
repository; it reuses the frozen ACC/breast implementations:
  analysis.method_strengthening   load_string_graph, build_restart_matrix,
                                  random_walk_with_restart, weighted_pagerank,
                                  build_association_matrix, compute_c_acc_matrix,
                                  percentile_columns, generate_degree_matched_seed_sets,
                                  partial_spearman
  analysis.normalization_sensitivity  _gene_variants, network_smooth_with_restart,
                                  symmetric_normalized_operator

Differences from the breast implementation, all prescribed by the protocol:
  * 124-drug primary universe (108 = S1); percentiles within each universe.
  * Five RNG streams x 2000 draws, default_rng(SeedSequence([stream, arm])).
  * Null percentiles kept in float64 (breast stored float32).
  * Top-k ties broken by drug name.

Usage (from the repository root):
  python -m analysis.melanoma_audit --project <revision_dir> --string-dir <dir> --out <dir> [--benchmark]
  python -m analysis.melanoma_audit --validate-breast --string-dir <dir> --out <dir>
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import scipy
from scipy import sparse

from analysis.method_strengthening import (
    build_association_matrix,
    build_restart_matrix,
    compute_c_acc_matrix,
    generate_degree_matched_seed_sets,
    load_string_graph,
    partial_spearman,
    percentile_columns,
    random_walk_with_restart,
    weighted_pagerank,
)
from analysis.normalization_sensitivity import (
    _gene_variants,
    network_smooth_with_restart,
    symmetric_normalized_operator,
)
from analysis import melanoma_rules as rules

ANALYSIS_VERSION = "melanoma-audit-stageB-v1"
VARIANTS = rules.VARIANTS
PRIMARY = rules.PRIMARY
REPO_ROOT = Path(__file__).resolve().parents[1]
BRAF_I = ("Dabrafenib", "Vemurafenib", "Encorafenib")
MEK_I = ("Trametinib", "Cobimetinib", "Binimetinib")


# ------------------------------------------------------------------ I/O
def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f"{v:.12g}" if isinstance(v, float) else v) for k, v in row.items()})


def jsonable(obj):
    if isinstance(obj, Fraction):
        return f"{obj.numerator}/{obj.denominator}"
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        return None if math.isnan(float(obj)) else float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def frozen_hashes(path: Path) -> dict[str, str]:
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split("  ")
        if len(parts) == 2 and len(parts[0]) == 64:
            out[parts[1].strip()] = parts[0]
    return out


def check_freezes(project: Path) -> dict[str, dict[str, str]]:
    """Verify stage-A, A1 addendum and stage-B freeze records before running."""
    protocol = project / "04_melanoma_protocol"
    checked = {}
    for record, base in (("FREEZE_STAGE_A.txt", protocol), ("FREEZE_STAGE_A_ADDENDUM.txt", protocol),
                         ("FREEZE_STAGE_B.txt", project)):
        listed = frozen_hashes(protocol / record)
        if not listed:
            raise RuntimeError(f"{record} lists no files")
        for name, digest in listed.items():
            if sha256(base / name) != digest:
                raise RuntimeError(f"file changed since freeze ({record}): {name}")
        checked[record] = listed
    config = json.loads((protocol / "protocol_config.json").read_text(encoding="utf-8"))
    if config.get("status") != "STAGE_A_FROZEN":
        raise RuntimeError("protocol config is not the frozen stage-A config")
    return checked


# ---------------------------------------------------------------- engine
class Engine:
    """STRING graph + 399-gene projection + 124-drug association matrix."""

    def __init__(self, info: Path, links: Path, assoc: Mapping[str, set[str]],
                 seed_genes: set[str], threshold: int, alpha: float):
        self.drugs = tuple(sorted(assoc))
        self.genes = tuple(sorted(set().union(*assoc.values())))
        self.graph = load_string_graph(info, links, required_nodes=set(self.genes) | set(seed_genes),
                                       threshold=threshold)
        self.threshold = threshold
        self.symmetric = symmetric_normalized_operator(self.graph.adjacency)
        self.gene_index = np.asarray([self.graph.node_index[g] for g in self.genes], dtype=int)
        self.A = build_association_matrix(self.drugs, self.genes, assoc)
        self.edge_counts = np.asarray(self.A.sum(axis=0)).ravel()
        self.n_assoc = np.asarray(self.A.sum(axis=1)).ravel()
        self.drug_index = {d: i for i, d in enumerate(self.drugs)}
        self._uniform: dict[float, np.ndarray] = {}
        self.alpha = alpha
        self.iterations: list[tuple[int, int]] = []

    def uniform(self, alpha: float) -> np.ndarray:
        if alpha not in self._uniform:
            n = len(self.graph.node_names)
            col, _, _ = random_walk_with_restart(self.graph, np.full((n, 1), 1.0 / n), alpha=alpha)
            self._uniform[alpha] = col[self.gene_index, :]
        return self._uniform[alpha]

    def gene_variants(self, weight_columns: Sequence[Mapping[str, float]], alpha: float | None = None):
        alpha = self.alpha if alpha is None else alpha
        restart = build_restart_matrix(self.graph, weight_columns)
        column, it_c, _ = random_walk_with_restart(self.graph, restart, alpha=alpha)
        symmetric, it_s, _ = network_smooth_with_restart(self.symmetric, restart, alpha=alpha)
        self.iterations.append((it_c, it_s))
        return _gene_variants(column[self.gene_index, :], symmetric[self.gene_index, :], self.uniform(alpha))

    def c_raw(self, gene_values: np.ndarray, k: float = 3.0) -> np.ndarray:
        return compute_c_acc_matrix(self.A, gene_values, pseudo_count=k)

    def background(self, gene_values: np.ndarray) -> np.ndarray:
        return (self.edge_counts @ gene_values) / self.edge_counts.sum()

    def c_subset(self, gene_values: np.ndarray, keep: sparse.csr_matrix, k: float = 3.0) -> np.ndarray:
        """C over a subset of each drug's genes, keeping the full-edge background mu."""
        mu = self.background(gene_values)
        n_sub = np.asarray(keep.sum(axis=1)).ravel()
        return (keep @ gene_values + k * mu[None, :]) / (n_sub[:, None] + k)


def pct(raw: np.ndarray, rows: np.ndarray) -> np.ndarray:
    return percentile_columns(raw[rows, :])


# ------------------------------------------------------------------ nulls
def null_seed_sets(engine: Engine, weights: Mapping[str, float], arm_index: int,
                   streams: Sequence[int], per_stream: int, bins: int) -> tuple[list[dict], list]:
    rows: list[dict] = []
    edges_used = []
    for s, base in enumerate(streams):
        draws, edges = generate_degree_matched_seed_sets(
            weights, engine.graph.node_names, engine.graph.degree, n_draws=per_stream,
            rng_seed=np.random.SeedSequence([int(base), int(arm_index)]), n_bins=bins)
        for row in draws:
            row["stream"] = s
            row["replicate"] = int(row["replicate"]) + s * per_stream
        rows.extend(draws)
        edges_used.append([float(e) for e in edges])
    return rows, edges_used


def null_raw_scores(engine: Engine, rows: Sequence[Mapping], n_null: int, batch: int, k: float = 3.0):
    by_rep: dict[int, dict[str, float]] = {}
    for row in rows:
        by_rep.setdefault(int(row["replicate"]), {})[row["null_seed"]] = float(row["weight"])
    out = {v: np.empty((len(engine.drugs), n_null), dtype=np.float64) for v in VARIANTS}
    for start in range(0, n_null, batch):
        reps = list(range(start, min(start + batch, n_null)))
        variants = engine.gene_variants([by_rep[r] for r in reps])
        for v in VARIANTS:
            out[v][:, reps] = engine.c_raw(variants[v], k)
    return out


# --------------------------------------------------------------- per arm
def arm_tables(engine: Engine, name: str, weights: Mapping[str, float], obs_raw: dict, null_raw: dict,
               universes: Mapping[str, np.ndarray], six: Sequence[str], n_streams: int = 5) -> dict:
    """Per-variant, per-universe drug statistics and H4 group statistics."""
    result: dict[str, Any] = {"drug_rows": [], "group": {}, "z": {}, "pct": {}}
    for uname, rows_idx in universes.items():
        names = [engine.drugs[i] for i in rows_idx]
        result["z"][uname], result["pct"][uname], result["group"][uname] = {}, {}, {}
        for v in VARIANTS:
            obs = pct(obs_raw[v], rows_idx)[:, 0]
            nul = pct(null_raw[v], rows_idx)
            n_null = nul.shape[1]
            p_vals, bs, zs = [], [], []
            for i in range(len(names)):
                b, p = rules.empirical_upper(obs[i], nul[i])
                bs.append(b)
                p_vals.append(p)
                zs.append(rules.z_score(obs[i], nul[i]))
            q_vals = rules.benjamini_hochberg(p_vals)
            ranks = rules.ordinal_ranks(names, obs)
            for i, d in enumerate(names):
                lo, hi = rules.clopper_pearson(bs[i], n_null)
                result["drug_rows"].append({
                    "arm": name, "universe": uname, "variant": v, "drug": d, "rank": ranks[d],
                    "C_raw": float(obs_raw[v][rows_idx[i], 0]), "percentile": float(obs[i]),
                    "null_mean": float(nul[i].mean()), "null_sd": float(nul[i].std(ddof=1)),
                    "z": zs[i], "exceedances": bs[i], "p_upper": p_vals[i],
                    "tail_cp95_low": lo, "tail_cp95_high": hi, "q_bh_within": q_vals[i],
                    "null_draws": n_null})
            result["z"][uname][v] = dict(zip(names, zs, strict=True))
            result["pct"][uname][v] = dict(zip(names, obs, strict=True))
            six_idx = [names.index(d) for d in six]
            obs_group = float(obs[six_idx].mean())
            null_group = nul[six_idx, :].mean(axis=0)
            b, p = rules.empirical_upper(obs_group, null_group)
            per = n_null // n_streams
            stream_p = [rules.empirical_upper(obs_group, null_group[s * per:(s + 1) * per])[1]
                        for s in range(n_streams)]
            result["group"][uname][v] = {"observed": obs_group, "null_mean": float(null_group.mean()),
                                         "p": p, "exceedances": b, "p_by_stream": stream_p}
    return result


def exposure_vectors(engine: Engine, weights: Mapping[str, float], k: float = 3.0):
    seeds = set(weights)
    total = sum(weights.values())
    lc, lw, hits = [], [], []
    for d in engine.drugs:
        a = {g for g, x in zip(engine.genes, engine.A[engine.drug_index[d]].toarray().ravel()) if x}
        h = sorted(a & seeds)
        hits.append(h)
        lc.append(len(h) / (len(a) + k))
        lw.append(sum(weights[g] / total for g in h) / (len(a) + k))
    return np.asarray(lc), np.asarray(lw), hits


def audit_metrics(engine: Engine, arm: str, weights: Mapping[str, float], obs_raw: dict, obs_variants: dict,
                  z: Mapping[str, Mapping[str, Mapping[str, float]]], universes: Mapping[str, np.ndarray],
                  centrality: Mapping[str, np.ndarray], top25: Mapping[str, int]) -> tuple[list, list]:
    lc, lw, hits = exposure_vectors(engine, weights)
    seed_set = set(weights)
    seed_mask = np.asarray([g in seed_set for g in engine.genes])
    keep_seed = (engine.A @ sparse.diags(seed_mask.astype(float))).tocsr()
    keep_non = (engine.A @ sparse.diags((~seed_mask).astype(float))).tocsr()
    metric_rows, decomposition_rows = [], []
    for uname, idx in universes.items():
        names = [engine.drugs[i] for i in idx]
        exposed = np.asarray([bool(hits[i]) for i in idx])
        for v in VARIANTS:
            c_pct = pct(obs_raw[v], idx)[:, 0]
            top = set(rules.top_k(names, c_pct, top25[uname]))
            top20_c = set(rules.top_k(names, c_pct, 20))
            top20_l = set(rules.top_k(names, lc[idx], 20))
            unexposed = [n for n, e in zip(names, exposed) if not e]
            u = Fraction(len(set(unexposed) & top), len(unexposed)) if unexposed else None
            zv = np.asarray([z[uname][v][n] for n in names])
            ok = ~np.isnan(zv)
            try:
                partial = partial_spearman(c_pct[ok], zv[ok], lc[idx][ok])
            except ValueError:
                partial = math.nan
            nonseed_pct = pct(engine.c_subset(obs_variants[v], keep_non), idx)[:, 0]
            seed_pct = pct(engine.c_subset(obs_variants[v], keep_seed), idx)[:, 0]
            ranks_c = rules.ordinal_ranks(names, c_pct)
            ranks_n = rules.ordinal_ranks(names, nonseed_pct)
            row = {"arm": arm, "universe": uname, "variant": v,
                   "exposed": int(exposed.sum()), "unexposed": len(unexposed),
                   "coverage_c": f"{int(exposed.sum())}/{len(names)}",
                   "unexposed_in_top25": len(set(unexposed) & top), "u": jsonable(u) if u is not None else "NA",
                   "max_unexposed_pct": float(max(c_pct[~exposed])) if (~exposed).any() else math.nan,
                   "rho_Lcount_C": rules.spearman(lc[idx], c_pct), "rho_Lweight_C": rules.spearman(lw[idx], c_pct),
                   "rho_nassoc_C": rules.spearman(engine.n_assoc[idx], c_pct),
                   "top20_intersection_Lcount": len(top20_c & top20_l),
                   "rho_C_z": rules.spearman(c_pct[ok], zv[ok]),
                   "partial_rho_C_z_given_Lcount": partial,
                   "rho_nonseed_C": rules.spearman(nonseed_pct, c_pct),
                   "top20_intersection_nonseed": len(top20_c & set(rules.top_k(names, nonseed_pct, 20))),
                   "drugs_without_nonseed_gene": sum(int(engine.n_assoc[i] - len(hits[i])) == 0 for i in idx)}
            for cname, cvec in centrality.items():
                row[f"rho_mean_{cname}_C"] = rules.spearman(cvec[idx], c_pct)
            metric_rows.append(row)
            for i, n in enumerate(names):
                decomposition_rows.append({"arm": arm, "universe": uname, "variant": v, "drug": n,
                                           "seed_hits": ";".join(hits[idx[i]]), "L_count": float(lc[idx[i]]),
                                           "L_weight": float(lw[idx[i]]), "pct_full": float(c_pct[i]),
                                           "pct_seed_only": float(seed_pct[i]), "pct_nonseed": float(nonseed_pct[i]),
                                           "rank_full": ranks_c[n], "rank_nonseed": ranks_n[n],
                                           "rank_shift_nonseed": ranks_n[n] - ranks_c[n]})
    return metric_rows, decomposition_rows


def drug_centrality(engine: Engine) -> dict[str, np.ndarray]:
    pagerank, _, _ = weighted_pagerank(engine.graph)
    gi = engine.gene_index
    per_gene = {"degree": engine.graph.degree[gi], "strength": engine.graph.strength[gi], "pagerank": pagerank[gi]}
    return {k: np.asarray(engine.A @ v).ravel() / engine.n_assoc for k, v in per_gene.items()}


def compare_rankings(names: Sequence[str], base: np.ndarray, other: np.ndarray, six: Sequence[str]) -> dict:
    rb, ro = rules.ordinal_ranks(names, base), rules.ordinal_ranks(names, other)
    shifts = np.asarray([abs(ro[n] - rb[n]) for n in names])
    tb, to = set(rules.top_k(names, base, 20)), set(rules.top_k(names, other, 20))
    out = {"spearman_vs_base": rules.spearman(base, other), "top20_jaccard": len(tb & to) / len(tb | to),
           "max_abs_shift": int(shifts.max()), "median_abs_shift": float(np.median(shifts)),
           "n_shift_ge5": int((shifts >= 5).sum()), "n_shift_ge10": int((shifts >= 10).sum())}
    for d in six:
        out[f"rank_{d}"] = ro[d]
    return out


# ------------------------------------------------------------------ main
def load_inputs(project: Path):
    net = REPO_ROOT / "data" / "bindex_network"
    assoc: dict[str, set[str]] = {}
    for e in read_csv(net / "bindex_edges_1304.csv"):
        assoc.setdefault(e["drug"], set()).add(e["gene"])
    u108 = sorted({r["drug"] for r in read_csv(REPO_ROOT / "results" / "positive_control" /
                                                "positive_control_drug_ranks_primary108.csv")})
    stage_a = project / "results" / "melanoma_stageA_seeds"
    branches = {}
    for path in sorted(stage_a.glob("*/seeds_frozen_string_checked.csv")):
        rows = read_csv(path)
        if any(r["in_STRING_v12"] != "yes" for r in rows):
            raise RuntimeError(f"unmapped seed in {path}; protocol forbids synthetic isolated seed nodes")
        arms: dict[str, dict[str, float]] = {}
        for r in rows:
            arms.setdefault(r["arm"], {})[r["gene"]] = float(r["normalized_weight"])
        branches[path.parent.name] = arms
    if "primary_metastasis" not in branches:
        raise RuntimeError("primary seeds_frozen_string_checked.csv not found")
    return assoc, u108, branches


def run(project: Path, string_dir: Path, out: Path, benchmark: bool = False) -> int:
    t0 = time.time()
    freezes = check_freezes(project)
    config = json.loads((project / "04_melanoma_protocol" / "protocol_config.json").read_text(encoding="utf-8"))
    stageb = json.loads((project / "results" / "melanoma_stageB_prep" / "stageB_exposure_anchor_summary.json")
                        .read_text(encoding="utf-8"))
    assoc, u108, branches = load_inputs(project)
    six = list(config["six_drugs"])
    alpha, k = float(config["alpha_restart"]), float(config["pseudo_count_k"])
    all_seeds = set().union(*(set(w) for arms in branches.values() for w in arms.values()))
    info = string_dir / "9606.protein.info.v12.0.txt.gz"
    links = string_dir / "9606.protein.links.v12.0.txt.gz"
    engine = Engine(info, links, assoc, all_seeds, int(config["string_threshold"]), alpha)
    universes = {"124": np.arange(len(engine.drugs)),
                 "108": np.asarray([engine.drug_index[d] for d in u108])}
    top25 = {"124": int(config["top25_cutoff_rank"]["124"]), "108": int(config["top25_cutoff_rank"]["108"])}
    primary = branches["primary_metastasis"]
    streams = [int(s) for s in config["rng_stream_seeds"]]
    per_stream = int(config["draws_per_rng_stream"])
    n_null = int(config["n_null_per_arm"])
    assert per_stream * len(streams) == n_null
    batch = int(config["null_batch_size"])

    if benchmark:
        rows, _ = null_seed_sets(engine, primary["M0"], 99, streams[:1], 100, int(config["requested_degree_bins"]))
        tb = time.time()
        null_raw_scores(engine, rows, 100, batch, k)
        sec = time.time() - tb
        report = {"benchmark_draws": 100, "seconds": sec, "projected_seconds_2_arms_10000": sec * 200,
                  "projected_with_30pct_margin": sec * 200 * 1.3, "note": "benchmark arm_index 99; no scientific output"}
        out.mkdir(parents=True, exist_ok=True)
        (out / "benchmark.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(json.dumps(report, indent=2))
        return 0

    centrality = drug_centrality(engine)
    arms_out, null_meta = {}, {}
    for arm_index, arm in enumerate(("M0", "M1")):
        weights = primary[arm]
        obs_variants = engine.gene_variants([weights])
        obs_raw = {v: engine.c_raw(obs_variants[v], k) for v in VARIANTS}
        rows, bin_edges = null_seed_sets(engine, weights, arm_index, streams, per_stream,
                                         int(config["requested_degree_bins"]))
        write_csv(out / f"null_seed_sets_{arm}.csv", rows)
        null_raw = null_raw_scores(engine, rows, n_null, batch, k)
        np.savez_compressed(out / f"null_C_raw_{arm}.npz", **null_raw)
        np.savez_compressed(out / f"observed_C_raw_{arm}.npz", drugs=np.asarray(engine.drugs),
                            **{v: obs_raw[v][:, 0] for v in VARIANTS})
        tables = arm_tables(engine, arm, weights, obs_raw, null_raw, universes, six, len(streams))
        metrics, decomposition = audit_metrics(engine, arm, weights, obs_raw, obs_variants, tables["z"],
                                               universes, centrality, top25)
        arms_out[arm] = {"tables": tables, "metrics": metrics, "decomposition": decomposition,
                         "obs_raw": obs_raw, "obs_variants": obs_variants}
        null_meta[arm] = {"degree_bin_edges_by_stream": bin_edges, "rng": [[s, arm_index] for s in streams]}

    write_csv(out / "drug_statistics.csv", arms_out["M0"]["tables"]["drug_rows"] + arms_out["M1"]["tables"]["drug_rows"])
    write_csv(out / "audit_metrics.csv", arms_out["M0"]["metrics"] + arms_out["M1"]["metrics"])
    write_csv(out / "seed_nonseed_decomposition.csv", arms_out["M0"]["decomposition"] + arms_out["M1"]["decomposition"])

    # cross-branch BH sensitivity: both arms x four variants x 124 drugs
    m0_rows = [r for r in arms_out["M0"]["tables"]["drug_rows"] + arms_out["M1"]["tables"]["drug_rows"]
               if r["universe"] == "124"]
    pooled_q = rules.benjamini_hochberg([r["p_upper"] for r in m0_rows])
    write_csv(out / "pooled_bh_992.csv", [{"arm": r["arm"], "variant": r["variant"], "drug": r["drug"],
                                           "p_upper": r["p_upper"], "q_bh_pooled": q}
                                          for r, q in zip(m0_rows, pooled_q, strict=True)])

    # H4 (124 primary; 108 = S1)
    h4 = {}
    for uname in ("124", "108"):
        g = arms_out["M0"]["tables"]["group"][uname]
        pc = arms_out["M0"]["tables"]["pct"][uname][PRIMARY]
        names = [engine.drugs[i] for i in universes[uname]]
        top = set(rules.top_k(names, [pc[n] for n in names], top25[uname]))
        h4[uname] = rules.evaluate_h4(g, {d: d in top for d in six}, BRAF_I, MEK_I)
        h4[uname]["group_by_variant"] = g
    # H6 (124 primary; 108 reported)
    x_set = list(stageb["X"])
    h6 = {}
    for uname in ("124", "108"):
        dz = {v: {d: arms_out["M1"]["tables"]["z"][uname][v][d] - arms_out["M0"]["tables"]["z"][uname][v][d]
                  for d in arms_out["M0"]["tables"]["z"][uname][v]} for v in VARIANTS}
        h6[uname] = rules.evaluate_h6(dz, x_set, six)
        h6[uname]["scale_free_rank_by_variant"] = {v: rules.scale_free_rank(dz[v], x_set) for v in VARIANTS}
        write_csv(out / f"H6_delta_z_{uname}.csv", [{"variant": v, "drug": d, "delta_z": dz[v][d],
                                                    "in_X": d in x_set, "six_drug": d in six}
                                                   for v in VARIANTS for d in sorted(dz[v])])
    # H5 (108 primary)
    anchors = {}
    for row in stageb["anchors"]:
        num, den = (int(x) for x in row["c"].split("/"))
        un, ud = (int(x) for x in row["u"].split("/"))
        anchors[row["disease"]] = {"c": Fraction(num, den), "rho": float(row["rho_L_count"]), "u": Fraction(un, ud)}
    lc, _, hits = exposure_vectors(engine, primary["M0"])
    h5 = {}
    for uname in ("108", "124"):
        idx = universes[uname]
        names = [engine.drugs[i] for i in idx]
        c_pct = pct(arms_out["M0"]["obs_raw"][PRIMARY], idx)[:, 0]
        exposed = [bool(hits[i]) for i in idx]
        unexposed = [n for n, e in zip(names, exposed) if not e]
        top = set(rules.top_k(names, c_pct, top25[uname]))
        c_mel = Fraction(sum(exposed), len(names))
        u_mel = Fraction(len(set(unexposed) & top), len(unexposed)) if unexposed else None
        rho = rules.spearman(lc[idx], c_pct)
        bounds = rules.h5_feasible_bounds(lc[idx], top25[uname])
        if uname == "108":
            h5[uname] = rules.evaluate_h5(c_mel, rho, u_mel, anchors, bounds)
        else:
            h5[uname] = {"verdict": "PENDING_124_ANCHORS", "note": "anchors in 124 are computed in task 5"}
        h5[uname].update({"c": c_mel, "rho_Lcount_C": rho, "u": u_mel, "feasible_bounds": bounds})

    # LOSO over M0 seeds (deterministic, no nulls)
    loso_rows = []
    base_pct = {v: pct(arms_out["M0"]["obs_raw"][v], universes["124"])[:, 0] for v in VARIANTS}
    names124 = list(engine.drugs)
    for seed in sorted(primary["M0"]):
        rest = {g: w for g, w in primary["M0"].items() if g != seed}
        if not rest:
            continue
        gv = engine.gene_variants([rest])
        direct = {d for d in names124 if seed in assoc[d]}
        for v in VARIANTS:
            other = pct(engine.c_raw(gv[v], k), universes["124"])[:, 0]
            row = {"omitted_seed": seed, "variant": v, **compare_rankings(names124, base_pct[v], other, six)}
            rb, ro = rules.ordinal_ranks(names124, base_pct[v]), rules.ordinal_ranks(names124, other)
            ds = [abs(ro[d] - rb[d]) for d in direct]
            us = [abs(ro[d] - rb[d]) for d in names124 if d not in direct]
            row.update({"n_direct": len(direct), "max_shift_direct": max(ds) if ds else "NA",
                        "max_shift_unassociated": max(us)})
            loso_rows.append(row)
    write_csv(out / "loso_M0.csv", loso_rows)

    # S10: alpha / k with fixed M0 seeds; STRING threshold 700 below
    s10_rows = []
    for a in config["sensitivity_branches"]["S10"]["alpha"]:
        if float(a) == alpha:
            continue
        gv = engine.gene_variants([primary["M0"]], alpha=float(a))
        for v in VARIANTS:
            other = pct(engine.c_raw(gv[v], k), universes["124"])[:, 0]
            s10_rows.append({"factor": "alpha", "value": a, "variant": v,
                             **compare_rankings(names124, base_pct[v], other, six)})
    for kk in config["sensitivity_branches"]["S10"]["k"]:
        if float(kk) == k:
            continue
        for v in VARIANTS:
            other = pct(engine.c_raw(arms_out["M0"]["obs_variants"][v], float(kk)), universes["124"])[:, 0]
            s10_rows.append({"factor": "k", "value": kk, "variant": v,
                             **compare_rankings(names124, base_pct[v], other, six)})

    # seed branches (deterministic)
    branch_rows = []
    for bname, arms in sorted(branches.items()):
        if bname == "primary_metastasis":
            continue
        for arm, weights in sorted(arms.items()):
            gv = engine.gene_variants([weights])
            lcb, _, hb = exposure_vectors(engine, weights)
            for v in VARIANTS:
                other = pct(engine.c_raw(gv[v], k), universes["124"])[:, 0]
                base = pct(arms_out[arm]["obs_raw"][v], universes["124"])[:, 0]
                exposed = np.asarray([bool(h) for h in hb])
                top = set(rules.top_k(names124, other, top25["124"]))
                unexp = [n for n, e in zip(names124, exposed) if not e]
                branch_rows.append({"branch": bname, "arm": arm, "variant": v, "n_seeds": len(weights),
                                    "exposed": int(exposed.sum()),
                                    "unexposed_in_top25": len(set(unexp) & top),
                                    "rho_Lcount_C": rules.spearman(lcb, other),
                                    **compare_rankings(names124, base, other, six)})
    write_csv(out / "seed_branches.csv", branch_rows)

    # STRING >= 700 (new graph, same seeds and alpha)
    engine700 = Engine(info, links, assoc, all_seeds, 700, alpha)
    gv = engine700.gene_variants([primary["M0"]])
    for v in VARIANTS:
        other = pct(engine700.c_raw(gv[v], k), universes["124"])[:, 0]
        s10_rows.append({"factor": "string_threshold", "value": 700, "variant": v,
                         **compare_rankings(names124, base_pct[v], other, six)})
    write_csv(out / "S10_parameter_sensitivity.csv", s10_rows)

    decisions = {"H4": h4, "H5": h5, "H6": h6, "X": x_set,
                 "H5b": "NOT_EVALUABLE (all six drugs unexposed; STAGE_B_PREDICTIONS.md section 3)"}
    (out / "decisions.json").write_text(json.dumps(jsonable(decisions), indent=2, ensure_ascii=False), encoding="utf-8")
    its = np.asarray(engine.iterations)
    manifest = {"analysis_version": ANALYSIS_VERSION, "freezes_verified": freezes,
                "code_sha256": {p.name: sha256(p) for p in (Path(__file__), Path(rules.__file__))},
                "string_sha256": {info.name: sha256(info), links.name: sha256(links)},
                "graph": {"nodes": len(engine.graph.node_names), "edge_rows_ge_threshold": engine.graph.source_edge_rows,
                          "nodes_700": len(engine700.graph.node_names), "edge_rows_700": engine700.graph.source_edge_rows},
                "null": null_meta, "n_null": n_null, "batch": batch,
                "iterations": {"column_max": int(its[:, 0].max()), "symmetric_max": int(its[:, 1].max())},
                "python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__,
                "platform": platform.platform(), "wall_seconds": time.time() - t0}
    (out / "run_manifest.json").write_text(json.dumps(jsonable(manifest), indent=2), encoding="utf-8")
    print(json.dumps(jsonable({"H4": {u: h4[u]["status"] for u in h4}, "H5_108": h5["108"]["verdict"],
                               "H5_annotation": h5["108"].get("annotation"),
                               "H6": {u: h6[u]["verdict"] for u in h6}}), indent=2))
    return 0


# ------------------------------------------------------- breast regression
def validate_breast(string_dir: Path, out: Path, tolerance: float = 5e-11) -> int:
    """Reproduce archived breast observed percentiles (108, four variants) with this engine."""
    net = REPO_ROOT / "data" / "bindex_network"
    assoc: dict[str, set[str]] = {}
    for e in read_csv(net / "bindex_edges_1304.csv"):
        assoc.setdefault(e["drug"], set()).add(e["gene"])
    archived = read_csv(REPO_ROOT / "results" / "positive_control" / "positive_control_drug_ranks_primary108.csv")
    raw = {r["gene"]: float(r["raw_weight"]) for r in
           read_csv(REPO_ROOT / "results" / "positive_control" / "positive_control_seed_frozen.csv")
           if r["include_primary"] == "yes"}
    engine = Engine(string_dir / "9606.protein.info.v12.0.txt.gz", string_dir / "9606.protein.links.v12.0.txt.gz",
                    assoc, set(raw), 400, 0.4)
    u108 = sorted({r["drug"] for r in archived})
    idx = np.asarray([engine.drug_index[d] for d in u108])
    variants = engine.gene_variants([raw])
    rows, worst = [], 0.0
    for v in VARIANTS:
        ours = dict(zip(u108, pct(engine.c_raw(variants[v]), idx)[:, 0], strict=True))
        for r in archived:
            if r["variant"] != v:
                continue
            diff = abs(ours[r["drug"]] - float(r["observed_C_ACC_percentile"]))
            worst = max(worst, diff)
            rows.append({"variant": v, "drug": r["drug"], "archived": float(r["observed_C_ACC_percentile"]),
                         "reproduced": float(ours[r["drug"]]), "abs_diff": diff})
    write_csv(out / "breast_regression.csv", rows)
    status = "PASS" if worst <= tolerance else "FAIL"
    print(json.dumps({"breast_regression": status, "max_abs_diff": worst, "rows": len(rows)}, indent=2))
    return 0 if status == "PASS" else 3


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path)
    parser.add_argument("--string-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--benchmark", action="store_true")
    parser.add_argument("--validate-breast", action="store_true")
    args = parser.parse_args(argv)
    if args.validate_breast:
        return validate_breast(args.string_dir, args.out)
    if args.project is None:
        parser.error("--project is required")
    return run(args.project, args.string_dir, args.out, args.benchmark)


if __name__ == "__main__":
    raise SystemExit(main())
