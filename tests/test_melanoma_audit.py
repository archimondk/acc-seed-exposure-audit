"""End-to-end and unit tests for analysis.melanoma_audit on a synthetic network.

The synthetic STRING files, association network, seeds, configuration and
freeze records are generated in a temporary directory, so the tests never
touch real data.
"""
import csv
import gzip
import hashlib
import json
import math
import random
from pathlib import Path

import numpy as np
import pytest

from analysis import melanoma_audit as ma

SIX = ["Dabrafenib", "Vemurafenib", "Encorafenib", "Trametinib", "Cobimetinib", "Binimetinib"]


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write_rows(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def build_world(root: Path, seed: int = 7):
    rng = random.Random(seed)
    repo, project, string_dir = root / "repo", root / "project", root / "string"
    assoc_genes = [f"A{i}" for i in range(28)] + ["ARID2", "PTEN", "BRAF"]
    nodes = assoc_genes + ["TP53", "CDKN2A"] + [f"N{i}" for i in range(50)]
    string_dir.mkdir(parents=True)
    with gzip.open(string_dir / "9606.protein.info.v12.0.txt.gz", "wt") as h:
        h.write("#string_protein_id\tpreferred_name\tprotein_size\tannotation\n")
        for i, g in enumerate(nodes):
            h.write(f"9606.P{i}\t{g}\t100\tx\n")
    with gzip.open(string_dir / "9606.protein.links.v12.0.txt.gz", "wt") as h:
        h.write("protein1 protein2 combined_score\n")
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                if rng.random() < 0.08 or j == i + 1:
                    s = rng.randint(150, 999)
                    h.write(f"9606.P{i} 9606.P{j} {s}\n9606.P{j} 9606.P{i} {s}\n")
    drugs = SIX + ["Dacarbazine", "Tioguanine"] + [f"D{i:02d}" for i in range(32)]
    edges = []
    for d in drugs:
        pool = [g for g in assoc_genes if g != "BRAF"]
        genes = set(rng.sample(pool, rng.randint(1, 6)))
        if d in ("Dacarbazine", "Tioguanine"):
            genes.add("BRAF")
        if d in SIX:
            genes -= {"ARID2", "PTEN"}
            genes = genes or {"A0"}
        edges += [{"gene": g, "NSC": "0", "drug": d} for g in sorted(genes)]
    net = repo / "data" / "bindex_network"
    write_rows(net / "bindex_edges_1304.csv", edges)
    u108 = drugs[:35]
    write_rows(repo / "results" / "positive_control" / "positive_control_drug_ranks_primary108.csv",
               [{"variant": "column_minmax", "drug": d} for d in u108])

    protocol = project / "04_melanoma_protocol"
    protocol.mkdir(parents=True)
    # keys read by melanoma_audit.run(); values mirror protocol_config.json except
    # the null size, cut-offs and S10 grid, which are reduced for the toy network
    config = {"status": "STAGE_A_FROZEN", "execution_allowed": "stage_A_only", "six_drugs": SIX,
              "alpha_restart": 0.4, "pseudo_count_k": 3, "string_threshold": 400,
              "rng_stream_seeds": [2026101801, 2026101802, 2026101803, 2026101804, 2026101805],
              "n_null_per_arm": 20, "draws_per_rng_stream": 4, "null_batch_size": 8,
              "top25_cutoff_rank": {"124": 10, "108": 9}, "requested_degree_bins": 3,
              "sensitivity_branches": {"S10": {"alpha": [0.3, 0.4], "string_threshold": [400, 700], "k": [1, 3]}}}
    (protocol / "protocol_config.json").write_text(json.dumps(config), encoding="utf-8")
    (protocol / "PROJECT_PLAN.md").write_text("toy plan", encoding="utf-8")
    (protocol / "prefreeze_network_check.csv").write_text("toy", encoding="utf-8")
    (protocol / "AMENDMENT_A1_braf_classification.md").write_text("toy A1", encoding="utf-8")
    (protocol / "FREEZE_STAGE_A.txt").write_text("\n".join(
        f"{sha(protocol / n)}  {n}" for n in ("PROJECT_PLAN.md", "protocol_config.json", "prefreeze_network_check.csv")) + "\n")
    (protocol / "FREEZE_STAGE_A_ADDENDUM.txt").write_text(
        f"{sha(protocol / 'AMENDMENT_A1_braf_classification.md')}  AMENDMENT_A1_braf_classification.md\n")
    seeds = project / "results" / "melanoma_stageA_seeds"
    m0 = {"ARID2": 0.3, "PTEN": 0.3, "TP53": 0.4}
    m1 = {"ARID2": 0.25, "PTEN": 0.25, "TP53": 0.3, "BRAF": 0.2}
    rows = [{"arm": a, "gene": g, "raw_weight": w, "normalized_weight": w, "in_STRING_v12": "yes"}
            for a, ws in (("M0", m0), ("M1", m1)) for g, w in ws.items()]
    write_rows(seeds / "primary_metastasis" / "seeds_frozen_string_checked.csv", rows)
    write_rows(seeds / "S8_uniform_weights" / "seeds_frozen_string_checked.csv",
               [{"arm": "M0", "gene": g, "raw_weight": 1, "normalized_weight": 1 / 3, "in_STRING_v12": "yes"} for g in m0])
    prep = project / "results" / "melanoma_stageB_prep"
    prep.mkdir(parents=True)
    (prep / "stageB_exposure_anchor_summary.json").write_text(json.dumps({
        "X": ["Dacarbazine", "Tioguanine"],
        "anchors": [{"disease": "ACC", "c": "46/108", "u": "0/1", "rho_L_count": 0.822},
                    {"disease": "Breast", "c": "14/108", "u": "13/94", "rho_L_count": 0.577}]}))
    rel = "results/melanoma_stageB_prep/stageB_exposure_anchor_summary.json"
    (protocol / "FREEZE_STAGE_B.txt").write_text(f"{sha(project / rel)}  {rel}\n")
    return repo, project, string_dir, drugs, u108, edges


@pytest.fixture
def world(tmp_path, monkeypatch):
    repo, project, string_dir, drugs, u108, edges = build_world(tmp_path)
    monkeypatch.setattr(ma, "REPO_ROOT", repo)
    return {"repo": repo, "project": project, "string": string_dir, "drugs": drugs, "u108": u108, "edges": edges,
            "tmp": tmp_path}


def engine_for(w, threshold=400):
    assoc = {}
    for e in w["edges"]:
        assoc.setdefault(e["drug"], set()).add(e["gene"])
    s = w["string"]
    return ma.Engine(s / "9606.protein.info.v12.0.txt.gz", s / "9606.protein.links.v12.0.txt.gz", assoc,
                     {"ARID2", "PTEN", "TP53", "BRAF"}, threshold, 0.4), assoc


# --------------------------------------------------------------- engine
def test_c_subset_full_equals_c_raw_and_empty_equals_mu(world):
    eng, _ = engine_for(world)
    gv = eng.gene_variants([{"ARID2": 0.5, "PTEN": 0.5}])
    for v in ma.VARIANTS:
        full = eng.c_subset(gv[v], eng.A)
        assert np.allclose(full, eng.c_raw(gv[v]))
        empty = eng.c_subset(gv[v], eng.A.multiply(0).tocsr())
        assert np.allclose(empty, eng.background(gv[v])[None, :])


def test_variants_in_unit_interval_and_restart_renormalised(world):
    eng, _ = engine_for(world)
    a = eng.gene_variants([{"ARID2": 1.0, "PTEN": 1.0}])
    b = eng.gene_variants([{"ARID2": 0.5, "PTEN": 0.5}])
    for v in ma.VARIANTS:
        assert a[v].min() >= 0 and a[v].max() <= 1
        assert np.allclose(a[v], b[v])


def test_threshold_700_drops_edges(world):
    e400, _ = engine_for(world, 400)
    e700, _ = engine_for(world, 700)
    assert e700.graph.source_edge_rows < e400.graph.source_edge_rows
    assert e700.graph.node_names == e400.graph.node_names


def test_null_streams_reproducible_arm_specific_and_matched(world):
    eng, _ = engine_for(world)
    w = {"ARID2": 0.3, "PTEN": 0.3, "TP53": 0.4}
    r1, _ = ma.null_seed_sets(eng, w, 0, [11, 12], 5, 3)
    r2, _ = ma.null_seed_sets(eng, w, 0, [11, 12], 5, 3)
    r3, _ = ma.null_seed_sets(eng, w, 1, [11, 12], 5, 3)
    assert r1 == r2
    assert [x["null_seed"] for x in r1] != [x["null_seed"] for x in r3]
    assert sorted({x["replicate"] for x in r1}) == list(range(10))
    assert all(x["null_seed"] not in w for x in r1)
    assert all(x["matched_seed_degree_bin"] == x["null_seed_degree_bin"] for x in r1)
    for rep in range(10):
        draw = [x for x in r1 if x["replicate"] == rep]
        assert sorted(x["weight"] for x in draw) == sorted(w.values())


# ------------------------------------------------------------ gates/I-O
def test_freeze_gate_detects_modified_file(world):
    ma.check_freezes(world["project"])
    p = world["project"] / "results" / "melanoma_stageB_prep" / "stageB_exposure_anchor_summary.json"
    p.write_text(p.read_text() + " ")
    with pytest.raises(RuntimeError):
        ma.check_freezes(world["project"])


def test_unmapped_seed_rejected(world):
    p = world["project"] / "results" / "melanoma_stageA_seeds" / "S8_uniform_weights" / "seeds_frozen_string_checked.csv"
    p.write_text(p.read_text().replace("yes", "no", 1))
    with pytest.raises(RuntimeError):
        ma.load_inputs(world["project"])


# ------------------------------------------------------------ end to end
def test_end_to_end_run(world):
    out = world["tmp"] / "out"
    assert ma.run(world["project"], world["string"], out) == 0
    for name in ("drug_statistics.csv", "audit_metrics.csv", "seed_nonseed_decomposition.csv", "pooled_bh_992.csv",
                 "H6_delta_z_124.csv", "H6_delta_z_108.csv", "loso_M0.csv", "S10_parameter_sensitivity.csv",
                 "seed_branches.csv", "decisions.json", "run_manifest.json", "null_C_raw_M0.npz", "null_seed_sets_M1.csv"):
        assert (out / name).exists(), name
    stats_rows = list(csv.DictReader((out / "drug_statistics.csv").open()))
    assert len(stats_rows) == 2 * (40 + 35) * 4
    for r in stats_rows:
        p = float(r["p_upper"])
        assert 1 / 21 - 1e-12 <= p <= 1
        assert int(r["null_draws"]) == 20
    # percentiles live in [0, 1] within each universe
    for u, n in (("124", 40), ("108", 35)):
        vals = [float(r["percentile"]) for r in stats_rows if r["universe"] == u and r["arm"] == "M0"
                and r["variant"] == "column_minmax"]
        assert len(vals) == n and min(vals) == 0 and max(vals) == 1
    decisions = json.loads((out / "decisions.json").read_text())
    assert set(decisions) >= {"H4", "H5", "H6", "X", "H5b"}
    assert decisions["H4"]["124"]["status"] in {"pass", "partial_recovery", "fail"}
    assert decisions["H6"]["124"]["verdict"] in {"SUPPORTED", "PARTIAL", "NOT_SUPPORTED", "NOT_EVALUABLE"}
    assert decisions["H5"]["124"]["verdict"] == "PENDING_124_ANCHORS"
    pooled = list(csv.DictReader((out / "pooled_bh_992.csv").open()))
    assert len(pooled) == 2 * 4 * 40
    loso = list(csv.DictReader((out / "loso_M0.csv").open()))
    assert len(loso) == 3 * 4
    s10 = list(csv.DictReader((out / "S10_parameter_sensitivity.csv").open()))
    assert {r["factor"] for r in s10} == {"alpha", "k", "string_threshold"}
    branches = list(csv.DictReader((out / "seed_branches.csv").open()))
    assert {r["branch"] for r in branches} == {"S8_uniform_weights"}
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["n_null"] == 20 and len(manifest["null"]["M0"]["rng"]) == 5


def test_end_to_end_is_deterministic(world):
    a, b = world["tmp"] / "a", world["tmp"] / "b"
    ma.run(world["project"], world["string"], a)
    ma.run(world["project"], world["string"], b)
    for name in ("drug_statistics.csv", "audit_metrics.csv", "decisions.json", "loso_M0.csv"):
        assert (a / name).read_bytes() == (b / name).read_bytes(), name


def test_z_matches_manual_computation(world):
    out = world["tmp"] / "out"
    ma.run(world["project"], world["string"], out)
    raw = np.load(out / "null_C_raw_M0.npz")["column_minmax"]
    obs_raw = np.load(out / "observed_C_raw_M0.npz")["column_minmax"]
    from analysis.method_strengthening import percentile_columns
    null_pct = percentile_columns(raw)
    obs_pct = percentile_columns(obs_raw)[:, 0]
    rows = [r for r in csv.DictReader((out / "drug_statistics.csv").open())
            if r["arm"] == "M0" and r["universe"] == "124" and r["variant"] == "column_minmax"]
    drugs = sorted(world["drugs"])
    for r in rows:
        i = drugs.index(r["drug"])
        obs = obs_pct[i]
        assert float(r["percentile"]) == pytest.approx(obs, abs=1e-11)
        b = int(np.count_nonzero(null_pct[i] >= obs))
        assert int(r["exceedances"]) == b
        if r["z"] not in ("", "nan"):
            z = (obs - null_pct[i].mean()) / null_pct[i].std(ddof=1)
            assert float(r["z"]) == pytest.approx(z, rel=1e-9)


def test_benchmark_writes_no_scientific_output(world):
    out = world["tmp"] / "bench"
    assert ma.run(world["project"], world["string"], out, benchmark=True) == 0
    assert sorted(p.name for p in out.iterdir()) == ["benchmark.json"]
