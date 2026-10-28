"""Gold-standard tests for stage-A melanoma seed generation."""
import csv
import json
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis"))
import melanoma_seeds as ms  # noqa: E402


def rec(cls, prot, gene="BRAF"):
    return {"Hugo_Symbol": gene, "Variant_Classification": cls, "HGVSp_Short": prot}


# ------------------------------------------------------------------ BH
def test_bh_matches_hand_calculation():
    p = [0.01, 0.04, 0.03, 0.5]
    # sorted 0.01,0.03,0.04,0.5 -> raw 0.04,0.06,0.0533,0.5 -> cummin from top
    q = ms.benjamini_hochberg(p)
    assert q == pytest.approx([0.04, 0.0533333333, 0.0533333333, 0.5])


def test_bh_backward_cummin_and_cap():
    q = ms.benjamini_hochberg([0.9, 0.95, 1.0])
    assert q == pytest.approx([1.0, 1.0, 1.0])
    q = ms.benjamini_hochberg([0.02, 0.02, 0.02])
    assert q == pytest.approx([0.02, 0.02, 0.02])


def test_bh_placeholder_p1_counts_in_family():
    q_small = ms.benjamini_hochberg([0.001])[0]
    q_family = ms.benjamini_hochberg([0.001] + [1.0] * 197)[0]
    assert q_small == pytest.approx(0.001)
    assert q_family == pytest.approx(0.198)


# ----------------------------------------------------------------- BRAF
@pytest.mark.parametrize("records,expected", [
    ([], ("REF", "", "")),
    ([rec("Missense_Mutation", "p.V600E")], ("EK", "E", "")),
    ([rec("Missense_Mutation", "p.V600K")], ("EK", "K", "")),
    ([rec("Missense_Mutation", "p.V600K"), rec("Silent", "p.V600=")], ("EK", "K", "")),
    ([rec("Missense_Mutation", "p.V600R")], ("EXCLUDE_OTHER", "p.V600R", "")),
    ([rec("In_Frame_Del", "p.V600_K601delinsE")], ("EXCLUDE_OTHER", "p.V600_K601delinsE", "")),
    ([rec("Missense_Mutation", "p.V600E"), rec("Missense_Mutation", "p.K601E")], ("EXCLUDE_OTHER", "p.K601E;p.V600E", "")),
    ([rec("Missense_Mutation", "p.V600E"), rec("Missense_Mutation", "p.V600K")], ("EXCLUDE_OTHER", "p.V600E;p.V600K", "")),
    ([rec("Silent", "p.S363=")], ("EXCLUDE_SILENT_ONLY", "Silent:p.S363=", "")),
    # Amendment A1: other-codon missense may co-occur with V600E/K
    ([rec("Missense_Mutation", "p.V600E"), rec("Missense_Mutation", "p.M438I")], ("EK", "E", "p.M438I")),
    ([rec("Missense_Mutation", "p.V600E"), rec("Missense_Mutation", "p.P318L"), rec("Missense_Mutation", "p.P318S")],
     ("EK", "E", "p.P318L;p.P318S")),
    ([rec("Missense_Mutation", "p.V600K"), rec("Missense_Mutation", "p.G606V")], ("EK", "K", "p.G606V")),
    ([rec("Missense_Mutation", "p.V600E"), rec("In_Frame_Del", "p.N486_P490delinsK")],
     ("EXCLUDE_OTHER", "p.N486_P490delinsK;p.V600E", "")),
])
def test_classify_braf_A1(records, expected):
    assert ms.classify_braf(records) == expected


@pytest.mark.parametrize("records,expected", [
    ([rec("Missense_Mutation", "p.V600E"), rec("Missense_Mutation", "p.M438I")], ("EXCLUDE_OTHER", "p.M438I;p.V600E", "")),
    ([rec("Missense_Mutation", "p.V600E"), rec("Silent", "p.V600=")], ("EK", "E", "")),
])
def test_classify_braf_strict(records, expected):
    assert ms.classify_braf(records, "strict") == expected


def test_classify_braf_rejects_unknown_mode():
    with pytest.raises(ValueError):
        ms.classify_braf([], "loose")


# --------------------------------------------------------------- cohort
def sample(pid, sid, stype="Metastasis"):
    return {"PATIENT_ID": pid, "SAMPLE_ID": sid, "SAMPLE_TYPE": stype}


def test_cohort_dedupes_to_first_eligible_sample_per_patient():
    samples = [sample("P1", "P1-07"), sample("P1", "P1-06"), sample("P2", "P2-06"),
               sample("P3", "P3-01", "Primary"), sample("P4", "P4-06")]
    patients = [{"PATIENT_ID": p, "HISTORY_NEOADJUVANT_TRTYN": "No"} for p in ("P1", "P2", "P3", "P4")]
    sequenced = {"P1-06", "P1-07", "P2-06", "P3-01"}
    rppa = {"P1-06", "P1-07", "P2-06", "P3-01", "P4-06"}
    braf = {"P1-06": [rec("Missense_Mutation", "p.V600E")], "P1-07": [rec("Missense_Mutation", "p.V600E")]}
    chosen, log = ms.build_cohort(samples, patients, sequenced, rppa, braf, "Metastasis")
    assert chosen == {"P1": {"sample_id": "P1-06", "group": "EK", "allele": "E"},
                      "P2": {"sample_id": "P2-06", "group": "REF", "allele": ""}}
    by_id = {r["sample_id"]: r for r in log}
    assert by_id["P1-07"]["selected"] == "no" and by_id["P1-07"]["eligible"] == "yes"
    assert "sample_type=Primary" in by_id["P3-01"]["exclusion_reasons"]
    assert "not_in_sequenced_case_list" in by_id["P4-06"]["exclusion_reasons"]


def test_cohort_requires_rppa_on_same_sample():
    samples = [sample("P1", "P1-06")]
    chosen, log = ms.build_cohort(samples, [], {"P1-06"}, set(), {}, "Metastasis")
    assert chosen == {} and "no_rppa_column" in log[0]["exclusion_reasons"]


def test_neoadjuvant_filter():
    samples = [sample("P1", "P1-06"), sample("P2", "P2-06")]
    patients = [{"PATIENT_ID": "P1", "HISTORY_NEOADJUVANT_TRTYN": "No"},
                {"PATIENT_ID": "P2", "HISTORY_NEOADJUVANT_TRTYN": "Yes"}]
    chosen, _ = ms.build_cohort(samples, patients, {"P1-06", "P2-06"}, {"P1-06", "P2-06"}, {},
                                "Metastasis", "neoadjuvant_no")
    assert set(chosen) == {"P1"}


# ------------------------------------------------------------- mutation
def write_maf(path, rows):
    with path.open("w", newline="") as handle:
        handle.write("#version 2.4\n")
        w = csv.DictWriter(handle, fieldnames=["Hugo_Symbol", "Variant_Classification", "HGVSp_Short",
                                               "Tumor_Sample_Barcode"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def test_mutation_thresholds_and_once_per_patient(tmp_path):
    ek = {f"S{i}": f"P{i}" for i in range(100)}
    rows = []
    for i in range(10):  # TP53: 10/100 -> eligible, R = 0.1/0.2 = 0.5
        rows.append({"Hugo_Symbol": "TP53", "Variant_Classification": "Missense_Mutation",
                     "HGVSp_Short": "x", "Tumor_Sample_Barcode": f"S{i}"})
        rows.append({"Hugo_Symbol": "TP53", "Variant_Classification": "Nonsense_Mutation",
                     "HGVSp_Short": "y", "Tumor_Sample_Barcode": f"S{i}"})
    for i in range(9):  # PTEN: 9 -> not eligible
        rows.append({"Hugo_Symbol": "PTEN", "Variant_Classification": "Missense_Mutation",
                     "HGVSp_Short": "x", "Tumor_Sample_Barcode": f"S{i}"})
    for i in range(20):  # NF1 silent only -> not counted
        rows.append({"Hugo_Symbol": "NF1", "Variant_Classification": "Silent",
                     "HGVSp_Short": "x", "Tumor_Sample_Barcode": f"S{i}"})
    rows.append({"Hugo_Symbol": "TP53", "Variant_Classification": "Missense_Mutation",
                 "HGVSp_Short": "x", "Tumor_Sample_Barcode": "OTHER"})
    write_maf(tmp_path / "m.txt", rows)
    out = {r["gene"]: r for r in ms.mutation_statistics(tmp_path / "m.txt", ek)}
    assert out["TP53"]["mutated_patients"] == 10 and out["TP53"]["G"] == 1.0
    assert out["TP53"]["R"] == pytest.approx(0.5)
    assert out["PTEN"]["eligible"] == "no" and out["PTEN"]["R"] == 0.0
    assert out["NF1"]["mutated_patients"] == 0


# ----------------------------------------------------------------- RPPA
def write_rppa(path, features, case, ref):
    samples = list(case) + list(ref)
    with path.open("w", newline="") as handle:
        handle.write("Composite.Element.REF\t" + "\t".join(samples) + "\n")
        for name, case_vals, ref_vals in features:
            handle.write(name + "\t" + "\t".join(case_vals + ref_vals) + "\n")


def make_features(n_total=198):
    case = [f"C{i}" for i in range(60)]
    ref = [f"R{i}" for i in range(60)]
    up = ([str(1.0 + i * 0.001) for i in range(60)], [str(-1.0 + i * 0.001) for i in range(60)])
    feats = [("GENEA|a", *up),                     # selected, diff 2.0 -> L = 1.0
             ("GENEB|b1", [str(0.6 + i * 1e-4) for i in range(60)], [str(i * 1e-4) for i in range(60)]),
             ("GENEB|b2", [str(0.8 + i * 1e-4) for i in range(60)], [str(i * 1e-4) for i in range(60)]),
             ("GENEC GENED|cd", *up),               # multi-gene
             ("GENEE|const", ["0"] * 60, ["0"] * 60),
             ("GENEF|sparse", ["1"] * 49 + ["NA"] * 11, ["0"] * 60)]
    while len(feats) < n_total:
        k = len(feats)
        feats.append((f"FILL{k}|f", [str((i * 7 + k) % 13) for i in range(60)],
                      [str((i * 5 + k) % 13) for i in range(60)]))
    return case, ref, feats


def test_rppa_rules(tmp_path):
    case, ref, feats = make_features()
    write_rppa(tmp_path / "r.txt", feats, case, ref)
    rows, selected = ms.rppa_statistics(tmp_path / "r.txt", case, ref)
    by = {r["feature"]: r for r in rows}
    assert by["GENEC GENED|cd"]["not_tested_reason"] == "multi_gene_annotation"
    assert by["GENEE|const"]["not_tested_reason"] == "constant_feature"
    assert by["GENEF|sparse"]["not_tested_reason"] == "fewer_than_50_finite_in_a_group"
    assert by["GENEC GENED|cd"]["p_two_sided"] == 1.0
    assert selected["GENEA"]["L"] == pytest.approx(min(abs(by["GENEA|a"]["median_difference"]) / 2, 1))
    assert selected["GENEB"]["feature"] == "GENEB|b2"   # larger effect wins
    assert "GENEC" not in selected


def test_rppa_split_mode_divides_L(tmp_path):
    case, ref, feats = make_features()
    write_rppa(tmp_path / "r.txt", feats, case, ref)
    _, selected = ms.rppa_statistics(tmp_path / "r.txt", case, ref, "split_equally")
    assert selected["GENEC"]["L"] == pytest.approx(selected["GENEA"]["L"] / 2)


def test_rppa_family_size_enforced(tmp_path):
    case, ref, feats = make_features(10)
    write_rppa(tmp_path / "r.txt", feats, case, ref)
    with pytest.raises(ValueError):
        ms.rppa_statistics(tmp_path / "r.txt", case, ref)


def test_rppa_fraction_rule(tmp_path):
    case = [f"C{i}" for i in range(120)]
    ref = [f"R{i}" for i in range(60)]
    feats = [("GENEA|a", [str(i) for i in range(55)] + ["NA"] * 65, [str(-i) for i in range(60)])]
    feats += [(f"F{k}|f", ["1"] * 120, ["0"] * 60) for k in range(197)]
    write_rppa(tmp_path / "r.txt", feats, case, ref)
    rows, _ = ms.rppa_statistics(tmp_path / "r.txt", case, ref)
    assert rows[0]["not_tested_reason"] == "finite_fraction_below_0.5_in_a_group"


# ---------------------------------------------------------------- seeds
def test_seed_weights_and_arms():
    mutation = [{"gene": g, "G": 0.0, "R": 0.0, "mutated_patients": 0, "fraction": 0.0} for g in ms.MUTATION_GENES]
    mutation = [dict(m, G=1.0, R=1.0, mutated_patients=106, fraction=1.0) if m["gene"] == "BRAF" else m for m in mutation]
    mutation = [dict(m, G=1.0, R=0.5, mutated_patients=10, fraction=0.1) if m["gene"] == "TP53" else m for m in mutation]
    rppa = {"MAP2K1": {"feature": "MAP2K1|MEK1", "effect": 1.0, "q": 0.001, "direction": "up", "L": 0.5},
            "XYZ": {"feature": "XYZ|x", "effect": 4.0, "q": 0.001, "direction": "down", "L": 1.0}}
    rows = {r["gene"]: r for r in ms.build_seed_rows(mutation, rppa)}
    assert rows["BRAF"]["raw_weight"] == pytest.approx(0.5)
    assert rows["BRAF"]["include_M0"] == "no" and rows["BRAF"]["include_M1"] == "yes"
    assert rows["MAP2K1"]["include_M0"] == "no" and rows["MAP2K1"]["include_M1"] == "yes"
    assert rows["MAP2K2"]["include_M1"] == "no"           # not eligible -> not forced in
    assert rows["TP53"]["raw_weight"] == pytest.approx(0.40)
    assert rows["XYZ"]["raw_weight"] == pytest.approx(0.10)
    m0 = ms.arm_weights(list(rows.values()), "M0")
    assert {r["gene"] for r in m0} == {"TP53", "XYZ"}
    assert math.fsum(r["normalized_weight"] for r in m0) == pytest.approx(1.0)
    m1 = ms.arm_weights(list(rows.values()), "M1")
    assert {r["gene"] for r in m1} == {"TP53", "XYZ", "BRAF", "MAP2K1"}
    uni = ms.arm_weights(list(rows.values()), "M0", uniform=True)
    assert all(r["normalized_weight"] == pytest.approx(0.5) for r in uni)


# ----------------------------------------------------------- freeze gate
def test_freeze_gate_rejects_modified_files(tmp_path):
    for name, text in [("PROJECT_PLAN.md", "plan"), ("prefreeze_network_check.csv", "x")]:
        (tmp_path / name).write_text(text)
    (tmp_path / "protocol_config.json").write_text(json.dumps(
        {"status": "STAGE_A_FROZEN", "execution_allowed": "stage_A_only"}))
    lines = [f"{ms.sha256(tmp_path / n)}  {n}" for n in
             ("PROJECT_PLAN.md", "protocol_config.json", "prefreeze_network_check.csv")]
    (tmp_path / "FREEZE_STAGE_A.txt").write_text("header\n" + "\n".join(lines) + "\n")
    assert len(ms.check_stage_a_freeze(tmp_path)) == 3
    (tmp_path / "PROJECT_PLAN.md").write_text("plan edited")
    with pytest.raises(RuntimeError):
        ms.check_stage_a_freeze(tmp_path)


def test_freeze_gate_rejects_unauthorised_config(tmp_path):
    for name in ("PROJECT_PLAN.md", "prefreeze_network_check.csv"):
        (tmp_path / name).write_text("x")
    (tmp_path / "protocol_config.json").write_text(json.dumps({"status": "DRAFT", "execution_allowed": False}))
    lines = [f"{ms.sha256(tmp_path / n)}  {n}" for n in
             ("PROJECT_PLAN.md", "protocol_config.json", "prefreeze_network_check.csv")]
    (tmp_path / "FREEZE_STAGE_A.txt").write_text("\n".join(lines) + "\n")
    with pytest.raises(RuntimeError):
        ms.check_stage_a_freeze(tmp_path)


def test_module_does_not_touch_stage_b_inputs():
    source = Path(ms.__file__).read_text()
    for token in ("bindex", "protein.links", "random_walk", "propagat" + "e(", "C_MEL"):
        assert token not in source
