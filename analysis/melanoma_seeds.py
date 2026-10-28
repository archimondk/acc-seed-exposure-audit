"""Stage-A seed generation for protocol melanoma_R1_v0.3.

Implements PROJECT_PLAN.md sections 1.2, 2.1-2.3 and the seed-producing
sensitivity branches of section 7 that are fully specified (S2, S3, S7, S9).
It never reads the drug-gene association network, STRING graph, drug lists
or any drug score: propagation, scoring, nulls and the exposure map belong to
stage B and are not implemented here.

Usage:
    python -I melanoma_seeds.py --data <raw_dir> --protocol <protocol_dir> --out <out_dir>
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import scipy
from scipy import stats

ANALYSIS_VERSION = "melanoma-seeds-stageA-v1"
PROTOCOL_ID = "melanoma_R1_v0.3"
MUTATION_GENES = (
    "ARID2", "BRAF", "CDKN2A", "DDX3X", "IDH1", "MAP2K1", "NF1",
    "NRAS", "PPP6C", "PTEN", "RAC1", "RB1", "TP53",
)
FORCED_EXCLUSIONS = ("BRAF", "MAP2K1", "MAP2K2")
NONSYNONYMOUS = frozenset({
    "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
    "Splice_Site", "Translation_Start_Site", "Frame_Shift_Del",
    "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
})
EK_PROTEIN = {"p.V600E": "E", "p.V600K": "K"}
MUT_MIN_PATIENTS = 10
MUT_MIN_FRACTION = 0.05
MUT_RECURRENCE_SCALE = 0.20
RPPA_FAMILY_N = 198
RPPA_MIN_N = 50
RPPA_MIN_FRACTION = 0.5
RPPA_Q = 0.01
RPPA_MIN_EFFECT = 0.5
RPPA_EFFECT_SCALE = 2.0
W_G, W_R, W_L = 0.30, 0.20, 0.10
MISSING_TOKENS = frozenset({"", "NA", "NaN", "nan", "null", "NULL"})


# ---------------------------------------------------------------- utilities
def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_cbio_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        lines = (line for line in handle if not line.startswith("#"))
        rows = list(csv.DictReader(lines, delimiter="\t"))
    if not rows:
        raise ValueError(f"no data rows in {path}")
    return rows


def read_case_list(path: Path) -> set[str]:
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("case_list_ids:"):
            return {s.strip() for s in line.split(":", 1)[1].split("\t") if s.strip()}
    raise ValueError(f"no case_list_ids in {path}")


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f"{v:.12g}" if isinstance(v, float) else v) for k, v in row.items()})


def benjamini_hochberg(p_values: Sequence[float]) -> list[float]:
    """Step-up BH with backward cumulative minimum, capped at 1."""
    m = len(p_values)
    order = sorted(range(m), key=lambda i: (p_values[i], i))
    q = [0.0] * m
    running = 1.0
    for rank in range(m, 0, -1):
        index = order[rank - 1]
        running = min(running, p_values[index] * m / rank)
        q[index] = min(running, 1.0)
    return q


def parse_value(token: str) -> float | None:
    token = token.strip()
    if token in MISSING_TOKENS:
        return None
    value = float(token)
    return value if math.isfinite(value) else None


# ------------------------------------------------------------- BRAF status
_POSITION = re.compile(r"^p\.[A-Z*]?(\d+)")


def _protein_start(protein: str) -> int | None:
    match = _POSITION.match(protein)
    return int(match.group(1)) if match else None


def classify_braf(records: Iterable[Mapping[str, str]], mode: str = "A1") -> tuple[str, str, str]:
    """Return (status, allele_or_detail, co_occurring) for one sample's BRAF records.

    mode 'A1' (Amendment A1, primary): E/K when exactly one of p.V600E/p.V600K
    is present, no other non-synonymous annotation at codon 600/601 and no
    delins; non-synonymous BRAF variants at other codons may co-occur and are
    returned in co_occurring. mode 'strict' (descriptive sensitivity): any
    other non-synonymous BRAF record excludes the sample.
    'REF' = no BRAF record of any class; 'EXCLUDE_SILENT_ONLY' = synonymous
    or other non-listed BRAF records only. Synonymous co-annotation alongside
    E/K is never a conflict.
    """
    if mode not in ("A1", "strict"):
        raise ValueError(mode)
    records = list(records)
    if not records:
        return "REF", "", ""
    nonsyn = [r for r in records if r["Variant_Classification"] in NONSYNONYMOUS]
    if not nonsyn:
        detail = ";".join(sorted({f"{r['Variant_Classification']}:{r.get('HGVSp_Short', '')}" for r in records}))
        return "EXCLUDE_SILENT_ONLY", detail, ""
    proteins = {r.get("HGVSp_Short", "").strip() for r in nonsyn}
    detail = ";".join(sorted(proteins))
    ek = proteins & set(EK_PROTEIN)
    others = proteins - ek
    if len(ek) != 1:
        return "EXCLUDE_OTHER", detail, ""
    if mode == "strict":
        return ("EK", EK_PROTEIN[next(iter(ek))], "") if not others else ("EXCLUDE_OTHER", detail, "")
    blocking = {p for p in others if "delins" in p or _protein_start(p) in (600, 601)}
    if blocking:
        return "EXCLUDE_OTHER", detail, ""
    return "EK", EK_PROTEIN[next(iter(ek))], ";".join(sorted(others))


def braf_records_by_sample(maf_path: Path) -> dict[str, list[dict[str, str]]]:
    out: dict[str, list[dict[str, str]]] = defaultdict(list)
    with maf_path.open("r", encoding="utf-8-sig", newline="") as handle:
        lines = (line for line in handle if not line.startswith("#"))
        for row in csv.DictReader(lines, delimiter="\t"):
            if row["Hugo_Symbol"] == "BRAF":
                out[row["Tumor_Sample_Barcode"].strip()].append(row)
    return out


# ------------------------------------------------------------------ cohort
def build_cohort(
    sample_rows: Sequence[Mapping[str, str]],
    patient_rows: Sequence[Mapping[str, str]],
    sequenced: set[str],
    rppa_samples: set[str],
    braf_by_sample: Mapping[str, list[dict[str, str]]],
    sample_type_filter: str | None,
    patient_filter: str | None = None,
    braf_mode: str = "A1",
) -> tuple[dict[str, dict[str, str]], list[dict[str, object]]]:
    """Return ({patient: chosen sample info}, full sample log)."""
    patient_info = {r["PATIENT_ID"].strip(): r for r in patient_rows}
    log: list[dict[str, object]] = []
    eligible: list[tuple[str, str, str, str]] = []
    for row in sorted(sample_rows, key=lambda r: (r["PATIENT_ID"].strip(), r["SAMPLE_ID"].strip())):
        patient, sample = row["PATIENT_ID"].strip(), row["SAMPLE_ID"].strip()
        stype = row["SAMPLE_TYPE"].strip()
        if sample in sequenced:
            status, detail, co_occurring = classify_braf(braf_by_sample.get(sample, []), braf_mode)
        else:
            status, detail, co_occurring = "NOT_SEQUENCED", "", ""
        reasons = []
        if sample_type_filter is not None and stype != sample_type_filter:
            reasons.append(f"sample_type={stype}")
        if sample not in sequenced:
            reasons.append("not_in_sequenced_case_list")
        if sample not in rppa_samples:
            reasons.append("no_rppa_column")
        if status not in ("EK", "REF"):
            reasons.append(f"braf_{status}")
        if patient_filter == "neoadjuvant_no":
            value = patient_info.get(patient, {}).get("HISTORY_NEOADJUVANT_TRTYN", "").strip()
            if value != "No":
                reasons.append(f"neoadjuvant={value or 'missing'}")
        log.append({"patient_id": patient, "sample_id": sample, "sample_type": stype,
                    "braf_status": status, "braf_detail": detail,
                    "braf_co_occurring_other_codon": co_occurring,
                    "eligible": "yes" if not reasons else "no",
                    "exclusion_reasons": ";".join(reasons), "selected": "no"})
        if not reasons:
            eligible.append((patient, sample, status, detail))
    chosen: dict[str, dict[str, str]] = {}
    for patient, sample, status, detail in eligible:  # already sorted by (patient, sample)
        if patient not in chosen:
            chosen[patient] = {"sample_id": sample, "group": status, "allele": detail}
    selected = {v["sample_id"] for v in chosen.values()}
    for entry in log:
        if entry["sample_id"] in selected:
            entry["selected"] = "yes"
    return chosen, log


# ---------------------------------------------------------------- mutation
def mutation_statistics(maf_path: Path, ek_samples: Mapping[str, str]) -> list[dict[str, object]]:
    """ek_samples: {sample_id: patient_id} for the E/K group only."""
    mutated: dict[str, set[str]] = defaultdict(set)
    genes = set(MUTATION_GENES)
    with maf_path.open("r", encoding="utf-8-sig", newline="") as handle:
        lines = (line for line in handle if not line.startswith("#"))
        for row in csv.DictReader(lines, delimiter="\t"):
            sample = row["Tumor_Sample_Barcode"].strip()
            if sample in ek_samples and row["Hugo_Symbol"] in genes and row["Variant_Classification"] in NONSYNONYMOUS:
                mutated[row["Hugo_Symbol"]].add(ek_samples[sample])
    denominator = len(set(ek_samples.values()))
    out = []
    for gene in MUTATION_GENES:
        count = len(mutated.get(gene, ()))
        fraction = count / denominator if denominator else math.nan
        eligible = count >= MUT_MIN_PATIENTS and fraction >= MUT_MIN_FRACTION
        out.append({"gene": gene, "mutated_patients": count, "denominator": denominator,
                    "fraction": fraction, "eligible": "yes" if eligible else "no",
                    "G": 1.0 if eligible else 0.0,
                    "R": min(fraction / MUT_RECURRENCE_SCALE, 1.0) if eligible else 0.0})
    return out


# -------------------------------------------------------------------- RPPA
def rppa_statistics(
    rppa_path: Path,
    case_samples: Sequence[str],
    ref_samples: Sequence[str],
    multi_gene_mode: str = "single_gene_only",
) -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    with rppa_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        raw = [row for row in reader if row]
    if len(raw) != RPPA_FAMILY_N:
        raise ValueError(f"expected {RPPA_FAMILY_N} RPPA features, found {len(raw)}")
    column = {name: i for i, name in enumerate(header)}
    case_idx = [column[s] for s in case_samples]
    ref_idx = [column[s] for s in ref_samples]
    n_case_group, n_ref_group = len(case_idx), len(ref_idx)
    rows: list[dict[str, object]] = []
    p_values: list[float] = []
    for record in raw:
        feature = record[0]
        genes = tuple(g.strip().upper() for g in feature.split("|", 1)[0].split() if g.strip())
        a = np.asarray([v for v in (parse_value(record[i]) for i in case_idx) if v is not None], dtype=float)
        b = np.asarray([v for v in (parse_value(record[i]) for i in ref_idx) if v is not None], dtype=float)
        reason = ""
        if not genes:
            reason = "no_gene_annotation"
        elif len(genes) > 1 and multi_gene_mode == "single_gene_only":
            reason = "multi_gene_annotation"
        elif len(a) < RPPA_MIN_N or len(b) < RPPA_MIN_N:
            reason = "fewer_than_50_finite_in_a_group"
        elif len(a) < RPPA_MIN_FRACTION * n_case_group or len(b) < RPPA_MIN_FRACTION * n_ref_group:
            reason = "finite_fraction_below_0.5_in_a_group"
        elif np.ptp(np.concatenate([a, b])) == 0:
            reason = "constant_feature"
        if reason:
            p, diff = 1.0, math.nan
        else:
            diff = float(np.median(a) - np.median(b))
            p = float(stats.mannwhitneyu(a, b, alternative="two-sided",
                                         use_continuity=True, method="asymptotic").pvalue)
        rows.append({"feature": feature, "genes": ";".join(genes), "case_n": len(a), "ref_n": len(b),
                     "median_difference": diff, "p_two_sided": p, "not_tested_reason": reason})
        p_values.append(p)
    q_values = benjamini_hochberg(p_values)
    selected: dict[str, dict[str, object]] = {}
    for row, q in zip(rows, q_values, strict=True):
        row["q_bh_198"] = q
        ok = (not row["not_tested_reason"]) and q < RPPA_Q and abs(float(row["median_difference"])) >= RPPA_MIN_EFFECT
        row["selected"] = "yes" if ok else "no"
        if not ok:
            continue
        effect = abs(float(row["median_difference"]))
        genes = str(row["genes"]).split(";")
        share = len(genes) if multi_gene_mode == "split_equally" else 1
        for gene in genes:
            lineage = min(effect / RPPA_EFFECT_SCALE, 1.0) / share
            current = selected.get(gene)
            # keep largest L; ties broken by feature name (lexicographically first)
            if (current is None or lineage > current["L"]
                    or (lineage == current["L"] and row["feature"] < current["feature"])):
                selected[gene] = {"feature": row["feature"], "effect": effect, "q": q,
                                  "direction": "up" if float(row["median_difference"]) > 0 else "down",
                                  "L": lineage}
    return rows, selected


# ------------------------------------------------------------------- seeds
def build_seed_rows(
    mutation_rows: Sequence[Mapping[str, object]],
    rppa_selected: Mapping[str, Mapping[str, object]],
) -> list[dict[str, object]]:
    mut = {r["gene"]: r for r in mutation_rows}
    candidates = sorted(set(MUTATION_GENES) | set(rppa_selected) | set(FORCED_EXCLUSIONS))
    rows = []
    for gene in candidates:
        g = float(mut[gene]["G"]) if gene in mut else 0.0
        r = float(mut[gene]["R"]) if gene in mut else 0.0
        rp = rppa_selected.get(gene)
        lineage = float(rp["L"]) if rp else 0.0
        raw = W_G * g + W_R * r + W_L * lineage
        forced = gene in FORCED_EXCLUSIONS
        reasons = []
        if raw <= 0:
            reasons.append("did_not_meet_locked_mutation_or_rppa_rule")
        if forced:
            reasons.append("M0_forced_exclusion_known_target")
        rows.append({
            "gene": gene, "genomic_driver": g, "recurrence": r, "lineage_biomarker": lineage,
            "core_pathway": 0.0, "prognostic_subtype": 0.0, "therapeutic_relevance": 0.0,
            "raw_weight": raw,
            "mutation_patients": int(mut[gene]["mutated_patients"]) if gene in mut else 0,
            "mutation_fraction": float(mut[gene]["fraction"]) if gene in mut else 0.0,
            "rppa_feature": rp["feature"] if rp else "",
            "rppa_abs_median_difference": float(rp["effect"]) if rp else 0.0,
            "rppa_direction": rp["direction"] if rp else "",
            "rppa_q_bh_198": float(rp["q"]) if rp else 1.0,
            "include_M0": "yes" if raw > 0 and not forced else "no",
            "include_M1": "yes" if raw > 0 else "no",
            "reasons": ";".join(reasons),
            "in_STRING_v12": "PENDING_LOCAL_CHECK",
        })
    return rows


def arm_weights(seed_rows: Sequence[Mapping[str, object]], arm: str, uniform: bool = False) -> list[dict[str, object]]:
    key = f"include_{arm}"
    kept = [r for r in seed_rows if r[key] == "yes"]
    if not kept:
        return []
    raw = [1.0 if uniform else float(r["raw_weight"]) for r in kept]
    total = math.fsum(raw)
    return [{"arm": arm, "gene": r["gene"], "raw_weight": float(r["raw_weight"]),
             "normalized_weight": w / total, "in_STRING_v12": r["in_STRING_v12"]}
            for r, w in zip(kept, raw, strict=True)]


# ------------------------------------------------------------- freeze gate
def check_stage_a_freeze(protocol_dir: Path) -> dict[str, str]:
    freeze = (protocol_dir / "FREEZE_STAGE_A.txt").read_text(encoding="utf-8")
    expected = {}
    for line in freeze.splitlines():
        parts = line.split("  ")
        if len(parts) == 2 and len(parts[0]) == 64:
            expected[parts[1].strip()] = parts[0]
    if set(expected) != {"PROJECT_PLAN.md", "protocol_config.json", "prefreeze_network_check.csv"}:
        raise RuntimeError("FREEZE_STAGE_A.txt does not list the three frozen files")
    for name, digest in expected.items():
        if sha256(protocol_dir / name) != digest:
            raise RuntimeError(f"frozen file changed since stage-A freeze: {name}")
    config = json.loads((protocol_dir / "protocol_config.json").read_text(encoding="utf-8"))
    if config.get("status") != "STAGE_A_FROZEN" or config.get("execution_allowed") != "stage_A_only":
        raise RuntimeError("config does not authorise stage-A execution")
    return expected


def check_amendments(protocol_dir: Path, required: Sequence[str]) -> dict[str, str]:
    """Verify confirmed stage-A amendments listed in FREEZE_STAGE_A_ADDENDUM.txt."""
    path = protocol_dir / "FREEZE_STAGE_A_ADDENDUM.txt"
    if not required:
        return {}
    if not path.exists():
        raise RuntimeError("required amendment not confirmed: FREEZE_STAGE_A_ADDENDUM.txt missing")
    listed = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split("  ")
        if len(parts) == 2 and len(parts[0]) == 64:
            listed[parts[1].strip()] = parts[0]
    for name in required:
        if name not in listed or sha256(protocol_dir / name) != listed[name]:
            raise RuntimeError(f"amendment missing or changed since confirmation: {name}")
    return {name: listed[name] for name in required}


REQUIRED_AMENDMENTS = ("AMENDMENT_A1_braf_classification.md",)


# --------------------------------------------------------------- branches
def run_branch(name: str, data: Path, inputs: Mapping[str, object], sample_type: str | None,
               patient_filter: str | None = None, allele_filter: str | None = None,
               multi_gene_mode: str = "single_gene_only", braf_mode: str = "A1") -> dict[str, object]:
    chosen, log = build_cohort(inputs["samples"], inputs["patients"], inputs["sequenced"],
                               inputs["rppa_samples"], inputs["braf"], sample_type, patient_filter, braf_mode)
    if allele_filter is not None:
        chosen = {p: v for p, v in chosen.items() if v["group"] == "REF" or v["allele"] == allele_filter}
    ek = {v["sample_id"]: p for p, v in chosen.items() if v["group"] == "EK"}
    ref = sorted(v["sample_id"] for v in chosen.values() if v["group"] == "REF")
    mutation_rows = mutation_statistics(data / "data_mutations.txt", ek)
    rppa_rows, rppa_selected = rppa_statistics(data / "data_rppa_zscores.txt", sorted(ek), ref, multi_gene_mode)
    seed_rows = build_seed_rows(mutation_rows, rppa_selected)
    counts = {
        "EK_patients": len(ek),
        "E_patients": sum(v["allele"] == "E" for v in chosen.values() if v["group"] == "EK"),
        "K_patients": sum(v["allele"] == "K" for v in chosen.values() if v["group"] == "EK"),
        "REF_patients": len(ref),
        "rppa_features_tested": sum(not r["not_tested_reason"] for r in rppa_rows),
        "rppa_features_selected": sum(r["selected"] == "yes" for r in rppa_rows),
        "mutation_genes_eligible": sum(r["eligible"] == "yes" for r in mutation_rows),
        "seeds_M0": sum(r["include_M0"] == "yes" for r in seed_rows),
        "seeds_M1": sum(r["include_M1"] == "yes" for r in seed_rows),
        "M1_restored": [r["gene"] for r in seed_rows if r["include_M1"] == "yes" and r["include_M0"] == "no"],
    }
    return {"name": name, "chosen": chosen, "log": log, "mutation": mutation_rows,
            "rppa": rppa_rows, "seeds": seed_rows, "counts": counts}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    frozen = check_stage_a_freeze(args.protocol)
    amendments = check_amendments(args.protocol, REQUIRED_AMENDMENTS)
    config = json.loads((args.protocol / "protocol_config.json").read_text(encoding="utf-8"))
    data = args.data
    input_files = ["data_clinical_patient.txt", "data_clinical_sample.txt", "data_mutations.txt",
                   "data_rppa_zscores.txt", "cases_sequenced.txt"]
    with (data / "data_rppa_zscores.txt").open(encoding="utf-8-sig") as handle:
        rppa_header = handle.readline().rstrip("\n").split("\t")
    inputs = {
        "samples": read_cbio_tsv(data / "data_clinical_sample.txt"),
        "patients": read_cbio_tsv(data / "data_clinical_patient.txt"),
        "sequenced": read_case_list(data / "cases_sequenced.txt"),
        "rppa_samples": set(rppa_header[1:]),
        "braf": braf_records_by_sample(data / "data_mutations.txt"),
    }
    expected = config["expected_feasibility_counts"]
    primary = run_branch("primary_metastasis", data, inputs, "Metastasis")
    c = primary["counts"]
    count_check = {"EK_expected": expected["V600E_K"], "EK_observed": c["EK_patients"],
                   "REF_expected": expected["no_reported_BRAF_variant"], "REF_observed": c["REF_patients"],
                   "E_expected": 92, "E_observed": c["E_patients"], "K_expected": 14, "K_observed": c["K_patients"]}
    count_check["pass"] = all(count_check[f"{k}_expected"] == count_check[f"{k}_observed"] for k in ("EK", "REF", "E", "K"))
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    manifest = {"analysis_version": ANALYSIS_VERSION, "protocol_id": PROTOCOL_ID,
                "stage": "A", "frozen_files_verified": frozen, "amendments_verified": amendments,
                "datahub_commit": config["datahub_commit"], "study": config["study"],
                "inputs_sha256": {f: sha256(data / f) for f in input_files},
                "code_sha256": sha256(Path(__file__)),
                "python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__,
                "platform": platform.platform(), "primary_count_check": count_check}
    if not count_check["pass"]:
        write_csv(out / "STOPPED_cohort_log_primary.csv", primary["log"])
        manifest["status"] = "STOPPED_COUNT_MISMATCH"
        (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
        print(json.dumps(count_check, indent=2))
        return 2
    branches = [primary,
                run_branch("S2_E_only", data, inputs, "Metastasis", allele_filter="E"),
                run_branch("S3_neoadjuvant_No", data, inputs, "Metastasis", patient_filter="neoadjuvant_no"),
                run_branch("S7_all_samples", data, inputs, None),
                run_branch("S9_multigene_split", data, inputs, "Metastasis", multi_gene_mode="split_equally"),
                run_branch("A1sens_strict_braf", data, inputs, "Metastasis", braf_mode="strict")]
    summary = {}
    for branch in branches:
        d = out / branch["name"]
        write_csv(d / "cohort_sample_log.csv", branch["log"])
        write_csv(d / "cohort_selected_patients.csv",
                  [{"patient_id": p, **v} for p, v in sorted(branch["chosen"].items())])
        write_csv(d / "mutation_statistics.csv", branch["mutation"])
        write_csv(d / "rppa_feature_tests.csv", branch["rppa"])
        write_csv(d / "seed_candidates.csv", branch["seeds"])
        weights = arm_weights(branch["seeds"], "M0") + arm_weights(branch["seeds"], "M1")
        if weights:
            write_csv(d / "seeds_frozen.csv", weights)
        summary[branch["name"]] = branch["counts"]
    task2_expected = {"S3_neoadjuvant_No": (100, 121), "S7_all_samples": (142, 152), "A1sens_strict_braf": (103, 125)}
    manifest["sensitivity_count_checks"] = {
        name: {"EK_expected": e, "REF_expected": r,
               "EK_observed": summary[name]["EK_patients"], "REF_observed": summary[name]["REF_patients"],
               "match": (summary[name]["EK_patients"], summary[name]["REF_patients"]) == (e, r)}
        for name, (e, r) in task2_expected.items()}
    s8 = arm_weights(primary["seeds"], "M0", uniform=True)
    write_csv(out / "S8_uniform_weights" / "seeds_frozen.csv", s8)
    summary["S8_uniform_weights"] = {"seeds_M0": len(s8), "note": "M0 membership of primary, weight 1/n"}
    summary["not_run_in_this_script"] = {
        "S4": "needs the two site-conflict sample IDs registered in task 2",
        "S5": "HC3 covariate-adjusted RPPA selection; separate implementation",
        "S6": "primary samples 36 vs 28: descriptive only, no reselection",
        "S1/S10": "use primary M0 seeds; stage B (propagation)"}
    manifest["status"] = "COMPLETED"
    manifest["branch_counts"] = summary
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
