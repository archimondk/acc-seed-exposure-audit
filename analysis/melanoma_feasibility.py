"""Outcome-blind data feasibility audit; no propagation, drug ranking or efficacy tests."""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "melanoma_feasibility_20261018"
OUT = ROOT / "results" / "melanoma_feasibility_20261018"
COMMIT = "58341090c8bf0368ebe03f7aa95ec5137a8def25"
BASE = f"https://raw.githubusercontent.com/cBioPortal/datahub/{COMMIT}/public/skcm_tcga/"
MEDIA = f"https://media.githubusercontent.com/media/cBioPortal/datahub/{COMMIT}/public/skcm_tcga/"


def digest(path: Path, algorithm: str = "sha256") -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(2**20), b""):
            h.update(chunk)
    return h.hexdigest()


def save_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")


def download(name: str, url: str, expected: str | None = None, algorithm: str = "sha256") -> dict:
    path = RAW / name
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        print(f"Fetching {name}", flush=True)
        # Keep partial downloads separate; never replace an existing frozen input.
        part = path.with_suffix(path.suffix + ".part")
        with requests.get(url, timeout=(20, 90), stream=True) as response:
            response.raise_for_status()
            with part.open("wb") as stream:
                for chunk in response.iter_content(2**20):
                    stream.write(chunk)
        if expected and digest(part, algorithm) != expected:
            raise ValueError(f"Checksum mismatch: {name}")
        os.replace(part, path)
    if expected and digest(path, algorithm) != expected:
        raise ValueError(f"Cached checksum mismatch: {name}")
    record = {"file": name, "url": url, "bytes": path.stat().st_size,
              "sha256": digest(path), "verified_provider_hash": expected,
              "provider_hash_algorithm": algorithm if expected else None}
    print(f"Verified {name}: {record['bytes']} bytes", flush=True)
    return record


def fetch_inputs() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    names = ["data_clinical_sample.txt", "data_clinical_patient.txt", "data_mutations.txt",
             "data_rppa_zscores.txt", "data_gene_panel_matrix.txt", "data_mutsig.txt",
             "case_lists/cases_sequenced.txt", "meta_study.txt", "meta_mutations.txt",
             "meta_rppa_zscores.txt", "LICENSE"]

    def fetch_tcga(name: str) -> dict:
        response = requests.get(BASE + name, timeout=40)
        response.raise_for_status()
        pointer = re.search(r"oid sha256:([a-f0-9]{64})", response.text)
        if pointer:
            return download("tcga/" + name, MEDIA + name, pointer.group(1))
        return download("tcga/" + name, BASE + name)

    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        records = list(pool.map(fetch_tcga, names))
    meta_path = RAW / "depmap24q4_article.json"
    if not meta_path.exists():
        response = requests.get("https://api.figshare.com/v2/articles/27993248", timeout=40)
        response.raise_for_status()
        save_json(meta_path, response.json())
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    wanted = {"Model.csv", "OmicsSomaticMutations.csv", "OmicsProfiles.csv",
              "OmicsDefaultModelProfiles.csv", "OmicsSomaticMutationsMatrixHotspot.csv", "README.txt"}
    files = [f for f in meta["files"] if f["name"] in wanted]
    if {f["name"] for f in files} != wanted:
        raise ValueError("Required DepMap files missing from article")
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        records += list(pool.map(lambda f: download("depmap24q4/" + f["name"], f["download_url"],
                                                    f["computed_md5"], "md5"), files))
    save_json(RAW / "source_manifest.json", {
        "retrieved_utc": datetime.now(timezone.utc).isoformat(), "datahub_commit": COMMIT,
        "depmap_article_id": meta["id"], "depmap_article_version": meta["version"],
        "depmap_license": meta.get("license"), "inputs": records,
        "scope": "Local research audit only; no raw-data redistribution"})


def fetch_pancan() -> None:
    """Version-aligned alternative to legacy split-codon annotations."""
    study = "skcm_tcga_pan_can_atlas_2018"
    base = f"https://raw.githubusercontent.com/cBioPortal/datahub/{COMMIT}/public/{study}/"
    media = f"https://media.githubusercontent.com/media/cBioPortal/datahub/{COMMIT}/public/{study}/"
    names = ["data_clinical_sample.txt", "data_clinical_patient.txt", "data_mutations.txt",
             "data_rppa_zscores.txt", "case_lists/cases_sequenced.txt", "meta_study.txt",
             "meta_mutations.txt", "meta_rppa_zscores.txt", "LICENSE"]

    def fetch_one(name: str) -> dict:
        response = requests.get(base + name, timeout=40)
        response.raise_for_status()
        pointer = re.search(r"oid sha256:([a-f0-9]{64})", response.text)
        return download("pancan/" + name, (media if pointer else base) + name,
                        pointer.group(1) if pointer else None)

    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        records = list(pool.map(fetch_one, names))
    save_json(RAW / "pancan_source_manifest.json", {
        "retrieved_utc": datetime.now(timezone.utc).isoformat(), "datahub_commit": COMMIT,
        "study": study, "inputs": records,
        "reason": "Audit legacy split-codon annotations; no drug outcomes inspected"})


def verify_prism_source() -> None:
    path = RAW / "prism20q2_article.json"
    if not path.exists():
        r = requests.get("https://api.figshare.com/v2/articles/20564034", timeout=40)
        r.raise_for_status()
        save_json(path, r.json())
    meta = json.loads(path.read_text(encoding="utf-8"))
    name = "prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv"
    remote = next(f for f in meta["files"] if f["name"] == name)
    local = ROOT / name
    assert digest(local, "md5") == remote["computed_md5"], "Local PRISM differs from official file"
    save_json(RAW / "prism_source_manifest.json", {
        "checked_utc": datetime.now(timezone.utc).isoformat(), "article_id": meta["id"],
        "article_version": meta["version"], "license": meta["license"], "file": name,
        "url": remote["download_url"], "bytes": local.stat().st_size,
        "provider_md5_verified": remote["computed_md5"], "sha256": digest(local)})


def classify_braf(changes: list[str], sequenced: bool) -> str:
    """Do not turn an unprofiled sample or split-codon event into a WT/E/K call."""
    if not sequenced:
        return "not_profiled"
    variants = {str(v).removeprefix("p.") for v in changes if pd.notna(v) and str(v)}
    if not changes:
        return "no_reported_BRAF_variant"
    if not variants:
        return "other_BRAF_variant"
    # V600E plus a synonymous V600= record is still an explicit E call;
    # a second nonsynonymous change (e.g. E+M) may encode a split MNV.
    codon600 = {v for v in variants if "600" in v and not v.endswith("=") and v != "V600V"}
    if len(codon600) > 1 or any("del" in v or "ins" in v for v in codon600):
        return "ambiguous_V600"
    if "V600E" in variants:
        return "V600E"
    if "V600K" in variants:
        return "V600K"
    return "other_BRAF_variant"


def normalize_name(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(value).upper())


def variant_text(changes: list[str]) -> str:
    return ";".join(sorted({str(v) if pd.notna(v) else "UNANNOTATED" for v in changes}))


def select_one_sample(frame: pd.DataFrame, scenario: str) -> pd.DataFrame:
    if scenario not in {"all", "primary", "metastatic"}:
        raise ValueError(scenario)
    eligible = frame.copy()
    if scenario != "all":
        sample_type = "Primary" if scenario == "primary" else "Metastasis"
        eligible = eligible[eligible.SAMPLE_TYPE.eq(sample_type)]
    # Scenario-specific deterministic selection, not an outcome/feature-driven choice.
    return eligible.sort_values(["PATIENT_ID", "SAMPLE_ID"]).drop_duplicates("PATIENT_ID")


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#", low_memory=False,
                       na_values=["[Not Available]", "[Not Applicable]", "[Unknown]", "[Not Evaluated]"])


def mutation_subset(path: Path, gene_col: str, columns: list[str]) -> tuple[pd.DataFrame, int, set[str]]:
    chunks, count, samples = [], 0, set()
    id_col = "Tumor_Sample_Barcode" if gene_col == "Hugo_Symbol" else "ModelID"
    for chunk in pd.read_csv(path, sep="\t" if gene_col == "Hugo_Symbol" else ",",
                             comment="#", usecols=columns, chunksize=100_000, low_memory=False):
        count += len(chunk)
        samples.update(chunk[id_col].dropna().astype(str))
        chunks.append(chunk.loc[chunk[gene_col].eq("BRAF")].copy())
    return pd.concat(chunks, ignore_index=True), count, samples


def audit_tcga(folder: str) -> tuple[pd.DataFrame, dict]:
    source = RAW / folder
    clinical = read_tsv(source / "data_clinical_sample.txt")
    patient = read_tsv(source / "data_clinical_patient.txt")
    rppa = read_tsv(source / "data_rppa_zscores.txt").set_index("Composite.Element.REF")
    assert not clinical.SAMPLE_ID.duplicated().any(), "Duplicated sample metadata"
    assert not patient.PATIENT_ID.duplicated().any(), "Duplicated patient metadata"
    assert not rppa.columns.duplicated().any(), "Duplicated RPPA sample columns"
    cases = (source / "case_lists/cases_sequenced.txt").read_text(encoding="utf-8")
    sequenced = set(re.search(r"case_list_ids:\s*(.*)", cases).group(1).strip().split("\t"))
    cols = ["Hugo_Symbol", "Tumor_Sample_Barcode", "HGVSp_Short", "HGVSc", "Chromosome",
            "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2",
            "Variant_Classification", "t_ref_count", "t_alt_count"]
    braf, mutation_rows, mutation_samples = mutation_subset(source / "data_mutations.txt", "Hugo_Symbol", cols)
    braf.to_csv(OUT / f"{folder}_BRAF_records.csv", index=False)
    grouped = braf.groupby("Tumor_Sample_Barcode").HGVSp_Short.apply(list).to_dict()
    clinical["sequenced"] = clinical.SAMPLE_ID.isin(sequenced)
    clinical["has_rppa"] = clinical.SAMPLE_ID.isin(rppa.columns)
    clinical["BRAF_changes"] = clinical.SAMPLE_ID.map(lambda s: variant_text(grouped.get(s, [])))
    clinical["BRAF_status"] = [classify_braf(grouped.get(s, []), s in sequenced) for s in clinical.SAMPLE_ID]
    clinical["same_sample_eligible"] = clinical.sequenced & clinical.has_rppa
    clinical["comparison_group"] = clinical.BRAF_status.replace({"V600E": "V600E_K", "V600K": "V600E_K"})
    clinical["analysis_candidate"] = clinical.same_sample_eligible & clinical.comparison_group.isin(
        ["V600E_K", "no_reported_BRAF_variant"])
    clinical["exclusion_reason"] = np.select(
        [~clinical.sequenced, ~clinical.has_rppa, clinical.BRAF_status.eq("ambiguous_V600"),
         clinical.BRAF_status.eq("other_BRAF_variant")],
        ["not_in_sequenced_list", "no_same_sample_RPPA", "ambiguous_codon600", "other_BRAF_variant"],
        default="eligible_candidate")
    clinical.to_csv(OUT / f"{folder}_sample_audit.csv", index=False)
    flow = [
        {"stage": "clinical_sample_records", "n": len(clinical)},
        {"stage": "sequenced_clinical_samples", "n": int(clinical.sequenced.sum())},
        {"stage": "same_sample_sequenced_and_RPPA", "n": int(clinical.same_sample_eligible.sum())},
        {"stage": "exact_EK_or_no_reported_BRAF_variant", "n": int(clinical.analysis_candidate.sum())},
        {"stage": "one_sample_per_patient_all", "n": len(select_one_sample(clinical[clinical.analysis_candidate], "all"))},
    ]
    pd.DataFrame(flow).to_csv(OUT / f"{folder}_inclusion_flow.csv", index=False)
    rows, feature_rows, covariate_rows = [], [], []
    candidates = clinical.loc[clinical.analysis_candidate].copy()
    duplicate = clinical[clinical.PATIENT_ID.duplicated(keep=False)]
    duplicate.to_csv(OUT / f"{folder}_duplicate_patients.csv", index=False)
    conflicts = clinical[clinical.sequenced].groupby("PATIENT_ID").comparison_group.nunique()
    conflict_patients = set(conflicts[conflicts > 1].index)
    clinical[clinical.PATIENT_ID.isin(conflict_patients)].to_csv(OUT / f"{folder}_genotype_discordant_patients.csv", index=False)
    for scenario in ["all", "primary", "metastatic"]:
        selected = select_one_sample(candidates, scenario)
        assert selected.PATIENT_ID.is_unique
        selected.to_csv(OUT / f"{folder}_candidate_{scenario}.csv", index=False)
        for group, part in selected.groupby("comparison_group"):
            values = rppa[part.SAMPLE_ID].apply(pd.to_numeric, errors="coerce")
            finite = np.isfinite(values)
            rows.append({"source": folder, "scenario": scenario, "group": group,
                         "patients": part.PATIENT_ID.nunique(), "samples": len(part),
                         "rppa_features": len(rppa), "finite_fraction": float(finite.to_numpy().mean()),
                         "features_with_at_least_50": int((finite.sum(axis=1) >= 50).sum()),
                         "features_with_at_least_20": int((finite.sum(axis=1) >= 20).sum()),
                         "primary_n": int(part.SAMPLE_TYPE.eq("Primary").sum()),
                         "metastatic_n": int(part.SAMPLE_TYPE.eq("Metastasis").sum())})
            for feature, n in finite.sum(axis=1).items():
                feature_rows.append({"scenario": scenario, "group": group, "feature": feature,
                                     "finite_n": int(n), "total_n": len(part)})
            merged = part.merge(patient, on="PATIENT_ID", how="left", suffixes=("", "_patient"), validate="many_to_one")
            fields = [c for c in merged.columns if c in {"SEX", "AGE", "HISTORY_NEOADJUVANT_TRTYN",
                       "TISSUE_SOURCE_SITE", "TUMOR_TISSUE_SITE", "METASTATIC_SITE_PATIENT",
                       "IS_FFPE", "SAMPLE_TYPE", "TUMOR_SITE", "HISTORY_NEOADJUVANT_TX_TYPE"}]
            for field in fields:
                if field == "AGE":
                    series = pd.to_numeric(merged[field], errors="coerce")
                    for stat, value in {"missing": series.isna().sum(), "median": series.median()}.items():
                        covariate_rows.append({"scenario": scenario, "group": group, "field": field,
                                               "value": stat, "n_or_value": value})
                else:
                    for value, n in merged[field].fillna("MISSING").value_counts().items():
                        covariate_rows.append({"scenario": scenario, "group": group, "field": field,
                                               "value": str(value), "n_or_value": int(n)})
    pd.DataFrame(rows).to_csv(OUT / f"{folder}_cohort_summary.csv", index=False)
    pd.DataFrame(feature_rows).to_csv(OUT / f"{folder}_RPPA_feature_coverage.csv", index=False)
    pd.DataFrame(covariate_rows).to_csv(OUT / f"{folder}_covariate_coverage.csv", index=False)
    metastatic = select_one_sample(candidates, "metastatic").merge(
        patient[["PATIENT_ID", "HISTORY_NEOADJUVANT_TRTYN"]], on="PATIENT_ID", how="left", validate="one_to_one")
    treatment_rows = []
    for group, part in metastatic.groupby("comparison_group"):
        treatment_rows.append({"group": group, "metastatic_n": len(part),
                               "recorded_no_neoadjuvant_n": int(part.HISTORY_NEOADJUVANT_TRTYN.eq("No").sum()),
                               "recorded_yes_neoadjuvant_n": int(part.HISTORY_NEOADJUVANT_TRTYN.eq("Yes").sum()),
                               "missing_neoadjuvant_n": int(part.HISTORY_NEOADJUVANT_TRTYN.isna().sum())})
    pd.DataFrame(treatment_rows).to_csv(OUT / f"{folder}_treatment_scenario_counts.csv", index=False)
    features = pd.DataFrame({"feature": rppa.index,
                             "gene_annotation": [str(x).split("|")[0] for x in rppa.index]})
    features["multi_gene"] = features.gene_annotation.str.split().str.len().gt(1)
    features.to_csv(OUT / f"{folder}_RPPA_feature_dictionary.csv", index=False)
    summary = {"clinical_samples": len(clinical), "clinical_patients": clinical.PATIENT_ID.nunique(),
               "sequenced_list_n": len(sequenced), "mutation_rows": mutation_rows,
               "mutation_samples": len(mutation_samples), "rppa_samples": len(rppa.columns),
               "rppa_features": len(rppa), "multi_gene_features": int(features.multi_gene.sum()),
               "rppa_duplicate_feature_names": int(rppa.index.duplicated().sum()),
               "same_sample_seq_rppa": int(clinical.same_sample_eligible.sum()),
               "BRAF_status_all": clinical.BRAF_status.value_counts().to_dict(),
               "duplicate_patient_samples": len(duplicate),
               "genotype_discordant_patients": len(conflict_patients),
               "cohorts": rows}
    return clinical, summary


DRUG_ALIASES = {"dabrafenib": "Dabrafenib", "vemurafenib": "Vemurafenib",
                "encorafenib": "Encorafenib", "lgx818": "Encorafenib",
                "trametinib": "Trametinib", "cobimetinib": "Cobimetinib",
                "binimetinib": "Binimetinib", "mek162": "Binimetinib"}


def external_model_mask(frame: pd.DataFrame) -> pd.Series:
    # After a left merge booleans may be object/nullable; bitwise ~ on Python
    # object booleans yields integers. Unknown overlap must not pass the gate.
    return frame.QC_eligible.eq(True) & frame.NCI60_overlap.eq(False)


def audit_prism() -> dict:
    model = pd.read_csv(RAW / "depmap24q4/Model.csv")
    assert model.ModelID.is_unique
    profiles = pd.read_csv(RAW / "depmap24q4/OmicsDefaultModelProfiles.csv")
    dna_models = set(profiles.loc[profiles.ProfileType.str.lower().eq("dna"), "ModelID"])
    hotspot = pd.read_csv(RAW / "depmap24q4/OmicsSomaticMutationsMatrixHotspot.csv", usecols=[0])
    matrix_models = set(hotspot.iloc[:, 0])
    cols = ["HugoSymbol", "ModelID", "ProteinChange", "DNAChange", "Chrom", "Pos", "Ref", "Alt",
            "AF", "DP", "GT", "PS", "Hotspot", "VariantInfo"]
    braf, mutation_rows, mutation_models = mutation_subset(RAW / "depmap24q4/OmicsSomaticMutations.csv", "HugoSymbol", cols)
    braf.to_csv(OUT / "depmap_BRAF_records.csv", index=False)
    braf_map = braf.groupby("ModelID").ProteinChange.apply(list).to_dict()
    model["mutation_profile_available"] = model.ModelID.isin(dna_models & matrix_models)
    model["BRAF_changes"] = model.ModelID.map(lambda s: variant_text(braf_map.get(s, [])))
    model["BRAF_status"] = [classify_braf(braf_map.get(s, []), s in dna_models & matrix_models) for s in model.ModelID]
    model["melanoma_broad"] = model.OncotreePrimaryDisease.eq("Melanoma") & model.OncotreeLineage.eq("Skin")
    model["cutaneous_explicit"] = model.OncotreeSubtype.eq("Cutaneous Melanoma")
    # Resolve NCI-60 names to stable IDs, then expand to shared RRID/PatientID.
    nci_path = ROOT / "data/nci60_raw/output/documentation/NCI60_CELL_LINE_METADATA.xls"
    nci = pd.read_excel(nci_path, header=7)
    nci_names = nci["Cell Line Name"].dropna().astype(str)
    nci_names = nci_names[nci_names.str.match(r"^[A-Z]+:")]
    mapping_rows, direct_ids = [], set()
    manual_alias = {"MDAMB435": "MDAMB435S", "LOXIMVI": "LOXIMVI"}
    for original in nci_names:
        clean = re.sub(r"\s+[a-z]$", "", original.split(":", 1)[1].strip())
        key = normalize_name(clean)
        alias = manual_alias.get(key, key)
        matched = model[model.StrippedCellLineName.map(normalize_name).eq(alias)
                        | model.CellLineName.map(normalize_name).eq(alias)]
        direct_ids.update(matched.ModelID)
        mapping_rows.append({"NCI60_name": original, "normalized_name": key, "match_name": alias,
                             "ModelIDs": ";".join(matched.ModelID), "RRIDs": ";".join(matched.RRID.dropna()),
                             "PatientIDs": ";".join(matched.PatientID.dropna()),
                             "mapping": "explicit_alias" if key != alias else "normalized_exact",
                             "matched_n": len(matched)})
    pd.DataFrame(mapping_rows).to_csv(OUT / "NCI60_DepMap_identity_map.csv", index=False)
    direct = model[model.ModelID.isin(direct_ids)]
    rrids = set(direct.RRID.dropna())
    patients = set(direct.PatientID.dropna())
    model["NCI60_overlap"] = model.ModelID.isin(direct_ids) | model.RRID.isin(rrids) | model.PatientID.isin(patients)
    model.to_csv(OUT / "depmap_model_genotype_audit.csv", index=False)
    small, network = [], []
    network_drugs = pd.read_csv(ROOT / "data/bindex_network/bindex_drugs_124.csv")
    drug_col = "drug" if "drug" in network_drugs else "name"
    aliases = {normalize_name(n): n for n in network_drugs[drug_col]}
    aliases.update({normalize_name(k): v for k, v in DRUG_ALIASES.items()})
    use = ["broad_id", "depmap_id", "ccle_name", "screen_id", "passed_str_profiling", "name", "auc", "r2"]
    total = 0
    for chunk in pd.read_csv(ROOT / "prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv",
                             usecols=use, chunksize=100_000, low_memory=False):
        total += len(chunk)
        chunk["canonical_drug"] = chunk.name.map(lambda n: aliases.get(normalize_name(n)))
        matching = chunk[chunk.canonical_drug.notna()].copy()
        matching["AUC_finite"] = np.isfinite(pd.to_numeric(matching.auc, errors="coerce"))
        matching["STR_pass"] = matching.passed_str_profiling.astype(str).str.upper().eq("TRUE")
        # Retain availability flags only; no drug-response values/rank comparisons are exported.
        matching = matching.drop(columns=["auc", "r2"])
        network.append(matching)
    screen = pd.concat(network, ignore_index=True).merge(model, left_on="depmap_id", right_on="ModelID",
                                                        how="left", validate="many_to_one")
    melanoma = screen[screen.melanoma_broad.fillna(False)].copy()
    melanoma["QC_eligible"] = melanoma.STR_pass & melanoma.AUC_finite
    melanoma["external_candidate"] = external_model_mask(melanoma)
    assert not melanoma.loc[melanoma.external_candidate, "NCI60_overlap"].any()
    melanoma["six_drug"] = melanoma.canonical_drug.isin(set(DRUG_ALIASES.values()))
    keep = ["canonical_drug", "name", "broad_id", "screen_id", "depmap_id", "CellLineName", "RRID", "PatientID",
            "OncotreeSubtype", "BRAF_status", "BRAF_changes", "mutation_profile_available", "STR_pass", "AUC_finite",
            "NCI60_overlap", "external_candidate", "cutaneous_explicit", "six_drug", "CulturedResistanceDrug",
            "EngineeredModel", "ModelTreatment"]
    melanoma[keep].to_csv(OUT / "PRISM_network_drug_model_availability.csv", index=False)
    for drug, part in melanoma.groupby("canonical_drug"):
        ext = part[part.external_candidate]
        ek = ext[ext.BRAF_status.isin(["V600E", "V600K"])]
        ref = ext[ext.BRAF_status.eq("no_reported_BRAF_variant")]
        small.append({"drug": drug, "six_drug": drug in set(DRUG_ALIASES.values()),
                      "available_models": part.depmap_id.nunique(),
                      "QC_models": part.loc[part.QC_eligible, "depmap_id"].nunique(),
                      "NCI60_overlap_models": part.loc[part.NCI60_overlap, "depmap_id"].nunique(),
                      "external_models": ext.depmap_id.nunique(),
                      "external_patientIDs": ext.PatientID.nunique(),
                      "external_EK_models": ek.depmap_id.nunique(), "external_EK_patients": ek.PatientID.nunique(),
                      "external_E_models": ek.loc[ek.BRAF_status.eq("V600E"), "depmap_id"].nunique(),
                      "external_K_models": ek.loc[ek.BRAF_status.eq("V600K"), "depmap_id"].nunique(),
                      "external_no_reported_variant_models": ref.depmap_id.nunique(),
                      "external_explicit_cutaneous_EK": ek.loc[ek.cutaneous_explicit, "depmap_id"].nunique(),
                      "duplicate_model_rows": int(part.depmap_id.duplicated().sum()),
                      "compound_ids": ";".join(sorted(part.broad_id.unique()))})
    coverage = pd.DataFrame(small).sort_values("drug")
    coverage.to_csv(OUT / "PRISM_drug_coverage_summary.csv", index=False)
    melanoma[melanoma.six_drug & melanoma.NCI60_overlap][keep].to_csv(OUT / "PRISM_six_drug_excluded_overlap.csv", index=False)
    ext_models = melanoma.loc[melanoma.external_candidate].drop_duplicates("ModelID")
    ext_models[ext_models.PatientID.duplicated(keep=False)][keep].to_csv(OUT / "PRISM_related_external_models.csv", index=False)
    return {"screen_rows_total": total, "depmap_mutation_rows": mutation_rows,
            "official_model_rows": len(model), "local_model_sha256": digest(ROOT / "Model.csv"),
            "official_model_sha256": digest(RAW / "depmap24q4/Model.csv"),
            "mutation_models_with_any_call": len(mutation_models),
            "network_compounds_matched_any_lineage": int(screen.canonical_drug.nunique()),
            "network_compounds_in_melanoma": len(coverage),
            "network_compounds_external_EK_at_least_3": int((coverage.external_EK_models >= 3).sum()),
            "network_compounds_external_EK_at_least_5": int((coverage.external_EK_models >= 5).sum()),
            "melanoma_model_subtypes": model.loc[model.melanoma_broad, "OncotreeSubtype"].value_counts().to_dict(),
            "six_drug_coverage": coverage[coverage.six_drug].to_dict(orient="records"),
            "unmapped_NCI60_names": [r["NCI60_name"] for r in mapping_rows if r["matched_n"] == 0]}


def audit_all() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    verify_prism_source()
    legacy, old_summary = audit_tcga("tcga")
    pancan, new_summary = audit_tcga("pancan")
    cross = legacy[["SAMPLE_ID", "BRAF_status", "BRAF_changes"]].merge(
        pancan[["SAMPLE_ID", "BRAF_status", "BRAF_changes"]], on="SAMPLE_ID", how="outer",
        suffixes=("_legacy", "_pancan"), indicator=True)
    cross.to_csv(OUT / "BRAF_cross_source_annotation_audit.csv", index=False)
    prism = audit_prism()
    save_json(OUT / "feasibility_metrics.json", {"legacy": old_summary, "pancan": new_summary, "prism": prism})
    local_inputs = [ROOT / "Model.csv", ROOT / "prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv",
                    ROOT / "data/nci60_raw/output/documentation/NCI60_CELL_LINE_METADATA.xls",
                    ROOT / "data/bindex_network/bindex_drugs_124.csv"]
    save_json(OUT / "run_manifest.json", {"completed_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version, "pandas": pd.__version__, "platform": platform.platform(),
        "command": "python -m analysis.melanoma_feasibility --audit", "code_sha256": digest(Path(__file__)),
        "scope": "Availability and quality only. No seed selection, differential-expression tests, propagation or drug-ranking outcomes.",
        "source_manifests": {p.name: digest(p) for p in RAW.glob("*source_manifest.json")},
        "local_inputs": [{"path": str(p.relative_to(ROOT)), "sha256": digest(p)} for p in local_inputs],
        "outputs": [{"file": p.name, "sha256": digest(p)} for p in sorted(OUT.glob("*.csv"))]})
    print(json.dumps({"legacy_cohorts": old_summary["cohorts"], "pancan_cohorts": new_summary["cohorts"],
                      "six_drug_coverage": prism["six_drug_coverage"]}, indent=2), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument("--fetch-pancan", action="store_true")
    parser.add_argument("--audit", action="store_true")
    args = parser.parse_args()
    if args.fetch:
        fetch_inputs()
    if args.fetch_pancan:
        fetch_pancan()
    if args.audit:
        audit_all()


if __name__ == "__main__":
    main()
