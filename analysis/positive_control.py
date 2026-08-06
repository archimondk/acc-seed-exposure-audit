"""Outcome-blind ER+/HER2- breast-cancer positive control.

The disease seed is built without drug names, targets, screens, treatment
response, or observed ranks. The locked ACC graph, association network,
108-drug universe, propagation settings, and degree-matched null are reused.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import scipy
from scipy import sparse, stats

from analysis.acc_primary_pipeline import compute_primary_analysis, load_inputs
from analysis.mechanism_enrichment import benjamini_hochberg
from analysis.method_strengthening import (
    CDK46_DRUGS,
    C_ACC_PSEUDO_COUNT,
    NULL_DRAWS,
    REQUESTED_DEGREE_BINS,
    RNG_SEED,
    STRING_THRESHOLD,
    build_association_matrix,
    build_restart_matrix,
    compute_c_acc_matrix,
    empirical_upper_p,
    generate_degree_matched_seed_sets,
    load_string_graph,
    minmax_columns,
    percentile_columns,
    random_walk_with_restart,
)
from analysis.normalization_sensitivity import (
    VARIANT_DESCRIPTIONS,
    _gene_variants,
    network_smooth_with_restart,
    symmetric_normalized_operator,
)


ANALYSIS_VERSION = "positive-control-erpos-her2neg-v1"
DATAHUB_COMMIT = "58341090c8bf0368ebe03f7aa95ec5137a8def25"
CONTROL_TARGETS = frozenset({"CDK4", "CDK6"})
DIRECT_OUTCOME_TOKENS = (
    "drug",
    "treatment",
    "response",
    "screen",
    "rank",
    "target",
    "sensitivity",
    "viability",
)
PRIMARY_VARIANT = "column_minmax"
TCGA_MAIN_TEXT_MUTATION_GENES = frozenset(
    {
        "AFF2",
        "AKT1",
        "CBFB",
        "CCND3",
        "CDH1",
        "CDKN1B",
        "GATA3",
        "KMT2C",
        "MAP2K4",
        "MAP3K1",
        "NF1",
        "PIK3CA",
        "PIK3R1",
        "PTEN",
        "PTPN22",
        "PTPRD",
        "RB1",
        "RUNX1",
        "SF3B1",
        "TBX3",
        "TP53",
    }
)
NONSYNONYMOUS_TYPES = frozenset(
    {
        "Missense_Mutation",
        "Nonsense_Mutation",
        "Frame_Shift_Del",
        "Frame_Shift_Ins",
        "In_Frame_Del",
        "In_Frame_Ins",
        "Splice_Site",
        "Translation_Start_Site",
        "Nonstop_Mutation",
    }
)
SEED_COLUMNS = (
    "gene",
    "genomic_driver",
    "recurrence",
    "core_pathway",
    "lineage_biomarker",
    "prognostic_subtype",
    "raw_weight",
    "include_primary",
    "exclusion_reason",
    "source_id",
    "source_version",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def classify_her2(ihc_her2: Any, fish_her2: Any) -> str:
    ihc = str(ihc_her2).strip().lower()
    fish = str(fish_her2).strip().lower()
    if ihc == "positive" or fish == "positive":
        return "positive"
    if fish == "negative":
        return "negative"
    unavailable = {
        "",
        "nan",
        "[not available]",
        "[not evaluated]",
        "[not applicable]",
    }
    if fish in unavailable and ihc == "negative":
        return "negative"
    return "unknown"


def classify_erpos_her2neg(
    er_status: Any,
    ihc_her2: Any,
    fish_her2: Any,
) -> bool:
    return (
        str(er_status).strip().lower() == "positive"
        and classify_her2(ihc_her2, fish_her2) == "negative"
    )


def _read_cbio_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader((line for line in stream if not line.startswith("#")), delimiter="\t"))
    if not rows:
        raise ValueError(f"No data rows in {path}")
    return rows


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: f"{value:.12g}" if isinstance(value, float) else value
                    for key, value in row.items()
                }
            )


def _patient_classes(
    clinical_patient_path: Path,
) -> tuple[set[str], set[str], dict[str, str]]:
    target: set[str] = set()
    comparator: set[str] = set()
    classes: dict[str, str] = {}
    for row in _read_cbio_tsv(clinical_patient_path):
        patient = row["PATIENT_ID"].strip()
        er = row["ER_STATUS_BY_IHC"].strip().lower()
        her2 = classify_her2(row["IHC_HER2"], row["HER2_FISH_STATUS"])
        if er not in {"positive", "negative"} or her2 == "unknown":
            classes[patient] = "excluded_unknown"
        elif er == "positive" and her2 == "negative":
            target.add(patient)
            classes[patient] = "erpos_her2neg"
        else:
            comparator.add(patient)
            classes[patient] = "other_definitive"
    return target, comparator, classes


def _primary_samples_by_patient(
    clinical_sample_path: Path,
) -> dict[str, list[str]]:
    samples: dict[str, list[str]] = defaultdict(list)
    for row in _read_cbio_tsv(clinical_sample_path):
        if row["SAMPLE_TYPE"].strip().lower() == "primary":
            samples[row["PATIENT_ID"].strip()].append(row["SAMPLE_ID"].strip())
    return {patient: sorted(ids) for patient, ids in samples.items()}


def _rppa_statistics(
    rppa_path: Path,
    target_patients: set[str],
    comparator_patients: set[str],
    primary_samples: Mapping[str, Sequence[str]],
) -> tuple[list[dict[str, Any]], dict[str, dict[str, float]], dict[str, int]]:
    target_sample_to_patient = {
        sample: patient
        for patient in target_patients
        for sample in primary_samples.get(patient, ())
    }
    comparator_sample_to_patient = {
        sample: patient
        for patient in comparator_patients
        for sample in primary_samples.get(patient, ())
    }
    with rppa_path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        raw_rows = list(reader)
    sample_columns = header[1:]

    def choose_one(mapping: Mapping[str, str]) -> list[tuple[int, str]]:
        by_patient: dict[str, list[tuple[int, str]]] = defaultdict(list)
        for index, sample in enumerate(sample_columns, start=1):
            patient = mapping.get(sample)
            if patient:
                by_patient[patient].append((index, sample))
        return [
            sorted(entries, key=lambda item: item[1])[0]
            for entries in by_patient.values()
        ]

    target_columns = choose_one(target_sample_to_patient)
    comparator_columns = choose_one(comparator_sample_to_patient)
    if len(target_columns) < 50 or len(comparator_columns) < 50:
        raise ValueError("RPPA cohort has fewer than 50 samples in a group")

    feature_rows: list[dict[str, Any]] = []
    p_values: list[float] = []
    for raw in raw_rows:
        feature = raw[0]
        target_values = np.asarray(
            [float(raw[index]) for index, _ in target_columns if raw[index] != ""],
            dtype=float,
        )
        comparator_values = np.asarray(
            [
                float(raw[index])
                for index, _ in comparator_columns
                if raw[index] != ""
            ],
            dtype=float,
        )
        if len(target_values) < 50 or len(comparator_values) < 50:
            p_value = 1.0
            difference = math.nan
        else:
            difference = float(
                np.median(target_values) - np.median(comparator_values)
            )
            p_value = float(
                stats.mannwhitneyu(
                    target_values,
                    comparator_values,
                    alternative="two-sided",
                ).pvalue
            )
        genes = tuple(
            gene.strip().upper()
            for gene in feature.split("|", 1)[0].split()
            if gene.strip()
        )
        feature_rows.append(
            {
                "feature": feature,
                "genes": ";".join(genes),
                "target_n": len(target_values),
                "comparator_n": len(comparator_values),
                "median_difference": difference,
                "p_two_sided": p_value,
            }
        )
        p_values.append(p_value)
    q_values = benjamini_hochberg(p_values)
    selected_by_gene: dict[str, dict[str, float]] = {}
    for row, q_value in zip(feature_rows, q_values, strict=True):
        row["q_bh_225"] = float(q_value)
        selected = (
            int(row["target_n"]) >= 50
            and int(row["comparator_n"]) >= 50
            and float(q_value) < 0.01
            and abs(float(row["median_difference"])) >= 0.5
        )
        row["selected"] = "yes" if selected else "no"
        if selected:
            for gene in str(row["genes"]).split(";"):
                candidate = {
                    "effect": abs(float(row["median_difference"])),
                    "q": float(q_value),
                }
                current = selected_by_gene.get(gene)
                if current is None or candidate["effect"] > current["effect"]:
                    selected_by_gene[gene] = candidate
    counts = {
        "target_rppa_patients": len(target_columns),
        "comparator_rppa_patients": len(comparator_columns),
        "rppa_features_tested": len(feature_rows),
        "rppa_features_selected": sum(
            row["selected"] == "yes" for row in feature_rows
        ),
        "rppa_genes_selected": len(selected_by_gene),
    }
    return feature_rows, selected_by_gene, counts


def _mutation_statistics(
    mutation_json_path: Path,
    sequenced_sample_ids_path: Path,
    target_patients: set[str],
    primary_samples: Mapping[str, Sequence[str]],
) -> tuple[dict[str, dict[str, float]], dict[str, int]]:
    sequenced_samples = set(
        json.loads(sequenced_sample_ids_path.read_text(encoding="utf-8"))
    )
    sequenced_target_patients = {
        patient
        for patient in target_patients
        if any(
            sample in sequenced_samples
            for sample in primary_samples.get(patient, ())
        )
    }
    mutations = json.loads(mutation_json_path.read_text(encoding="utf-8"))
    mutated_patients: dict[str, set[str]] = defaultdict(set)
    for row in mutations:
        patient = str(row["patientId"])
        if patient not in sequenced_target_patients:
            continue
        if str(row.get("mutationType", "")) not in NONSYNONYMOUS_TYPES:
            continue
        gene = str(row["gene"]["hugoGeneSymbol"]).upper()
        if gene in TCGA_MAIN_TEXT_MUTATION_GENES:
            mutated_patients[gene].add(patient)
    denominator = len(sequenced_target_patients)
    if denominator < 50:
        raise ValueError("Too few sequenced positive-control patients")
    selected: dict[str, dict[str, float]] = {}
    for gene in sorted(TCGA_MAIN_TEXT_MUTATION_GENES):
        count = len(mutated_patients.get(gene, set()))
        prevalence = count / denominator
        if count >= 10 and prevalence >= 0.05:
            selected[gene] = {
                "mutated_patients": count,
                "prevalence": prevalence,
            }
    counts = {
        "target_sequenced_patients": denominator,
        "mutation_candidates": len(TCGA_MAIN_TEXT_MUTATION_GENES),
        "mutation_genes_selected": len(selected),
    }
    return selected, counts


def build_seed_table(
    data_dir: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    target, comparator, patient_classes = _patient_classes(
        data_dir / "data_clinical_patient.txt"
    )
    primary = _primary_samples_by_patient(
        data_dir / "data_clinical_sample.txt"
    )
    rppa_rows, rppa_selected, rppa_counts = _rppa_statistics(
        data_dir / "data_rppa_zscores.txt",
        target,
        comparator,
        primary,
    )
    mutation_selected, mutation_counts = _mutation_statistics(
        data_dir / "mutations_api.json",
        data_dir / "sequenced_sample_ids.json",
        target,
        primary,
    )
    candidates = sorted(
        TCGA_MAIN_TEXT_MUTATION_GENES | set(rppa_selected) | CONTROL_TARGETS
    )
    rows: list[dict[str, Any]] = []
    for gene in candidates:
        mutation = mutation_selected.get(gene)
        rppa = rppa_selected.get(gene)
        genomic_driver = 1.0 if mutation else 0.0
        recurrence = (
            min(float(mutation["prevalence"]) / 0.20, 1.0)
            if mutation
            else 0.0
        )
        lineage = (
            min(float(rppa["effect"]) / 2.0, 1.0) if rppa else 0.0
        )
        raw_weight = 0.30 * genomic_driver + 0.20 * recurrence + 0.10 * lineage
        include = raw_weight > 0 and gene not in CONTROL_TARGETS
        reasons: list[str] = []
        if gene in CONTROL_TARGETS:
            reasons.append("forced_exclusion_direct_positive_control_target")
        if raw_weight <= 0:
            reasons.append("did_not_meet_locked_mutation_or_rppa_threshold")
        sources: list[str] = []
        if mutation:
            sources.append("TCGA2012_main_text_SMG+brca_tcga_mutations")
        if rppa:
            sources.append("brca_tcga_RPPA_locked_contrast")
        rows.append(
            {
                "gene": gene,
                "genomic_driver": genomic_driver,
                "recurrence": recurrence,
                "core_pathway": 0.0,
                "lineage_biomarker": lineage,
                "prognostic_subtype": 0.0,
                "raw_weight": raw_weight,
                "include_primary": "yes" if include else "no",
                "exclusion_reason": ";".join(reasons),
                "source_id": ";".join(sources) or "candidate_audit_only",
                "source_version": (
                    f"PMID23000897;DataHub:{DATAHUB_COMMIT}"
                ),
                "mutation_prevalence": (
                    float(mutation["prevalence"]) if mutation else 0.0
                ),
                "mutation_patient_n": (
                    int(mutation["mutated_patients"]) if mutation else 0
                ),
                "rppa_abs_median_difference": (
                    float(rppa["effect"]) if rppa else 0.0
                ),
                "rppa_q_bh_225": float(rppa["q"]) if rppa else 1.0,
            }
        )
    metrics = {
        "analysis_version": ANALYSIS_VERSION,
        "datahub_commit": DATAHUB_COMMIT,
        "patient_total": len(patient_classes),
        "target_patients": len(target),
        "comparator_patients": len(comparator),
        "excluded_unknown_patients": sum(
            value == "excluded_unknown" for value in patient_classes.values()
        ),
        **rppa_counts,
        **mutation_counts,
        "candidate_gene_rows": len(rows),
        "retained_seed_genes": sum(
            row["include_primary"] == "yes" for row in rows
        ),
        "forced_target_exclusions": sorted(CONTROL_TARGETS),
        "outcome_blind": True,
    }
    return rows, rppa_rows, metrics


def load_positive_control_seed_weights(path: Path) -> dict[str, float]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        fieldnames = tuple(reader.fieldnames or ())
        suspicious = [
            field
            for field in fieldnames
            if any(token in field.lower() for token in DIRECT_OUTCOME_TOKENS)
        ]
        if suspicious:
            raise ValueError(
                f"Seed table contains outcome-derived fields: {suspicious}"
            )
        missing = set(SEED_COLUMNS) - set(fieldnames)
        if missing:
            raise ValueError(f"Seed table missing columns: {sorted(missing)}")
        rows = list(reader)
    weights: dict[str, float] = {}
    for row in rows:
        if not _truthy(row["include_primary"]):
            continue
        gene = row["gene"].strip().upper()
        if gene in CONTROL_TARGETS:
            raise ValueError(
                f"Seed includes direct positive-control target: {gene}"
            )
        weight = float(row["raw_weight"])
        if not np.isfinite(weight) or weight <= 0:
            raise ValueError(f"Invalid retained seed weight for {gene}")
        weights[gene] = weight
    if not weights:
        raise ValueError("Positive-control seed set is empty")
    total = sum(weights.values())
    return {gene: weight / total for gene, weight in weights.items()}


def evaluate_success_criteria(
    *,
    primary_group_p: float,
    primary_group_q: float,
    top_quartile_flags: Mapping[str, bool],
    concordant_variants: int,
) -> dict[str, Any]:
    criteria = {
        "primary_group_p_lt_0_05": primary_group_p < 0.05,
        "primary_group_q_lt_0_05": primary_group_q < 0.05,
        "at_least_two_cdk46_drugs_top_quartile": (
            sum(bool(top_quartile_flags.get(drug, False)) for drug in CDK46_DRUGS)
            >= 2
        ),
        "direction_concordant_at_least_three_variants": (
            concordant_variants >= 3
        ),
    }
    passed = sum(criteria.values())
    status = "pass" if passed == 4 else ("partial_recovery" if passed >= 2 else "fail")
    return {"status": status, "criteria": criteria}


def freeze_seed_inputs(
    project_root: Path,
    data_dir: Path,
    output_dir: Path,
) -> tuple[Path, dict[str, Any]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    seed_rows, rppa_rows, metrics = build_seed_table(data_dir)
    seed_path = output_dir / "positive_control_seed_frozen.csv"
    rppa_path = output_dir / "rppa_feature_audit.csv"
    cohort_path = output_dir / "cohort_and_seed_quality.json"
    manifest_path = output_dir / "seed_freeze_manifest.json"
    _write_csv(seed_path, seed_rows)
    _write_csv(rppa_path, rppa_rows)
    cohort_path.write_text(
        json.dumps(metrics, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    input_names = (
        "data_clinical_patient.txt",
        "data_clinical_sample.txt",
        "data_rppa_zscores.txt",
        "mutations_api.json",
        "sequenced_sample_ids.json",
        "PMC3465532_fullText.xml",
    )
    manifest = {
        "freeze_timestamp_local": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "protocol": str(
            (project_root / "experiments" / "positive_control_protocol.md")
            .relative_to(project_root)
            .as_posix()
        ),
        "protocol_sha256": _sha256(
            project_root / "experiments" / "positive_control_protocol.md"
        ),
        "datahub_commit": DATAHUB_COMMIT,
        "inputs": {
            name: {
                "sha256": _sha256(data_dir / name),
                "bytes": (data_dir / name).stat().st_size,
            }
            for name in input_names
        },
        "seed_csv": {
            "path": seed_path.name,
            "sha256": _sha256(seed_path),
            "bytes": seed_path.stat().st_size,
        },
        "rppa_feature_audit": {
            "path": rppa_path.name,
            "sha256": _sha256(rppa_path),
            "bytes": rppa_path.stat().st_size,
        },
        "quality_metrics": metrics,
        "drug_ranking_inspected_before_freeze": False,
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return seed_path, manifest


def _drug_percentiles(
    association_matrix: sparse.csr_matrix,
    primary_indices: np.ndarray,
    gene_values: np.ndarray,
) -> np.ndarray:
    c_acc_all = compute_c_acc_matrix(
        association_matrix,
        gene_values,
        pseudo_count=C_ACC_PSEUDO_COUNT,
    )
    return percentile_columns(c_acc_all[primary_indices, :])


def run_positive_control_analysis(
    project_root: Path,
    seed_path: Path,
    n_null: int = NULL_DRAWS,
    null_batch_size: int = 64,
    seed_weights_override: Mapping[str, float] | None = None,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, np.ndarray],
    list[dict[str, Any]],
    dict[str, Any],
]:
    root = project_root.resolve()
    inputs = load_inputs(root)
    primary = compute_primary_analysis(inputs)
    primary_drugs = tuple(str(row["drug"]) for row in primary.primary_rows)
    all_drugs = tuple(sorted(inputs.associations))
    associated_genes = tuple(sorted(set().union(*inputs.associations.values())))
    seed_weights = (
        dict(seed_weights_override)
        if seed_weights_override is not None
        else load_positive_control_seed_weights(seed_path)
    )
    if not seed_weights:
        raise ValueError("Seed-weight override must not be empty")
    if any(
        not math.isfinite(float(weight)) or float(weight) <= 0
        for weight in seed_weights.values()
    ):
        raise ValueError("Seed-weight override must contain finite positive weights")
    graph = load_string_graph(
        root / "9606.protein.info.v12.0.txt.gz",
        root / "9606.protein.links.v12.0.txt.gz",
        required_nodes=set(associated_genes) | set(seed_weights),
        threshold=STRING_THRESHOLD,
    )
    symmetric_operator = symmetric_normalized_operator(graph.adjacency)
    gene_indices = np.asarray(
        [graph.node_index[gene] for gene in associated_genes], dtype=int
    )
    association_matrix = build_association_matrix(
        all_drugs, associated_genes, inputs.associations
    )
    all_drug_index = {drug: index for index, drug in enumerate(all_drugs)}
    primary_indices = np.asarray(
        [all_drug_index[drug] for drug in primary_drugs], dtype=int
    )

    observed_restart = build_restart_matrix(graph, (seed_weights,))
    observed_column, column_iterations, column_delta = random_walk_with_restart(
        graph, observed_restart
    )
    observed_symmetric, symmetric_iterations, symmetric_delta = (
        network_smooth_with_restart(symmetric_operator, observed_restart)
    )
    uniform_restart = np.full(
        (len(graph.node_names), 1),
        1.0 / len(graph.node_names),
        dtype=float,
    )
    uniform_column, uniform_iterations, uniform_delta = (
        random_walk_with_restart(graph, uniform_restart)
    )
    observed_gene_variants = _gene_variants(
        observed_column[gene_indices, :],
        observed_symmetric[gene_indices, :],
        uniform_column[gene_indices, :],
    )
    observed_drug_variants = {
        name: _drug_percentiles(
            association_matrix, primary_indices, values
        )[:, 0]
        for name, values in observed_gene_variants.items()
    }

    matched_rows, degree_bin_edges = generate_degree_matched_seed_sets(
        seed_weights,
        graph.node_names,
        graph.degree,
        n_draws=n_null,
        rng_seed=RNG_SEED,
        n_bins=REQUESTED_DEGREE_BINS,
    )
    rows_by_replicate: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in matched_rows:
        rows_by_replicate[int(row["replicate"])].append(row)
    null_matrices = {
        name: np.empty((len(primary_drugs), n_null), dtype=np.float32)
        for name in VARIANT_DESCRIPTIONS
    }
    null_iteration_ranges = {
        "column": [],
        "symmetric": [],
        "column_delta": [],
        "symmetric_delta": [],
    }
    for batch_start in range(0, n_null, null_batch_size):
        replicates = list(
            range(batch_start, min(batch_start + null_batch_size, n_null))
        )
        batch_weights = [
            {
                row["null_seed"]: float(row["weight"])
                for row in rows_by_replicate[replicate]
            }
            for replicate in replicates
        ]
        restart = build_restart_matrix(graph, batch_weights)
        column, column_iteration, column_final_delta = random_walk_with_restart(
            graph, restart
        )
        symmetric, symmetric_iteration, symmetric_final_delta = (
            network_smooth_with_restart(symmetric_operator, restart)
        )
        null_iteration_ranges["column"].append(column_iteration)
        null_iteration_ranges["symmetric"].append(symmetric_iteration)
        null_iteration_ranges["column_delta"].append(column_final_delta)
        null_iteration_ranges["symmetric_delta"].append(symmetric_final_delta)
        variants = _gene_variants(
            column[gene_indices, :],
            symmetric[gene_indices, :],
            uniform_column[gene_indices, :],
        )
        for name, values in variants.items():
            null_matrices[name][:, replicates] = _drug_percentiles(
                association_matrix, primary_indices, values
            ).astype(np.float32)

    primary_index = {drug: index for index, drug in enumerate(primary_drugs)}
    cdk_indices = np.asarray(
        [primary_index[drug] for drug in CDK46_DRUGS], dtype=int
    )
    summary_rows: list[dict[str, Any]] = []
    drug_rows: list[dict[str, Any]] = []
    group_null_rows: list[dict[str, Any]] = []
    for name, description in VARIANT_DESCRIPTIONS.items():
        observed = observed_drug_variants[name]
        null = null_matrices[name].astype(float)
        null_mean = null.mean(axis=1)
        null_sd = null.std(axis=1, ddof=1)
        p_values = np.asarray(
            [
                empirical_upper_p(observed[index], null[index])
                for index in range(len(primary_drugs))
            ],
            dtype=float,
        )
        q_values = np.asarray(
            benjamini_hochberg(p_values.tolist()), dtype=float
        )
        ranks = stats.rankdata(-observed, method="min").astype(int)
        observed_group = float(observed[cdk_indices].mean())
        null_group = null[cdk_indices, :].mean(axis=0)
        group_p = empirical_upper_p(observed_group, null_group)
        summary_rows.append(
            {
                "variant": name,
                "description": description,
                "CDK46_observed_mean_percentile": observed_group,
                "CDK46_null_mean_percentile": float(null_group.mean()),
                "CDK46_empirical_p_upper": group_p,
                "direction_observed_gt_null_mean": (
                    "yes" if observed_group > float(null_group.mean()) else "no"
                ),
                "n_drugs_q_lt_0_05": int(np.count_nonzero(q_values < 0.05)),
                "minimum_drug_q_bh_108": float(q_values.min()),
            }
        )
        for index, drug in enumerate(primary_drugs):
            drug_rows.append(
                {
                    "variant": name,
                    "drug": drug,
                    "rank_descending": int(ranks[index]),
                    "top_quartile": (
                        "yes"
                        if int(ranks[index])
                        <= math.ceil(len(primary_drugs) * 0.25)
                        else "no"
                    ),
                    "observed_C_ACC_percentile": float(observed[index]),
                    "null_mean_C_ACC_percentile": float(null_mean[index]),
                    "null_sd_C_ACC_percentile": float(null_sd[index]),
                    "z_degree_matched": (
                        float((observed[index] - null_mean[index]) / null_sd[index])
                        if null_sd[index] > 0
                        else math.nan
                    ),
                    "empirical_p_upper": float(p_values[index]),
                    "q_bh_108": float(q_values[index]),
                    "null_draws": n_null,
                }
            )
        for replicate, value in enumerate(null_group):
            group_null_rows.append(
                {
                    "variant": name,
                    "replicate": replicate,
                    "CDK46_null_mean_percentile": float(value),
                }
            )
    group_q = benjamini_hochberg(
        [float(row["CDK46_empirical_p_upper"]) for row in summary_rows]
    )
    for row, q_value in zip(summary_rows, group_q, strict=True):
        row["CDK46_q_bh_across_variants"] = float(q_value)

    primary_summary = next(
        row for row in summary_rows if row["variant"] == PRIMARY_VARIANT
    )
    primary_drug_rows = {
        str(row["drug"]): row
        for row in drug_rows
        if row["variant"] == PRIMARY_VARIANT
        and row["drug"] in CDK46_DRUGS
    }
    top_flags = {
        drug: primary_drug_rows[drug]["top_quartile"] == "yes"
        for drug in CDK46_DRUGS
    }
    concordant = sum(
        row["direction_observed_gt_null_mean"] == "yes"
        for row in summary_rows
    )
    decision = evaluate_success_criteria(
        primary_group_p=float(primary_summary["CDK46_empirical_p_upper"]),
        primary_group_q=float(primary_summary["CDK46_q_bh_across_variants"]),
        top_quartile_flags=top_flags,
        concordant_variants=concordant,
    )
    metrics = {
        "analysis_version": ANALYSIS_VERSION,
        "positive_control_status": decision["status"],
        "criteria": decision["criteria"],
        "primary_variant": PRIMARY_VARIANT,
        "null_draws": n_null,
        "null_rng_seed": RNG_SEED,
        "null_batch_size": null_batch_size,
        "empirical_p_minimum_resolution": 1.0 / (n_null + 1.0),
        "primary_universe_n": len(primary_drugs),
        "primary_drug_order": list(primary_drugs),
        "associated_gene_n": len(associated_genes),
        "positive_control_seed_n": len(seed_weights),
        "degree_bins_effective": len(degree_bin_edges) - 1,
        "CDK46_primary_top_quartile": top_flags,
        "direction_concordant_variants": concordant,
        "observed_iterations": {
            "column": column_iterations,
            "column_final_delta": column_delta,
            "uniform": uniform_iterations,
            "uniform_final_delta": uniform_delta,
            "symmetric": symmetric_iterations,
            "symmetric_final_delta": symmetric_delta,
        },
        "null_iterations": {
            "column_min": min(null_iteration_ranges["column"]),
            "column_max": max(null_iteration_ranges["column"]),
            "column_final_delta_max": max(
                null_iteration_ranges["column_delta"]
            ),
            "symmetric_min": min(null_iteration_ranges["symmetric"]),
            "symmetric_max": max(null_iteration_ranges["symmetric"]),
            "symmetric_final_delta_max": max(
                null_iteration_ranges["symmetric_delta"]
            ),
        },
        "variants": {
            str(row["variant"]): {
                key: value
                for key, value in row.items()
                if key not in {"variant", "description"}
            }
            for row in summary_rows
        },
    }
    drug_rows.sort(key=lambda row: (row["variant"], row["rank_descending"], row["drug"]))
    return (
        summary_rows,
        drug_rows,
        group_null_rows,
        null_matrices,
        matched_rows,
        metrics,
    )


def write_analysis_outputs(
    project_root: Path,
    output_dir: Path,
    seed_path: Path,
    summary_rows: Sequence[Mapping[str, Any]],
    drug_rows: Sequence[Mapping[str, Any]],
    group_null_rows: Sequence[Mapping[str, Any]],
    null_matrices: Mapping[str, np.ndarray],
    matched_rows: Sequence[Mapping[str, Any]],
    metrics: Mapping[str, Any],
    wall_seconds: float,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "summary": output_dir / "positive_control_variant_summary.csv",
        "drugs": output_dir / "positive_control_drug_ranks_primary108.csv",
        "group_null": output_dir / "positive_control_CDK46_group_null_10000.csv",
        "drug_null": output_dir / "positive_control_drug_null_10000.npz",
        "matched_seeds": output_dir / "positive_control_degree_matched_seeds.csv",
        "metrics": output_dir / "positive_control_metrics.json",
        "report": output_dir / "positive_control_report.md",
        "manifest": output_dir / "run_manifest.json",
    }
    _write_csv(paths["summary"], summary_rows)
    _write_csv(paths["drugs"], drug_rows)
    _write_csv(paths["group_null"], group_null_rows)
    _write_csv(paths["matched_seeds"], matched_rows)
    np.savez_compressed(
        paths["drug_null"],
        drug_names=np.asarray(metrics["primary_drug_order"], dtype="U"),
        **{name: matrix for name, matrix in null_matrices.items()},
    )
    metrics_path_payload = dict(metrics)
    metrics_path_payload["wall_clock_seconds"] = wall_seconds
    paths["metrics"].write_text(
        json.dumps(metrics_path_payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    primary = next(
        row for row in summary_rows if row["variant"] == PRIMARY_VARIANT
    )
    cdk_primary = [
        row
        for row in drug_rows
        if row["variant"] == PRIMARY_VARIANT and row["drug"] in CDK46_DRUGS
    ]
    report = [
        "# ER+/HER2-negative breast-cancer positive-control report",
        "",
        f"- Prespecified decision: **{metrics['positive_control_status']}**.",
        f"- Frozen disease seeds: {metrics['positive_control_seed_n']}.",
        f"- Degree-matched null draws: {metrics['null_draws']:,}.",
        f"- Primary CDK4/6 group P: {primary['CDK46_empirical_p_upper']:.4g}.",
        f"- Primary CDK4/6 group q across variants: {primary['CDK46_q_bh_across_variants']:.4g}.",
        "",
        "## Primary CDK4/6 ranks",
        "",
        "| Drug | Rank / 108 | Top quartile | C_ACC percentile |",
        "|---|---:|:---:|---:|",
        *[
            f"| {row['drug']} | {row['rank_descending']} | "
            f"{row['top_quartile']} | {row['observed_C_ACC_percentile']:.3f} |"
            for row in cdk_primary
        ],
        "",
        "## Prespecified criteria",
        "",
        *[
            f"- {'PASS' if passed else 'FAIL'}: `{name}`"
            for name, passed in metrics["criteria"].items()
        ],
        "",
        "Recovery tests implementation transportability for one established "
        "disease-mechanism pair. It does not validate ACC candidates or "
        "establish efficacy.",
        "",
    ]
    paths["report"].write_text("\n".join(report), encoding="utf-8")
    input_paths = (
        seed_path,
        project_root / "data" / "bindex_network" / "bindex_edges_1304.csv",
        project_root / "data" / "bindex_network" / "rACC_399_fullSTRING.csv",
        project_root / "9606.protein.info.v12.0.txt.gz",
        project_root / "9606.protein.links.v12.0.txt.gz",
        project_root / "experiments" / "positive_control_protocol.md",
    )
    manifest = {
        "analysis_version": ANALYSIS_VERSION,
        "command": (
            "python -m analysis.positive_control --project-root . "
            "--n-null 10000"
        ),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "logical_cpu_count": os.cpu_count(),
        "wall_clock_seconds": wall_seconds,
        "rng_seed": metrics["null_rng_seed"],
        "null_draws": metrics["null_draws"],
        "inputs": {
            path.relative_to(project_root).as_posix(): _sha256(path)
            for path in input_paths
        },
        "outputs": {
            key: {
                "path": path.name,
                "sha256": _sha256(path),
                "bytes": path.stat().st_size,
            }
            for key, path in paths.items()
            if key != "manifest"
        },
    }
    paths["manifest"].write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return paths


def run(
    project_root: Path,
    data_dir: Path | None = None,
    output_dir: Path | None = None,
    n_null: int = NULL_DRAWS,
    null_batch_size: int = 64,
    rebuild_seed: bool = True,
) -> dict[str, Any]:
    root = project_root.resolve()
    data = (
        data_dir.resolve()
        if data_dir is not None
        else root / "data" / "positive_control" / "erpos_her2neg_tcga"
    )
    target = (
        output_dir.resolve()
        if output_dir is not None
        else root / "results" / "positive_control"
    )
    seed_path = target / "positive_control_seed_frozen.csv"
    if rebuild_seed:
        seed_path, _ = freeze_seed_inputs(root, data, target)
    elif not seed_path.is_file():
        raise FileNotFoundError(seed_path)
    started = time.perf_counter()
    (
        summary_rows,
        drug_rows,
        group_null_rows,
        null_matrices,
        matched_rows,
        metrics,
    ) = run_positive_control_analysis(
        root,
        seed_path,
        n_null=n_null,
        null_batch_size=null_batch_size,
    )
    wall_seconds = time.perf_counter() - started
    write_analysis_outputs(
        root,
        target,
        seed_path,
        summary_rows,
        drug_rows,
        group_null_rows,
        null_matrices,
        matched_rows,
        metrics,
        wall_seconds,
    )
    return metrics


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--n-null", type=int, default=NULL_DRAWS)
    parser.add_argument("--null-batch-size", type=int, default=64)
    parser.add_argument(
        "--reuse-frozen-seed",
        action="store_true",
        help="Use the seed CSV already frozen in the output directory.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    metrics = run(
        args.project_root,
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        n_null=args.n_null,
        null_batch_size=args.null_batch_size,
        rebuild_seed=not args.reuse_frozen_seed,
    )
    print(
        json.dumps(
            {
                "status": metrics["positive_control_status"],
                "criteria": metrics["criteria"],
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
