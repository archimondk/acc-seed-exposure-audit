"""Locked primary-model pipeline for the ACC-PHARMA-NET revision.

This module is the single source of truth for the 108-drug complete-case
primary and evidence-informed rankings. It intentionally does not import any
legacy ADRS ranking, so old weights and hidden confidence transforms cannot
leak into the revised analysis.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import scipy
from scipy.stats import rankdata


@dataclass(frozen=True)
class AnalysisInputs:
    project_root: Path
    input_paths: Mapping[str, Path]
    associations: Mapping[str, frozenset[str]]
    edge_count: int
    r_acc: Mapping[str, float]
    mipe_mean_zauc: Mapping[str, float | None]
    nci60_potency: Mapping[str, float]
    external_score: Mapping[str, float]


@dataclass(frozen=True)
class AnalysisResult:
    primary_rows: tuple[dict[str, Any], ...]
    evidence_rows: tuple[dict[str, Any], ...]
    context_only_rows: tuple[dict[str, Any], ...]
    metrics: dict[str, Any]


def _read_csv(path: Path, required: Iterable[str]) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Required input is missing: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"Input has no data rows: {path}")
    missing = set(required) - set(rows[0])
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
    return rows


def load_inputs(project_root: Path) -> AnalysisInputs:
    root = project_root.resolve()
    network_dir = root / "data" / "bindex_network"
    input_paths = {
        "associations": network_dir / "bindex_edges_1304.csv",
        "r_acc": network_dir / "rACC_399_fullSTRING.csv",
        "mipe": network_dir / "Sactivity_124_v1.csv",
        "nci60": network_dir / "NCI60_potency_124.csv",
        "external": network_dir / "S_external_curated.csv",
    }

    edge_rows = _read_csv(input_paths["associations"], {"drug", "gene"})
    associations_mutable: dict[str, set[str]] = defaultdict(set)
    unique_edges: set[tuple[str, str]] = set()
    for row in edge_rows:
        drug = row["drug"].strip()
        gene = row["gene"].strip()
        if not drug or not gene:
            raise ValueError("Association rows must have non-empty drug and gene names")
        edge = (drug, gene)
        if edge in unique_edges:
            raise ValueError(f"Duplicate drug-gene association: {edge}")
        unique_edges.add(edge)
        associations_mutable[drug].add(gene)
    associations = {
        drug: frozenset(genes) for drug, genes in associations_mutable.items()
    }

    r_acc_rows = _read_csv(input_paths["r_acc"], {"gene", "rACC_full"})
    r_acc = {row["gene"].strip(): float(row["rACC_full"]) for row in r_acc_rows}
    associated_genes = set().union(*associations.values())
    missing_r_acc = associated_genes - set(r_acc)
    if missing_r_acc:
        raise ValueError(
            f"ACC relevance is missing for {len(missing_r_acc)} associated genes"
        )

    mipe_rows = _read_csv(input_paths["mipe"], {"drug", "mean_ZAUC"})
    mipe_mean_zauc = {
        row["drug"].strip(): (
            float(row["mean_ZAUC"]) if row["mean_ZAUC"].strip() else None
        )
        for row in mipe_rows
    }

    nci_rows = _read_csv(
        input_paths["nci60"], {"drug", "NCI60_mean_neglogGI50"}
    )
    nci60_potency = {
        row["drug"].strip(): float(row["NCI60_mean_neglogGI50"])
        for row in nci_rows
    }

    external_rows = _read_csv(input_paths["external"], {"drug", "S_external"})
    external_score = {drug: 0.5 for drug in associations}
    for row in external_rows:
        drug = row["drug"].strip()
        if drug in external_score:
            score = float(row["S_external"])
            if not 0.0 <= score <= 1.0:
                raise ValueError(f"External score outside [0, 1] for {drug}: {score}")
            external_score[drug] = score

    drugs = set(associations)
    missing_mipe_rows = drugs - set(mipe_mean_zauc)
    missing_nci_rows = drugs - set(nci60_potency)
    if missing_mipe_rows:
        raise ValueError(f"MIPE table is missing drugs: {sorted(missing_mipe_rows)}")
    if missing_nci_rows:
        raise ValueError(f"NCI-60 table is missing drugs: {sorted(missing_nci_rows)}")

    return AnalysisInputs(
        project_root=root,
        input_paths=input_paths,
        associations=associations,
        edge_count=len(unique_edges),
        r_acc=r_acc,
        mipe_mean_zauc=mipe_mean_zauc,
        nci60_potency=nci60_potency,
        external_score=external_score,
    )


def percentile_average(values: Mapping[str, float]) -> dict[str, float]:
    """Return average-rank percentiles scaled to [0, 1].

    The drug names are sorted before ranking so the returned mapping is fully
    deterministic. Tied values receive their average rank.
    """

    keys = sorted(values)
    if len(keys) < 2:
        raise ValueError("At least two values are required for percentile ranking")
    numeric = np.asarray([values[key] for key in keys], dtype=float)
    if not np.isfinite(numeric).all():
        raise ValueError("Percentile inputs must all be finite")
    percentiles = (rankdata(numeric, method="average") - 1.0) / (len(keys) - 1.0)
    return {key: float(value) for key, value in zip(keys, percentiles, strict=True)}


def _ordinal_ranks(scores: Mapping[str, float]) -> dict[str, int]:
    ordered = sorted(scores, key=lambda drug: (-scores[drug], drug))
    return {drug: index for index, drug in enumerate(ordered, start=1)}


def compute_primary_analysis(inputs: AnalysisInputs) -> AnalysisResult:
    drugs_all = sorted(inputs.associations)
    associated_genes = set().union(*inputs.associations.values())

    association_weighted_mean = float(
        np.mean(
            [
                inputs.r_acc[gene]
                for drug in drugs_all
                for gene in inputs.associations[drug]
            ]
        )
    )
    pseudo_count = 3.0
    c_acc: dict[str, float] = {}
    for drug in drugs_all:
        genes = inputs.associations[drug]
        n_assoc = len(genes)
        mean_drug = float(np.mean([inputs.r_acc[gene] for gene in genes]))
        c_acc[drug] = (
            n_assoc * mean_drug + pseudo_count * association_weighted_mean
        ) / (n_assoc + pseudo_count)

    primary_drugs = sorted(
        drug
        for drug in drugs_all
        if inputs.mipe_mean_zauc[drug] is not None
        and drug in inputs.nci60_potency
    )
    context_only_drugs = sorted(set(drugs_all) - set(primary_drugs))

    c_acc_pct = percentile_average({drug: c_acc[drug] for drug in primary_drugs})
    acc_potency_pct = percentile_average(
        {
            drug: -float(inputs.mipe_mean_zauc[drug])
            for drug in primary_drugs
            if inputs.mipe_mean_zauc[drug] is not None
        }
    )
    nci60_pct = percentile_average(
        {drug: inputs.nci60_potency[drug] for drug in primary_drugs}
    )

    x = np.asarray([nci60_pct[drug] for drug in primary_drugs], dtype=float)
    y = np.asarray([acc_potency_pct[drug] for drug in primary_drugs], dtype=float)
    slope = float(np.cov(x, y, bias=True)[0, 1] / np.var(x))
    intercept = float(y.mean() - slope * x.mean())
    residual = {
        drug: acc_potency_pct[drug] - (intercept + slope * nci60_pct[drug])
        for drug in primary_drugs
    }
    residual_pct = percentile_average(residual)

    primary_score = {
        drug: 0.5 * c_acc_pct[drug] + 0.5 * residual_pct[drug]
        for drug in primary_drugs
    }
    evidence_score = {
        drug: (
            0.4 * c_acc_pct[drug]
            + 0.4 * residual_pct[drug]
            + 0.2 * inputs.external_score[drug]
        )
        for drug in primary_drugs
    }
    primary_rank = _ordinal_ranks(primary_score)
    evidence_rank = _ordinal_ranks(evidence_score)

    primary_rows = tuple(
        {
            "rank_comp": primary_rank[drug],
            "drug": drug,
            "n_assoc": len(inputs.associations[drug]),
            "C_ACC": c_acc[drug],
            "C_ACC_pct": c_acc_pct[drug],
            "MIPE_potency_pct": acc_potency_pct[drug],
            "NCI60_potency_pct": nci60_pct[drug],
            "residual": residual[drug],
            "residual_pct": residual_pct[drug],
            "ADRS_comp": primary_score[drug],
        }
        for drug in sorted(primary_drugs, key=lambda item: primary_rank[item])
    )
    evidence_rows = tuple(
        {
            "rank_evidence_informed": evidence_rank[drug],
            "rank_comp": primary_rank[drug],
            "drug": drug,
            "C_ACC_pct": c_acc_pct[drug],
            "residual_pct": residual_pct[drug],
            "S_external": inputs.external_score[drug],
            "ADRS_EI": evidence_score[drug],
        }
        for drug in sorted(primary_drugs, key=lambda item: evidence_rank[item])
    )

    c_acc_pct_all = percentile_average(c_acc)
    c_acc_rank_all = _ordinal_ranks(c_acc)
    context_only_rows = tuple(
        {
            "drug": drug,
            "n_assoc": len(inputs.associations[drug]),
            "C_ACC": c_acc[drug],
            "C_ACC_pct_all124": c_acc_pct_all[drug],
            "rank_C_ACC_all124": c_acc_rank_all[drug],
            "missing_reason": "MIPE activity unavailable",
        }
        for drug in sorted(context_only_drugs, key=lambda item: c_acc_rank_all[item])
    )

    metrics: dict[str, Any] = {
        "model_version": "primary-108-v2",
        "n_edges": inputs.edge_count,
        "n_genes": len(associated_genes),
        "n_drugs_all": len(drugs_all),
        "n_primary": len(primary_drugs),
        "n_context_only": len(context_only_drugs),
        "c_acc_pseudo_count": pseudo_count,
        "c_acc_background_mean": association_weighted_mean,
        "percentile_method": "average rank scaled to [0,1] within analysis universe",
        "ols_slope": slope,
        "ols_intercept": intercept,
        "primary_formula": "0.50*C_ACC_pct + 0.50*residual_pct",
        "evidence_informed_formula": (
            "0.40*C_ACC_pct + 0.40*residual_pct + 0.20*S_external"
        ),
    }
    return AnalysisResult(
        primary_rows=primary_rows,
        evidence_rows=evidence_rows,
        context_only_rows=context_only_rows,
        metrics=metrics,
    )


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write an empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            formatted = {
                key: f"{value:.12g}" if isinstance(value, float) else value
                for key, value in row.items()
            }
            writer.writerow(formatted)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_outputs(
    inputs: AnalysisInputs,
    result: AnalysisResult,
    output_dir: Path,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "primary": output_dir / "ADRS_comp_primary_108.csv",
        "evidence": output_dir / "ADRS_evidence_informed_108.csv",
        "context_only": output_dir / "ADRS_context_only_16.csv",
        "metrics": output_dir / "primary_metrics.json",
        "manifest": output_dir / "run_manifest.md",
    }
    _write_csv(paths["primary"], result.primary_rows)
    _write_csv(paths["evidence"], result.evidence_rows)
    _write_csv(paths["context_only"], result.context_only_rows)
    with paths["metrics"].open("w", encoding="utf-8") as stream:
        json.dump(result.metrics, stream, indent=2, ensure_ascii=False)
        stream.write("\n")

    input_lines = "\n".join(
        f"- `{path.relative_to(inputs.project_root).as_posix()}`: `{_sha256(path)}`"
        for path in inputs.input_paths.values()
    )
    output_lines = "\n".join(
        f"- `{path.name}`: `{_sha256(path)}`"
        for key, path in paths.items()
        if key != "manifest"
    )
    manifest = f"""# Primary analysis run manifest

## Locked model

- Model version: `{result.metrics['model_version']}`
- Primary universe: `{result.metrics['n_primary']}` complete-case drugs
- C_ACC shrinkage pseudo-count: `{result.metrics['c_acc_pseudo_count']}`
- Percentile method: {result.metrics['percentile_method']}
- Primary formula: `{result.metrics['primary_formula']}`
- Evidence-informed formula: `{result.metrics['evidence_informed_formula']}`
- Clinical efficacy benchmark: retired after the C4 evidence-label audit because
  the strict subset contains two positive and zero negative clinical comparators.

## Environment

- Python: `{platform.python_version()}`
- NumPy: `{np.__version__}`
- SciPy: `{scipy.__version__}`
- Platform: `{platform.platform()}`
- Command: `python -m analysis.acc_primary_pipeline`

## Input SHA-256

{input_lines}

## Output SHA-256

{output_lines}

## Interpretation guardrail

No ROC-AUC or PR-AUC is emitted. The heterogeneous legacy labels are retained
only in the C4 evidence audit and are not a clinical-efficacy validation set.
"""
    paths["manifest"].write_text(manifest, encoding="utf-8")
    return paths


def run_pipeline(project_root: Path, output_dir: Path | None = None) -> AnalysisResult:
    inputs = load_inputs(project_root)
    result = compute_primary_analysis(inputs)
    destination = (
        output_dir
        if output_dir is not None
        else inputs.project_root / "results" / "primary_analysis"
    )
    write_outputs(inputs, result, destination)
    return result


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the locked 108-drug ACC-PHARMA-NET primary analysis."
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Project root containing data/bindex_network.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory; defaults to results/primary_analysis.",
    )
    return parser


def main() -> int:
    args = _build_parser().parse_args()
    result = run_pipeline(args.project_root, args.output_dir)
    print(
        json.dumps(
            {
                "status": "ok",
                "model_version": result.metrics["model_version"],
                "n_primary": result.metrics["n_primary"],
                "clinical_benchmark": "retired_in_C4",
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
