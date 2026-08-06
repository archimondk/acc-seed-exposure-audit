"""Quantify C_ACC pseudo-count sensitivity and MIPE-missingness patterns."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import scipy
from scipy import stats
from scipy.stats import rankdata

from analysis.acc_primary_pipeline import load_inputs


ANALYSIS_VERSION = "shrinkage-missingness-audit-v1"
PSEUDO_COUNTS = (1.0, 3.0, 5.0, 10.0)

MISSING_DRUG_CLASSES: dict[str, tuple[str, str]] = {
    "Pipobroman": ("Alkylator/platinum", "broad_cytotoxic"),
    "tepotinib": ("ALK/MET", "targeted"),
    "Arsenic trioxide": ("Pleiotropic/other", "other"),
    "Trilaciclib": ("CDK4/6", "targeted"),
    "Zanubrutinib": ("BTK", "targeted"),
    "Avapritinib": ("KIT/PDGFRA", "targeted"),
    "Talazoparib": ("PARP", "targeted"),
    "Selpercatinib": ("RET", "targeted"),
    "Fedratinib": ("JAK2", "targeted"),
    "Triethylenemelamine": ("Alkylator/platinum", "broad_cytotoxic"),
    "Sotorasib": ("KRAS G12C", "targeted"),
    "Uracil mustard": ("Alkylator/platinum", "broad_cytotoxic"),
    "Pemigatinib": ("FGFR", "targeted"),
    "Pralsetinib": ("RET", "targeted"),
    "Nitrogen mustard": ("Alkylator/platinum", "broad_cytotoxic"),
    "Ifosfamide": ("Alkylator/platinum", "broad_cytotoxic"),
}

BROAD_CYTOTOXIC_CLASSES = {
    "Alkylator/platinum",
    "Antimetabolite",
    "Topo/anthracycline",
    "Tubulin",
}
TARGETED_CLASSES = {
    "ALK/MET",
    "BTK",
    "CDK4/6",
    "EGFR/HER",
    "HDAC",
    "MEK",
    "Multikinase/VEGFR",
    "PARP",
    "Proteasome",
}


def classify_missing_drug(drug: str) -> tuple[str, str]:
    try:
        return MISSING_DRUG_CLASSES[drug]
    except KeyError as error:
        raise ValueError(f"Missing drug lacks an explicit class audit: {drug}") from error


def top_k_jaccard_from_ranks(
    first: Mapping[str, int],
    second: Mapping[str, int],
    k: int,
) -> float:
    if set(first) != set(second):
        raise ValueError("Rank maps must cover the same drugs")
    first_top = {drug for drug, rank in first.items() if rank <= k}
    second_top = {drug for drug, rank in second.items() if rank <= k}
    return len(first_top & second_top) / len(first_top | second_top)


def _ordinal_ranks(scores: Mapping[str, float]) -> dict[str, int]:
    ordered = sorted(scores, key=lambda drug: (-scores[drug], drug))
    return {drug: index for index, drug in enumerate(ordered, start=1)}


def _percentiles(scores: Mapping[str, float]) -> dict[str, float]:
    drugs = sorted(scores)
    values = np.asarray([scores[drug] for drug in drugs], dtype=float)
    ranked = (rankdata(values, method="average") - 1.0) / (len(values) - 1.0)
    return {
        drug: float(value)
        for drug, value in zip(drugs, ranked, strict=True)
    }


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"Input has no rows: {path}")
    return rows


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
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


def _coarse_group(mechanism_class: str) -> str:
    if mechanism_class in BROAD_CYTOTOXIC_CLASSES:
        return "broad_cytotoxic"
    if mechanism_class in TARGETED_CLASSES:
        return "targeted"
    return "other"


def run_analysis(
    project_root: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    root = project_root.resolve()
    inputs = load_inputs(root)
    all_drugs = tuple(sorted(inputs.associations))
    primary_drugs = tuple(
        drug for drug in all_drugs if inputs.mipe_mean_zauc[drug] is not None
    )
    missing_drugs = tuple(
        drug for drug in all_drugs if inputs.mipe_mean_zauc[drug] is None
    )
    if set(missing_drugs) != set(MISSING_DRUG_CLASSES):
        raise ValueError("The explicit 16-drug missingness audit is out of sync")

    background = float(
        np.mean(
            [
                inputs.r_acc[gene]
                for drug in all_drugs
                for gene in inputs.associations[drug]
            ]
        )
    )
    c_acc_by_k: dict[float, dict[str, float]] = {}
    ranks_all_by_k: dict[float, dict[str, int]] = {}
    ranks_primary_by_k: dict[float, dict[str, int]] = {}
    percentiles_primary_by_k: dict[float, dict[str, float]] = {}
    for pseudo_count in PSEUDO_COUNTS:
        c_acc: dict[str, float] = {}
        for drug in all_drugs:
            genes = inputs.associations[drug]
            target_sum = sum(inputs.r_acc[gene] for gene in genes)
            c_acc[drug] = (
                target_sum + pseudo_count * background
            ) / (len(genes) + pseudo_count)
        c_acc_by_k[pseudo_count] = c_acc
        ranks_all_by_k[pseudo_count] = _ordinal_ranks(c_acc)
        primary_scores = {drug: c_acc[drug] for drug in primary_drugs}
        ranks_primary_by_k[pseudo_count] = _ordinal_ranks(primary_scores)
        percentiles_primary_by_k[pseudo_count] = _percentiles(primary_scores)

    reference_pct = percentiles_primary_by_k[3.0]
    reference_ranks = ranks_primary_by_k[3.0]
    summary_rows: list[dict[str, Any]] = []
    for pseudo_count in PSEUDO_COUNTS:
        current_pct = percentiles_primary_by_k[pseudo_count]
        current_ranks = ranks_primary_by_k[pseudo_count]
        shifts = np.asarray(
            [
                abs(current_ranks[drug] - reference_ranks[drug])
                for drug in primary_drugs
            ],
            dtype=float,
        )
        summary_rows.append(
            {
                "pseudo_count_k": pseudo_count,
                "spearman_C_ACC_pct_vs_k3": float(
                    stats.spearmanr(
                        [reference_pct[drug] for drug in primary_drugs],
                        [current_pct[drug] for drug in primary_drugs],
                    ).statistic
                ),
                "top20_jaccard_vs_k3": top_k_jaccard_from_ranks(
                    reference_ranks,
                    current_ranks,
                    20,
                ),
                "median_absolute_rank_shift_vs_k3": float(np.median(shifts)),
                "maximum_absolute_rank_shift_vs_k3": int(shifts.max()),
            }
        )

    drug_rows: list[dict[str, Any]] = []
    for drug in all_drugs:
        row: dict[str, Any] = {
            "drug": drug,
            "MIPE_available": inputs.mipe_mean_zauc[drug] is not None,
            "n_assoc": len(inputs.associations[drug]),
        }
        for pseudo_count in PSEUDO_COUNTS:
            label = str(int(pseudo_count))
            row[f"C_ACC_k{label}"] = c_acc_by_k[pseudo_count][drug]
            row[f"rank_all124_k{label}"] = ranks_all_by_k[pseudo_count][drug]
            row[f"rank_primary108_k{label}"] = (
                ranks_primary_by_k[pseudo_count].get(drug, "")
            )
        drug_rows.append(row)

    primary_members = _read_csv(
        root / "results/mechanism_enrichment/mechanism_members_primary108.csv"
    )
    primary_class = {
        row["drug"]: row["mechanism_class"] for row in primary_members
    }
    missing_rows = []
    for drug in sorted(
        missing_drugs,
        key=lambda item: ranks_all_by_k[3.0][item],
    ):
        mechanism_class, coarse_group = classify_missing_drug(drug)
        missing_rows.append(
            {
                "drug": drug,
                "mechanism_class": mechanism_class,
                "coarse_group": coarse_group,
                "n_assoc": len(inputs.associations[drug]),
                "C_ACC_k3": c_acc_by_k[3.0][drug],
                "rank_C_ACC_all124_k3": ranks_all_by_k[3.0][drug],
                "missing_reason": "MIPE activity unavailable",
                "classification_basis": "explicit drug-identity audit",
            }
        )
    primary_groups = {
        drug: _coarse_group(primary_class[drug]) for drug in primary_drugs
    }
    missing_group_counts = Counter(row["coarse_group"] for row in missing_rows)
    primary_group_counts = Counter(primary_groups.values())
    contingency = np.asarray(
        [
            [
                missing_group_counts["targeted"],
                len(missing_drugs) - missing_group_counts["targeted"],
            ],
            [
                primary_group_counts["targeted"],
                len(primary_drugs) - primary_group_counts["targeted"],
            ],
        ],
        dtype=int,
    )
    fisher_odds, fisher_p = stats.fisher_exact(contingency)
    missing_ranks = np.asarray(
        [ranks_all_by_k[3.0][drug] for drug in missing_drugs],
        dtype=float,
    )
    primary_ranks = np.asarray(
        [ranks_all_by_k[3.0][drug] for drug in primary_drugs],
        dtype=float,
    )
    rank_test = stats.mannwhitneyu(
        missing_ranks,
        primary_ranks,
        alternative="two-sided",
    )
    missing_n_assoc = np.asarray(
        [len(inputs.associations[drug]) for drug in missing_drugs],
        dtype=float,
    )
    primary_n_assoc = np.asarray(
        [len(inputs.associations[drug]) for drug in primary_drugs],
        dtype=float,
    )
    assoc_test = stats.mannwhitneyu(
        missing_n_assoc,
        primary_n_assoc,
        alternative="two-sided",
    )
    metrics = {
        "analysis_version": ANALYSIS_VERSION,
        "pseudo_counts": list(PSEUDO_COUNTS),
        "all_drug_n": len(all_drugs),
        "primary_drug_n": len(primary_drugs),
        "MIPE_missing_n": len(missing_drugs),
        "missing_group_counts": dict(missing_group_counts),
        "primary_group_counts": dict(primary_group_counts),
        "targeted_missingness_fisher_exact": {
            "contingency_missing_vs_primary_targeted_other": contingency.tolist(),
            "odds_ratio": float(fisher_odds),
            "p_two_sided": float(fisher_p),
        },
        "C_ACC_rank_missing_vs_primary_mannwhitney": {
            "missing_median_rank": float(np.median(missing_ranks)),
            "primary_median_rank": float(np.median(primary_ranks)),
            "U": float(rank_test.statistic),
            "p_two_sided": float(rank_test.pvalue),
        },
        "association_count_missing_vs_primary_mannwhitney": {
            "missing_median_n_assoc": float(np.median(missing_n_assoc)),
            "primary_median_n_assoc": float(np.median(primary_n_assoc)),
            "U": float(assoc_test.statistic),
            "p_two_sided": float(assoc_test.pvalue),
        },
        "interpretation_boundary": (
            "Mechanism groups for the 16 missing drugs are an explicit "
            "drug-identity audit. These tests characterize missingness and "
            "cannot prove MCAR, MAR or MNAR."
        ),
    }
    return summary_rows, drug_rows, missing_rows, metrics


def write_outputs(
    project_root: Path,
    output_dir: Path,
    summary_rows: Sequence[Mapping[str, Any]],
    drug_rows: Sequence[Mapping[str, Any]],
    missing_rows: Sequence[Mapping[str, Any]],
    metrics: Mapping[str, Any],
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "shrinkage_summary": output_dir / "shrinkage_k_sensitivity_summary.csv",
        "shrinkage_drugs": output_dir / "shrinkage_k_sensitivity_drug124.csv",
        "missing_drugs": output_dir / "MIPE_missing16_audit.csv",
        "metrics": output_dir / "shrinkage_missingness_metrics.json",
        "report": output_dir / "shrinkage_missingness_report.md",
        "manifest": output_dir / "run_manifest.md",
    }
    _write_csv(paths["shrinkage_summary"], summary_rows)
    _write_csv(paths["shrinkage_drugs"], drug_rows)
    _write_csv(paths["missing_drugs"], missing_rows)
    paths["metrics"].write_text(
        json.dumps(metrics, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    lines = [
        "# Shrinkage and MIPE-missingness audit",
        "",
        "## Pseudo-count sensitivity",
        "",
        "| k | Spearman vs k=3 | Top-20 Jaccard | Median absolute rank shift | Maximum absolute rank shift |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['pseudo_count_k']:g} | "
            f"{row['spearman_C_ACC_pct_vs_k3']:.3f} | "
            f"{row['top20_jaccard_vs_k3']:.3f} | "
            f"{row['median_absolute_rank_shift_vs_k3']:.1f} | "
            f"{row['maximum_absolute_rank_shift_vs_k3']} |"
        )
    fisher = metrics["targeted_missingness_fisher_exact"]
    rank_test = metrics["C_ACC_rank_missing_vs_primary_mannwhitney"]
    lines.extend(
        [
            "",
            "## MIPE missingness",
            "",
            f"- Missing drugs: {metrics['MIPE_missing_n']}/"
            f"{metrics['all_drug_n']}.",
            f"- Targeted-class missingness Fisher exact odds ratio: "
            f"{fisher['odds_ratio']:.3f}; two-sided P={fisher['p_two_sided']:.4f}.",
            f"- Median all-124 C_ACC rank: missing "
            f"{rank_test['missing_median_rank']:.1f}, observed "
            f"{rank_test['primary_median_rank']:.1f}; Mann-Whitney "
            f"P={rank_test['p_two_sided']:.4f}.",
            "",
            "These descriptive tests do not establish a missingness mechanism. "
            "The 16-drug table is retained so readers can assess the excluded "
            "high-context candidates directly.",
            "",
        ]
    )
    paths["report"].write_text("\n".join(lines), encoding="utf-8")
    paths["manifest"].write_text(
        "\n".join(
            [
                "# Shrinkage and missingness run manifest",
                "",
                f"- Analysis version: `{ANALYSIS_VERSION}`",
                "- Command: `python -m analysis.shrinkage_missingness_audit "
                "--project-root .`",
                f"- Python: `{platform.python_version()}`",
                f"- NumPy: `{np.__version__}`",
                f"- SciPy: `{scipy.__version__}`",
                "",
                "## Outputs",
                "",
                *[
                    f"- `{path.name}`"
                    for key, path in paths.items()
                    if key != "manifest"
                ],
                "",
            ]
        ),
        encoding="utf-8",
    )
    return paths


def run(
    project_root: Path,
    output_dir: Path | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    root = project_root.resolve()
    summary, drugs, missing, metrics = run_analysis(root)
    target = (
        output_dir.resolve()
        if output_dir is not None
        else root / "results" / "reviewer_minor_audits"
    )
    write_outputs(root, target, summary, drugs, missing, metrics)
    return summary, drugs, missing, metrics


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    summary, _, _, metrics = run(
        args.project_root,
        output_dir=args.output_dir,
    )
    print(
        json.dumps(
            {
                "status": "ok",
                "analysis_version": ANALYSIS_VERSION,
                "MIPE_missing_n": metrics["MIPE_missing_n"],
                "shrinkage_summary": summary,
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
