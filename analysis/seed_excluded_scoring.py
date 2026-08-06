"""Reproduce Protocol Amendment 5 seed-excluded drug scoring.

The implementation follows the frozen protocol literally: it reads the
archived 108-drug primary table, removes the 45 disease-biology seeds from
each fixed drug association set, retains the locked pseudo-count/reference
mean, and combines the new context percentile with the archived residual
percentile.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from scipy.stats import spearmanr

from analysis.acc_primary_pipeline import load_inputs, percentile_average
from analysis.method_strengthening import load_disease_seed_weights


ANALYSIS_VERSION = "seed-excluded-scoring-v1"
PROTOCOL_PATH = Path("experiments/amendment5_seed_excluded_scoring_protocol_v1.md")
PRIMARY_PATH = Path("results/primary_analysis/ADRS_comp_primary_108.csv")
OUTPUT_DIR = Path("results/seed_excluded_scoring")
PSEUDO_COUNT = 3.0
TOP_K = 20


@dataclass(frozen=True)
class SeedExcludedResult:
    rows: tuple[dict[str, Any], ...]
    metrics: dict[str, Any]


def _read_csv(path: Path, required: Iterable[str]) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Required input is missing: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"Input has no rows: {path}")
    missing = set(required) - set(rows[0])
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    return rows


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _ordinal_ranks(scores: Mapping[str, float]) -> dict[str, int]:
    ordered = sorted(scores, key=lambda drug: (-scores[drug], drug))
    return {drug: index for index, drug in enumerate(ordered, start=1)}


def _top_k(ranks: Mapping[str, int], k: int = TOP_K) -> set[str]:
    return {drug for drug, rank in ranks.items() if rank <= k}


def _jaccard(first: set[str], second: set[str]) -> float:
    return len(first & second) / len(first | second)


def run_analysis(project_root: Path) -> SeedExcludedResult:
    root = project_root.resolve()
    inputs = load_inputs(root)
    seeds = frozenset(
        load_disease_seed_weights(
            root / "data" / "ACC_P0.5C_gene_weights_v1.csv"
        )
    )
    primary_rows = _read_csv(
        root / PRIMARY_PATH,
        (
            "rank_comp",
            "drug",
            "C_ACC",
            "C_ACC_pct",
            "residual_pct",
            "ADRS_comp",
        ),
    )
    primary_by_drug = {row["drug"].strip(): row for row in primary_rows}
    drugs = tuple(sorted(primary_by_drug))
    if len(drugs) != 108:
        raise ValueError(f"Expected 108 locked drugs, got {len(drugs)}")
    if not set(drugs).issubset(inputs.associations):
        raise ValueError("Primary drug universe is not contained in the association map")

    mu_0 = float(
        np.mean(
            [
                inputs.r_acc[gene]
                for drug in sorted(inputs.associations)
                for gene in inputs.associations[drug]
            ]
        )
    )
    if not math.isfinite(mu_0):
        raise ValueError("Association-weighted reference mean is not finite")

    nonseed_scores: dict[str, float] = {}
    seed_gene_count: dict[str, int] = {}
    nonseed_gene_count: dict[str, int] = {}
    for drug in drugs:
        genes = inputs.associations[drug]
        retained = genes - seeds
        seed_gene_count[drug] = len(genes & seeds)
        nonseed_gene_count[drug] = len(retained)
        nonseed_scores[drug] = (
            sum(inputs.r_acc[gene] for gene in retained) + PSEUDO_COUNT * mu_0
        ) / (len(retained) + PSEUDO_COUNT)

    nonseed_context_pct = percentile_average(nonseed_scores)
    nonseed_composite = {
        drug: 0.5 * nonseed_context_pct[drug]
        + 0.5 * float(primary_by_drug[drug]["residual_pct"])
        for drug in drugs
    }
    rank_nonseed_context = _ordinal_ranks(nonseed_scores)
    rank_nonseed_composite = _ordinal_ranks(nonseed_composite)
    rank_locked_context = _ordinal_ranks(
        {drug: float(primary_by_drug[drug]["C_ACC"]) for drug in drugs}
    )
    rank_locked_composite = {
        drug: int(primary_by_drug[drug]["rank_comp"]) for drug in drugs
    }

    rows: list[dict[str, Any]] = []
    for drug in sorted(drugs, key=lambda item: rank_nonseed_composite[item]):
        rows.append(
            {
                "drug": drug,
                "seed_gene_count": seed_gene_count[drug],
                "nonseed_gene_count": nonseed_gene_count[drug],
                "C_ACC_nonseed": nonseed_scores[drug],
                "C_ACC_nonseed_pct": nonseed_context_pct[drug],
                "ADRS_nonseed": nonseed_composite[drug],
                "rank_locked_context": rank_locked_context[drug],
                "rank_nonseed_context": rank_nonseed_context[drug],
                "rank_locked_composite": rank_locked_composite[drug],
                "rank_nonseed_composite": rank_nonseed_composite[drug],
                "composite_rank_change": (
                    rank_nonseed_composite[drug] - rank_locked_composite[drug]
                ),
            }
        )

    locked_context = np.asarray(
        [float(primary_by_drug[drug]["C_ACC"]) for drug in drugs], dtype=float
    )
    remedial_context = np.asarray([nonseed_scores[drug] for drug in drugs])
    locked_composite = np.asarray(
        [float(primary_by_drug[drug]["ADRS_comp"]) for drug in drugs], dtype=float
    )
    remedial_composite = np.asarray([nonseed_composite[drug] for drug in drugs])
    context_rho = float(spearmanr(locked_context, remedial_context).statistic)
    composite_rho = float(
        spearmanr(locked_composite, remedial_composite).statistic
    )

    locked_context_top = _top_k(rank_locked_context)
    nonseed_context_top = _top_k(rank_nonseed_context)
    locked_composite_top = _top_k(rank_locked_composite)
    nonseed_composite_top = _top_k(rank_nonseed_composite)
    context_intersection = len(locked_context_top & nonseed_context_top)
    composite_intersection = len(locked_composite_top & nonseed_composite_top)
    shifts = {
        drug: abs(rank_nonseed_composite[drug] - rank_locked_composite[drug])
        for drug in drugs
    }
    exposed = [drug for drug in drugs if seed_gene_count[drug] > 0]
    unexposed = [drug for drug in drugs if seed_gene_count[drug] == 0]
    metrics: dict[str, Any] = {
        "analysis_version": ANALYSIS_VERSION,
        "drug_count": len(drugs),
        "seed_count": len(seeds),
        "pseudo_count": PSEUDO_COUNT,
        "mu_0": mu_0,
        "context_spearman": context_rho,
        "context_top20_intersection": context_intersection,
        "context_top20_jaccard": _jaccard(
            locked_context_top, nonseed_context_top
        ),
        "composite_spearman": composite_rho,
        "composite_top20_intersection": composite_intersection,
        "composite_top20_jaccard": _jaccard(
            locked_composite_top, nonseed_composite_top
        ),
        "entered_top20": sorted(nonseed_composite_top - locked_composite_top),
        "left_top20": sorted(locked_composite_top - nonseed_composite_top),
        "directly_exposed_drug_count": len(exposed),
        "unexposed_drug_count": len(unexposed),
        "zero_nonseed_gene_count": sum(
            count == 0 for count in nonseed_gene_count.values()
        ),
        "one_nonseed_gene_count": sum(
            count == 1 for count in nonseed_gene_count.values()
        ),
        "at_least_two_nonseed_gene_count": sum(
            count >= 2 for count in nonseed_gene_count.values()
        ),
        "maximum_composite_rank_shift": max(shifts.values()),
        "exposed_median_abs_shift": float(
            np.median([shifts[drug] for drug in exposed])
        ),
        "exposed_max_abs_shift": max(shifts[drug] for drug in exposed),
        "unexposed_median_abs_shift": float(
            np.median([shifts[drug] for drug in unexposed])
        ),
        "unexposed_max_abs_shift": max(shifts[drug] for drug in unexposed),
        "protocol_sha256": _sha256(root / PROTOCOL_PATH),
        "input_sha256": {
            str(path.as_posix()): _sha256(root / path)
            for path in (
                Path("data/bindex_network/bindex_edges_1304.csv"),
                Path("data/bindex_network/rACC_399_fullSTRING.csv"),
                Path("data/ACC_P0.5C_gene_weights_v1.csv"),
                PRIMARY_PATH,
            )
        },
    }
    return SeedExcludedResult(rows=tuple(rows), metrics=metrics)


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("Cannot write an empty CSV")
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


def _format_drug_list(drugs: Sequence[str]) -> str:
    return ", ".join(drugs) if drugs else "None"


def _render_report(result: SeedExcludedResult) -> str:
    m = result.metrics
    rows = {row["drug"]: row for row in result.rows}
    top20 = sorted(result.rows, key=lambda row: row["rank_nonseed_composite"])[:20]
    largest = sorted(
        result.rows,
        key=lambda row: (-abs(row["composite_rank_change"]), row["drug"]),
    )[:10]
    focal = [rows[name] for name in ("Abemaciclib", "Palbociclib", "Ribociclib")]
    lines = [
        "# Protocol Amendment 5: seed-excluded scoring audit",
        "",
        f"- Analysis implementation: `{m['analysis_version']}`.",
        f"- Protocol SHA-256: `{m['protocol_sha256']}`.",
        "- Status: result-known post-hoc; protocol frozen before outputs.",
        "",
        "## Primary results",
        "",
        "| Comparison with locked ranking | Spearman rho | Top-20 intersection | Top-20 Jaccard |",
        "|---|---:|---:|---:|",
        f"| Seed-excluded C_ACC versus locked C_ACC | {m['context_spearman']:.4f} | {m['context_top20_intersection']}/20 | {m['context_top20_jaccard']:.4f} |",
        f"| ADRS_nonseed versus locked ADRS_comp | {m['composite_spearman']:.4f} | {m['composite_top20_intersection']}/20 | {m['composite_top20_jaccard']:.4f} |",
        "",
        f"Entered: {_format_drug_list(m['entered_top20'])}.",
        "",
        f"Left: {_format_drug_list(m['left_top20'])}.",
        "",
        "## Coverage and focal drugs",
        "",
        f"Among {m['drug_count']} drugs, {m['directly_exposed_drug_count']} had at least one direct seed association and {m['unexposed_drug_count']} had none. After exclusion, {m['zero_nonseed_gene_count']} drugs had no remaining gene, {m['one_nonseed_gene_count']} had one and {m['at_least_two_nonseed_gene_count']} had at least two. The locked reference mean was mu_0 = {m['mu_0']:.10f}.",
        "",
        "| Drug | Seed genes removed | Non-seed genes retained | Locked ADRS rank | ADRS_nonseed rank | Locked C_ACC rank | Seed-excluded C_ACC rank |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in focal:
        lines.append(
            f"| {row['drug']} | {row['seed_gene_count']} | {row['nonseed_gene_count']} | {row['rank_locked_composite']} | {row['rank_nonseed_composite']} | {row['rank_locked_context']} | {row['rank_nonseed_context']} |"
        )
    lines.extend(
        [
            "",
            "## Seed-excluded Top 20",
            "",
            "| New rank | Drug | Locked rank | Seed genes removed | Non-seed genes retained | C_ACC,nonseed | ADRS_nonseed | Rank change |",
            "|---:|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in top20:
        lines.append(
            f"| {row['rank_nonseed_composite']} | {row['drug']} | {row['rank_locked_composite']} | {row['seed_gene_count']} | {row['nonseed_gene_count']} | {row['C_ACC_nonseed']:.6f} | {row['ADRS_nonseed']:.6f} | {row['composite_rank_change']:+d} |"
        )
    lines.extend(
        [
            "",
            "## Largest composite-rank movements",
            "",
            "| Drug | Locked rank | ADRS_nonseed rank | Absolute shift | Seed genes removed | Non-seed genes retained |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in largest:
        lines.append(
            f"| {row['drug']} | {row['rank_locked_composite']} | {row['rank_nonseed_composite']} | {abs(row['composite_rank_change'])} | {row['seed_gene_count']} | {row['nonseed_gene_count']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "This deterministic sensitivity arm shows that direct seed contributions can be removed computationally. It does not establish that the remaining score is more biologically valid or clinically predictive, and it must be reported alongside rather than silently substituted for the locked analysis.",
            "",
            "## Input SHA-256",
            "",
        ]
    )
    for path, digest in m["input_sha256"].items():
        lines.append(f"- `{path}`: `{digest}`")
    return "\n".join(lines) + "\n"


def write_outputs(
    project_root: Path,
    result: SeedExcludedResult,
    output_dir: Path | None = None,
) -> tuple[Path, Path, Path]:
    root = project_root.resolve()
    target = (root / OUTPUT_DIR) if output_dir is None else output_dir.resolve()
    target.mkdir(parents=True, exist_ok=True)
    csv_path = target / "seed_excluded_scores.csv"
    metrics_path = target / "seed_excluded_metrics.json"
    report_path = target / "seed_excluded_scoring_audit.md"
    _write_csv(csv_path, result.rows)
    metrics_path.write_text(
        json.dumps(result.metrics, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    report_path.write_text(_render_report(result), encoding="utf-8")
    return csv_path, metrics_path, report_path


def run(project_root: Path, output_dir: Path | None = None) -> SeedExcludedResult:
    result = run_analysis(project_root)
    write_outputs(project_root, result, output_dir)
    return result


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
    result = run(args.project_root, args.output_dir)
    print(json.dumps(result.metrics, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
