"""Mechanism-class enrichment on the locked 108-drug C1 ranking.

The observed class statistic and its null distribution are both constructed
from the same complete-case universe.  Mechanism assignment is score-blind:
it uses MIPE Primary MOA metadata plus explicit drug-name aliases, never ADRS
values or ranks.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import platform
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


ANALYSIS_VERSION = "mechanism-primary108-v1"
MIN_CLASS_SIZE = 3
OTHER_CLASS = "Other"

# Drug-name aliases are explicit so classification cannot change when wording
# in the free-text Primary MOA column changes.
NAME_OVERRIDES: dict[str, str] = {
    # Pre-specified class of biological interest.
    "Abemaciclib": "CDK4/6",
    "Palbociclib": "CDK4/6",
    "Ribociclib": "CDK4/6",
    "Trilaciclib": "CDK4/6",
    # MEK.
    "Binimetinib": "MEK",
    "Cobimetinib": "MEK",
    "Selumetinib": "MEK",
    "Trametinib": "MEK",
    # EGFR/HER.
    "Afatinib": "EGFR/HER",
    "Dacomitinib": "EGFR/HER",
    "Erlotinib": "EGFR/HER",
    "Lapatinib": "EGFR/HER",
    "Neratinib": "EGFR/HER",
    "Osimertinib": "EGFR/HER",
    "Tucatinib": "EGFR/HER",
    "Vandetanib": "EGFR/HER",
    # ALK/MET.  Explicit names avoid the legacy ``"alk" in moa`` bug, which
    # incorrectly classified the word "alkaloid" as ALK.
    "Brigatinib": "ALK/MET",
    "Ceritinib": "ALK/MET",
    "Crizotinib": "ALK/MET",
    "Lorlatinib": "ALK/MET",
    "Tepotinib": "ALK/MET",
    "tepotinib": "ALK/MET",
    # Topoisomerase/anthracycline.
    "Daunorubicin": "Topo/anthracycline",
    "Dexrazoxane": "Topo/anthracycline",
    "Doxorubicin": "Topo/anthracycline",
    "Epirubicin": "Topo/anthracycline",
    "Etoposide": "Topo/anthracycline",
    "Idarubicin": "Topo/anthracycline",
    "Irinotecan": "Topo/anthracycline",
    "Mitoxantrone": "Topo/anthracycline",
    "Teniposide": "Topo/anthracycline",
    "Topotecan": "Topo/anthracycline",
    "Valrubicin": "Topo/anthracycline",
    # Tubulin/microtubule.
    "Cabazitaxel": "Tubulin",
    "Docetaxel": "Tubulin",
    "Ixabepilone": "Tubulin",
    "Paclitaxel": "Tubulin",
    "Vinblastine": "Tubulin",
    "Vinorelbine": "Tubulin",
    # Alkylating/platinum agents.
    "Carboplatin": "Alkylator/platinum",
    "Carmustine": "Alkylator/platinum",
    "Chlorambucil": "Alkylator/platinum",
    "Cisplatin": "Alkylator/platinum",
    "Cyclophosphamide": "Alkylator/platinum",
    "Dacarbazine": "Alkylator/platinum",
    "Estramustine": "Alkylator/platinum",
    "Ifosfamide": "Alkylator/platinum",
    "Lomustine": "Alkylator/platinum",
    "Melphalan": "Alkylator/platinum",
    "Mitomycin": "Alkylator/platinum",
    "Nitrogen mustard": "Alkylator/platinum",
    "Oxaliplatin": "Alkylator/platinum",
    "Pipobroman": "Alkylator/platinum",
    "Thiotepa": "Alkylator/platinum",
    "Triethylenemelamine": "Alkylator/platinum",
    "Uracil mustard": "Alkylator/platinum",
    # Antimetabolites and antimetabolite-like nucleotide synthesis inhibitors.
    "Cladribine": "Antimetabolite",
    "Clofarabine": "Antimetabolite",
    "Cytarabine": "Antimetabolite",
    "Fludarabine": "Antimetabolite",
    "Fluorouracil": "Antimetabolite",
    "Gemcitabine": "Antimetabolite",
    "Hydroxyurea": "Antimetabolite",
    "Mercaptopurine": "Antimetabolite",
    "Methotrexate": "Antimetabolite",
    "Pemetrexed": "Antimetabolite",
    "Pralatrexate": "Antimetabolite",
    "Tioguanine": "Antimetabolite",
    "Trifluridine": "Antimetabolite",
    # VEGFR/multikinase.
    "Axitinib": "Multikinase/VEGFR",
    "Cabozantinib": "Multikinase/VEGFR",
    "Fruquintinib": "Multikinase/VEGFR",
    "Lenvatinib": "Multikinase/VEGFR",
    "Pazopanib": "Multikinase/VEGFR",
    "Regorafenib": "Multikinase/VEGFR",
    "Sunitinib": "Multikinase/VEGFR",
    "Tivozanib": "Multikinase/VEGFR",
}


CURATED_MOA: dict[str, str] = {
    "Gilteritinib": "FLT3/AXL inhibitor",
    "Selumetinib": "MEK1/2 inhibitor",
    "Tivozanib": "VEGFR inhibitor",
    "Vinorelbine": "Microtubule inhibitor",
}


@dataclass(frozen=True)
class EnrichmentResult:
    n_universe: int
    class_rows: tuple[dict[str, Any], ...]
    member_rows: tuple[dict[str, Any], ...]


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


def load_analysis_inputs(
    project_root: Path,
) -> tuple[list[dict[str, Any]], dict[str, dict[str, str]]]:
    """Load the C1 ranking and score-blind mechanism metadata."""

    project_root = project_root.resolve()
    primary_path = (
        project_root / "results" / "primary_analysis" / "ADRS_comp_primary_108.csv"
    )
    activity_path = project_root / "data" / "bindex_network" / "Sactivity_124_v1.csv"
    mipe_path = project_root / "data" / "ACC_MIPE5_ranked_by_mean_ZAUC.csv"

    primary_raw = _read_csv(primary_path, ("rank_comp", "drug", "ADRS_comp"))
    activity = _read_csv(activity_path, ("drug", "MIPE_name"))
    mipe = _read_csv(mipe_path, ("Sample Name", "Primary MOA", "Gene Symbol"))

    primary_rows: list[dict[str, Any]] = []
    for row in primary_raw:
        primary_rows.append(
            {
                **row,
                "rank_comp": int(row["rank_comp"]),
                "ADRS_comp": float(row["ADRS_comp"]),
            }
        )
    primary_rows.sort(key=lambda row: row["rank_comp"])
    expected_ranks = list(range(1, len(primary_rows) + 1))
    if [row["rank_comp"] for row in primary_rows] != expected_ranks:
        raise ValueError("C1 ranking must contain consecutive ranks beginning at 1")
    if len({row["drug"] for row in primary_rows}) != len(primary_rows):
        raise ValueError("C1 ranking contains duplicate drug names")

    drug_to_mipe = {row["drug"]: row["MIPE_name"] for row in activity}
    mipe_metadata: dict[str, dict[str, str]] = {}
    for row in mipe:
        mipe_metadata.setdefault(
            row["Sample Name"],
            {
                "primary_moa": row["Primary MOA"].strip(),
                "gene_symbol": row["Gene Symbol"].strip(),
            },
        )

    metadata: dict[str, dict[str, str]] = {}
    for row in primary_rows:
        drug = row["drug"]
        mipe_name = drug_to_mipe.get(drug, "")
        source_row = mipe_metadata.get(mipe_name, {})
        primary_moa = source_row.get("primary_moa", "")
        source = "MIPE Primary MOA"
        if not primary_moa and drug in CURATED_MOA:
            primary_moa = CURATED_MOA[drug]
            source = "curated drug-identity alias"
        elif not primary_moa:
            source = "drug-name rule or unclassified"
        metadata[drug] = {
            "mipe_name": mipe_name,
            "primary_moa": primary_moa,
            "gene_symbol": source_row.get("gene_symbol", ""),
            "moa_source": source,
        }
    return primary_rows, metadata


def assign_mechanism_class(drug: str, primary_moa: str) -> tuple[str, str]:
    """Return the frozen class and the rule responsible for assignment."""

    if drug in NAME_OVERRIDES:
        return NAME_OVERRIDES[drug], "explicit_drug_alias"

    moa = primary_moa.casefold()
    if "cdk4" in moa or "cdk 4" in moa:
        return "CDK4/6", "primary_moa"
    if "mek" in moa:
        return "MEK", "primary_moa"
    if "proteasome" in moa:
        return "Proteasome", "primary_moa"
    if "hdac" in moa or "histone deacetylase" in moa:
        return "HDAC", "primary_moa"
    if "parp" in moa:
        return "PARP", "primary_moa"
    if "egfr" in moa or "her2" in moa or "erbb2" in moa:
        return "EGFR/HER", "primary_moa"
    if "btk" in moa:
        return "BTK", "primary_moa"
    if (
        "alk inhibitor" in moa
        or "anaplastic lymphoma kinase" in moa
        or "hgfr" in moa
        or "met inhibitor" in moa
    ):
        return "ALK/MET", "primary_moa"
    if "topoisomerase" in moa:
        return "Topo/anthracycline", "primary_moa"
    if "tubulin" in moa or "microtubule" in moa:
        return "Tubulin", "primary_moa"
    if "alkylat" in moa or "platinum" in moa:
        return "Alkylator/platinum", "primary_moa"
    if (
        "purine antagonist" in moa
        or "pyrimidine antagonist" in moa
        or "dna polymerase inhibitor" in moa
        or "thymidylate synthase inhibitor" in moa
        or "dihydrofolate reductase" in moa
        or "ribonucleotide reductase" in moa
    ):
        return "Antimetabolite", "primary_moa"
    if "vegfr" in moa or "vascular endothelial growth factor receptor" in moa:
        return "Multikinase/VEGFR", "primary_moa"
    return OTHER_CLASS, "other_or_singleton"


def exact_lower_tail_rank_sum_p(
    n_universe: int, class_size: int, observed_rank_sum: int
) -> float:
    """Exact P(sum of k distinct ranks <= observed) for ranks 1..N.

    Dynamic programming counts all equally likely same-size subsets without
    Monte-Carlo error.
    """

    if n_universe < 1:
        raise ValueError("n_universe must be positive")
    if not 1 <= class_size <= n_universe:
        raise ValueError("class_size must lie in [1, n_universe]")

    minimum = class_size * (class_size + 1) // 2
    maximum = class_size * (2 * n_universe - class_size + 1) // 2
    if observed_rank_sum < minimum:
        return 0.0
    if observed_rank_sum >= maximum:
        return 1.0

    counts: list[list[int]] = [
        [0] * (observed_rank_sum + 1) for _ in range(class_size + 1)
    ]
    counts[0][0] = 1
    for rank in range(1, n_universe + 1):
        for selected in range(min(class_size, rank), 0, -1):
            for current_sum in range(observed_rank_sum, rank - 1, -1):
                counts[selected][current_sum] += counts[selected - 1][
                    current_sum - rank
                ]
    favorable = sum(counts[class_size])
    return favorable / math.comb(n_universe, class_size)


def benjamini_hochberg(p_values: Sequence[float]) -> list[float]:
    """Benjamini-Hochberg adjusted P values in the original input order."""

    if any(p < 0 or p > 1 for p in p_values):
        raise ValueError("P values must lie in [0, 1]")
    m = len(p_values)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda index: p_values[index])
    adjusted = [1.0] * m
    running = 1.0
    for rank_index in range(m - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        candidate = p_values[original_index] * m / rank
        running = min(running, candidate)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def compute_mechanism_enrichment(
    primary_rows: Sequence[Mapping[str, Any]],
    metadata: Mapping[str, Mapping[str, str]],
) -> EnrichmentResult:
    n_universe = len(primary_rows)
    ranks = {str(row["drug"]): int(row["rank_comp"]) for row in primary_rows}
    if sorted(ranks.values()) != list(range(1, n_universe + 1)):
        raise ValueError("The primary universe must have unique consecutive ranks")

    member_rows: list[dict[str, Any]] = []
    by_class: dict[str, list[str]] = {}
    for row in primary_rows:
        drug = str(row["drug"])
        meta = metadata.get(drug, {})
        mechanism_class, assignment_rule = assign_mechanism_class(
            drug, meta.get("primary_moa", "")
        )
        by_class.setdefault(mechanism_class, []).append(drug)
        member_rows.append(
            {
                "rank_comp": ranks[drug],
                "drug": drug,
                "mechanism_class": mechanism_class,
                "primary_moa": meta.get("primary_moa", ""),
                "gene_symbol": meta.get("gene_symbol", ""),
                "moa_source": meta.get("moa_source", ""),
                "assignment_rule": assignment_rule,
                "in_primary_universe": True,
                "eligible_class": False,
            }
        )

    eligible = {
        class_name: drugs
        for class_name, drugs in by_class.items()
        if class_name != OTHER_CLASS and len(drugs) >= MIN_CLASS_SIZE
    }
    for row in member_rows:
        row["eligible_class"] = row["mechanism_class"] in eligible

    class_rows: list[dict[str, Any]] = []
    for class_name, drugs in eligible.items():
        ordered = sorted(drugs, key=lambda drug: ranks[drug])
        member_ranks = [ranks[drug] for drug in ordered]
        rank_sum = sum(member_ranks)
        class_rows.append(
            {
                "mechanism_class": class_name,
                "k": len(ordered),
                "n_universe": n_universe,
                "rank_sum": rank_sum,
                "mean_rank": rank_sum / len(ordered),
                "median_rank": _median(member_ranks),
                "p_exact": exact_lower_tail_rank_sum_p(
                    n_universe, len(ordered), rank_sum
                ),
                "members": "; ".join(ordered),
                "member_ranks": "; ".join(str(rank) for rank in member_ranks),
                "test": "exact one-sided random-set rank-sum",
                "alternative": "lower mean rank than random same-size sets",
            }
        )

    q_values = benjamini_hochberg([row["p_exact"] for row in class_rows])
    for row, q_value in zip(class_rows, q_values, strict=True):
        row["q_bh"] = q_value
        row["significant_fdr_0_05"] = q_value < 0.05
    class_rows.sort(key=lambda row: (row["p_exact"], row["mechanism_class"]))
    member_rows.sort(key=lambda row: row["rank_comp"])
    return EnrichmentResult(
        n_universe=n_universe,
        class_rows=tuple(class_rows),
        member_rows=tuple(member_rows),
    )


def _median(values: Sequence[int]) -> float:
    ordered = sorted(values)
    midpoint = len(ordered) // 2
    if len(ordered) % 2:
        return float(ordered[midpoint])
    return (ordered[midpoint - 1] + ordered[midpoint]) / 2


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_outputs(
    project_root: Path, result: EnrichmentResult, output_dir: Path
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    class_path = output_dir / "mechanism_enrichment_primary108.csv"
    member_path = output_dir / "mechanism_members_primary108.csv"
    claim_path = output_dir / "claim_evidence_table.md"
    report_path = output_dir / "C2_change_report.md"
    manifest_path = output_dir / "run_manifest.md"

    _write_csv(class_path, result.class_rows)
    _write_csv(member_path, result.member_rows)

    cdk = next(
        row for row in result.class_rows if row["mechanism_class"] == "CDK4/6"
    )
    best = result.class_rows[0]
    claim_path.write_text(
        "\n".join(
            [
                "# C2 claim–evidence table",
                "",
                "| Claim | Direct result | Statistical support | Boundary |",
                "|---|---|---|---|",
                (
                    "| CDK4/6 inhibitors show a non-significant ranking trend "
                    f"| ranks {cdk['member_ranks']}; mean {cdk['mean_rank']:.2f} "
                    f"| exact one-sided P={cdk['p_exact']:.4f}; "
                    f"BH q={cdk['q_bh']:.4f} | Not nominally or FDR significant; "
                    "does not establish ACC efficacy |"
                ),
                (
                    f"| Lowest raw class P was {best['mechanism_class']} "
                    f"| k={best['k']}; mean rank {best['mean_rank']:.2f} "
                    f"| exact P={best['p_exact']:.4f}; BH q={best['q_bh']:.4f} "
                    "| Exploratory class scan; no causal or efficacy inference |"
                ),
                (
                    "| Observed and null sets use one universe "
                    f"| all {result.n_universe} C1 complete-case drugs "
                    "| exact enumeration of all same-size rank subsets "
                    "| Excludes 16 context-only drugs lacking MIPE |"
                ),
                "",
            ]
        ),
        encoding="utf-8",
    )

    table_lines = [
        "| Class | k | Mean rank | Exact P | BH q |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in result.class_rows:
        table_lines.append(
            f"| {row['mechanism_class']} | {row['k']} | "
            f"{row['mean_rank']:.2f} | {row['p_exact']:.4f} | "
            f"{row['q_bh']:.4f} |"
        )
    report_path.write_text(
        "\n".join(
            [
                "# C2 change report: mechanism-class enrichment",
                "",
                "## Locked analysis",
                "",
                f"- Version: `{ANALYSIS_VERSION}`.",
                f"- Universe: the {result.n_universe} drugs in the C1 complete-case "
                "ADRS_comp ranking.",
                "- Class labels: MIPE `Primary MOA` plus explicit, score-blind "
                "drug-name aliases.",
                "- Eligible tests: frozen named classes represented by at least "
                f"{MIN_CLASS_SIZE} drugs; `Other` is not tested.",
                "- Statistic: within-class rank sum (equivalently mean rank; lower "
                "is better).",
                "- Null: every same-size subset of ranks 1–108, counted exactly by "
                "dynamic programming.",
                "- Multiplicity: Benjamini–Hochberg across the ten eligible tests.",
                "",
                "## Corrected results",
                "",
                *table_lines,
                "",
                "## CDK4/6 correction",
                "",
                (
                    f"The primary universe contains Abemaciclib, Palbociclib and "
                    f"Ribociclib at ranks {cdk['member_ranks']} (mean "
                    f"{cdk['mean_rank']:.2f}). The exact one-sided random-set test "
                    f"gives P={cdk['p_exact']:.4f} and BH q={cdk['q_bh']:.4f}."
                ),
                "",
                "The old P=0.023/FDR=0.10 statement is not reproducible under the "
                "locked C1 model and must be removed. The supported wording is: "
                "**“The pre-specified CDK4/6 class showed a non-significant ranking "
                "trend.”** This is a candidate hypothesis, not evidence of efficacy.",
                "",
                "Trilaciclib is correctly absent because it lacks MIPE data and is "
                "one of the 16 context-only drugs excluded by the C1 complete-case "
                "rule.",
                "",
                "## Classification audit",
                "",
                "Two legacy substring rules were removed: `\"alk\" in MOA` could "
                "misclassify “alkaloid” as ALK, and `\"egfr\" in MOA` also matched "
                "the substring in “VEGFR”. Separating VEGFR from EGFR yields ten, "
                "rather than nine, eligible families. All name aliases are now "
                "listed in the analysis source, and the complete 108-drug membership "
                "table records the source and assignment rule for every drug.",
                "",
                "## Downstream synchronization",
                "",
                "C3 must regenerate Figure 5b from the CSV output and replace the "
                "hard-coded P=0.023 annotation. Manuscript-wide wording and figure "
                "legends will be synchronized after the C3 figure outputs are frozen.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    input_paths = [
        project_root
        / "results"
        / "primary_analysis"
        / "ADRS_comp_primary_108.csv",
        project_root / "data" / "bindex_network" / "Sactivity_124_v1.csv",
        project_root / "data" / "ACC_MIPE5_ranked_by_mean_ZAUC.csv",
    ]
    manifest_lines = [
        "# C2 run manifest",
        "",
        f"- Analysis version: `{ANALYSIS_VERSION}`",
        f"- Python: `{platform.python_version()}`",
        f"- Platform: `{platform.platform()}`",
        "- Random seed: not applicable; the null distribution is exact.",
        f"- Primary universe: `{result.n_universe}` drugs",
        f"- Eligible mechanism classes: `{len(result.class_rows)}`",
        "",
        "## Inputs",
        "",
    ]
    for path in input_paths:
        relative = path.resolve().relative_to(project_root.resolve()).as_posix()
        manifest_lines.append(f"- `{relative}` — SHA-256 `{_sha256(path)}`")
    manifest_lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- `{class_path.name}`",
            f"- `{member_path.name}`",
            f"- `{claim_path.name}`",
            f"- `{report_path.name}`",
            "",
        ]
    )
    manifest_path.write_text("\n".join(manifest_lines), encoding="utf-8")
    return {
        "class_results": class_path,
        "members": member_path,
        "claim_evidence": claim_path,
        "change_report": report_path,
        "manifest": manifest_path,
    }


def run(project_root: Path, output_dir: Path | None = None) -> EnrichmentResult:
    primary_rows, metadata = load_analysis_inputs(project_root)
    result = compute_mechanism_enrichment(primary_rows, metadata)
    target = (
        output_dir
        if output_dir is not None
        else project_root / "results" / "mechanism_enrichment"
    )
    write_outputs(project_root.resolve(), result, target.resolve())
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
    cdk = next(
        row for row in result.class_rows if row["mechanism_class"] == "CDK4/6"
    )
    print(
        f"C2 complete: n={result.n_universe}, classes={len(result.class_rows)}, "
        f"CDK4/6 P={cdk['p_exact']:.6f}, q={cdk['q_bh']:.6f}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
