"""Quantify drug-resolution limits imposed by the locked association network."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
from typing import Sequence

import pandas as pd


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def compute_resolution(
    edges_path: Path,
    universe_path: Path,
) -> tuple[pd.DataFrame, dict[str, object]]:
    edges = pd.read_csv(edges_path)
    required_edges = {"drug", "gene"}
    if not required_edges.issubset(edges.columns):
        raise ValueError(f"edge table must contain {sorted(required_edges)}")

    universe_table = pd.read_csv(universe_path)
    if "drug" not in universe_table.columns:
        raise ValueError("universe table must contain a drug column")

    normalized = pd.DataFrame(
        {
            "drug": edges["drug"].astype(str).str.strip().str.lower(),
            "gene": edges["gene"].astype(str).str.strip().str.upper(),
        }
    ).drop_duplicates()
    universe = set(
        universe_table["drug"].astype(str).str.strip().str.lower()
    )
    gene_sets = {
        drug: set(group["gene"])
        for drug, group in normalized.groupby("drug")
        if drug in universe
    }
    missing = sorted(universe - set(gene_sets))
    if missing:
        raise ValueError(f"{len(missing)} locked-universe drugs lack association sets")

    rows: list[dict[str, object]] = []
    for drug_a, drug_b in itertools.combinations(sorted(gene_sets), 2):
        genes_a = gene_sets[drug_a]
        genes_b = gene_sets[drug_b]
        intersection = genes_a & genes_b
        union = genes_a | genes_b
        rows.append(
            {
                "drug_a": drug_a,
                "drug_b": drug_b,
                "intersection_n": len(intersection),
                "union_n": len(union),
                "jaccard": len(intersection) / len(union),
                "n_a": len(genes_a),
                "n_b": len(genes_b),
                "shared_genes": ";".join(sorted(intersection)),
            }
        )

    pairs = pd.DataFrame(rows).sort_values(
        ["jaccard", "drug_a", "drug_b"],
        ascending=[False, True, True],
        kind="mergesort",
    )
    identical = pairs.loc[pairs["jaccard"] == 1.0]
    near_duplicate = pairs.loc[pairs["jaccard"] >= 0.8]
    singleton_drugs = sorted(
        drug for drug, genes in gene_sets.items() if len(genes) == 1
    )

    summary: dict[str, object] = {
        "analysis": "locked_drug_association_set_resolution",
        "n_locked_drugs": len(gene_sets),
        "n_pairwise_comparisons": len(pairs),
        "n_singleton_drugs": len(singleton_drugs),
        "singleton_drugs": singleton_drugs,
        "n_identical_pairs": int(len(identical)),
        "identical_pairs": identical[
            ["drug_a", "drug_b", "n_a", "shared_genes"]
        ].to_dict(orient="records"),
        "n_pairs_jaccard_ge_0_8": int(len(near_duplicate)),
        "n_pairs_jaccard_ge_0_5": int((pairs["jaccard"] >= 0.5).sum()),
        "edges_sha256": _sha256(edges_path),
        "universe_sha256": _sha256(universe_path),
    }
    return pairs, summary


def build(project_root: Path) -> dict[str, Path]:
    root = project_root.resolve()
    edges_path = root / "data" / "bindex_network" / "bindex_edges_1304.csv"
    universe_path = root / "leakage_audit" / "observational" / "Ld_108.csv"
    output_dir = root / "leakage_audit" / "observational"
    output_dir.mkdir(parents=True, exist_ok=True)

    pairs, summary = compute_resolution(edges_path, universe_path)
    pairs_path = output_dir / "pairwise_association_jaccard.csv"
    summary_path = output_dir / "network_resolution_summary.json"
    pairs.to_csv(pairs_path, index=False, encoding="utf-8-sig")
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return {"pairs": pairs_path, "summary": summary_path}


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    outputs = build(args.project_root)
    for label, path in outputs.items():
        print(f"{label}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
