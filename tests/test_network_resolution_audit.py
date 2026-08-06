from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.network_resolution_audit import compute_resolution


def test_compute_resolution_counts_identical_and_near_duplicate_pairs(
    tmp_path: Path,
) -> None:
    edges = pd.DataFrame(
        {
            "drug": ["A", "A", "B", "B", "C", "D"],
            "gene": ["G1", "G2", "G1", "G2", "G1", "G3"],
        }
    )
    universe = pd.DataFrame({"drug": ["A", "B", "C", "D"]})
    edges_path = tmp_path / "edges.csv"
    universe_path = tmp_path / "universe.csv"
    edges.to_csv(edges_path, index=False)
    universe.to_csv(universe_path, index=False)

    pairs, summary = compute_resolution(edges_path, universe_path)

    assert len(pairs) == 6
    assert summary["n_locked_drugs"] == 4
    assert summary["n_singleton_drugs"] == 2
    assert summary["n_identical_pairs"] == 1
    assert summary["identical_pairs"][0]["drug_a"] == "a"
    assert summary["identical_pairs"][0]["drug_b"] == "b"
    assert summary["n_pairs_jaccard_ge_0_8"] == 1


def test_compute_resolution_rejects_missing_universe_drug(tmp_path: Path) -> None:
    edges = pd.DataFrame({"drug": ["A"], "gene": ["G1"]})
    universe = pd.DataFrame({"drug": ["A", "B"]})
    edges_path = tmp_path / "edges.csv"
    universe_path = tmp_path / "universe.csv"
    edges.to_csv(edges_path, index=False)
    universe.to_csv(universe_path, index=False)

    try:
        compute_resolution(edges_path, universe_path)
    except ValueError as error:
        assert "locked-universe drugs lack association sets" in str(error)
    else:
        raise AssertionError("expected missing universe drug to be rejected")
