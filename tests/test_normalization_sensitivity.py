import numpy as np
import pytest
from scipy import sparse
import csv
import json
from pathlib import Path

from analysis.normalization_sensitivity import (
    network_smooth_with_restart,
    top_k_jaccard,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_network_smooth_with_restart_converges_on_symmetric_operator() -> None:
    operator = sparse.csr_matrix(
        np.asarray(
            [
                [0.0, 1.0],
                [1.0, 0.0],
            ]
        )
    )
    restart = np.asarray([1.0, 0.0])

    propagated, iterations, delta = network_smooth_with_restart(
        operator,
        restart,
        alpha=0.4,
    )

    assert propagated.shape == (2, 1)
    assert iterations < 500
    assert delta < 1e-10
    assert propagated[:, 0] == pytest.approx([0.625, 0.375], abs=1e-9)


def test_top_k_jaccard_is_deterministic_under_ties() -> None:
    names = ("b", "a", "c", "d")
    first = np.asarray([1.0, 1.0, 0.5, 0.0])
    second = np.asarray([0.5, 1.0, 1.0, 0.0])

    assert top_k_jaccard(names, first, second, k=2) == pytest.approx(1 / 3)


def test_formal_normalization_sensitivity_outputs() -> None:
    output_dir = PROJECT_ROOT / "results/normalization_sensitivity"
    metrics = json.loads(
        (output_dir / "normalization_sensitivity_metrics.json").read_text(
            encoding="utf-8"
        )
    )
    with (output_dir / "normalization_sensitivity_summary.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        rows = list(csv.DictReader(stream))
    by_variant = {row["variant"]: row for row in rows}

    assert metrics["null_draws"] == 10_000
    assert metrics["BH_resolution_adequate"] is True
    assert set(by_variant) == {
        "column_minmax",
        "column_gene_rank",
        "uniform_ratio_gene_rank",
        "symmetric_gene_rank",
    }
    assert float(
        by_variant["uniform_ratio_gene_rank"]["gene_rho_degree"]
    ) == pytest.approx(0.5168782805267153)
    assert float(
        by_variant["symmetric_gene_rank"][
            "CDK46_q_bh_across_variants"
        ]
    ) == pytest.approx(0.034663200346666666)
