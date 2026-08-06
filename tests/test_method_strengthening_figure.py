from pathlib import Path

from analysis.method_strengthening_figure import load_figure_inputs


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_method_strengthening_figure_uses_frozen_source_tables() -> None:
    baselines, genes, null_rows, metrics = load_figure_inputs(PROJECT_ROOT)

    assert len(baselines) == 8
    assert len(genes) == 399
    assert len(null_rows) == 108
    assert metrics["null_draws"] == 10_000
    assert (
        metrics["hypothesis_decisions"]["H2_disease_context_beyond_centrality"][
            "status"
        ]
        == "partially_supported"
    )
