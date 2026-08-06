from pathlib import Path

from analysis.supplementary_legacy_figures import (
    build_weight_grid,
    compute_lomco,
    generate_outputs,
)


ROOT = Path(__file__).resolve().parents[1]


def test_four_component_grid_has_203_prespecified_settings() -> None:
    grid = build_weight_grid()
    assert len(grid) == 203
    assert all(abs(sum(weights) - 1.0) < 1e-12 for weights in grid)
    assert all(weights[0] > 0 or weights[1] > 0 for weights in grid)


def test_lomco_reproduces_legacy_negative_result() -> None:
    rows, metrics = compute_lomco(ROOT)
    assert len(rows) >= 100
    assert round(metrics["spearman_rho"], 3) == -0.111


def test_supplementary_assets_are_generated(tmp_path: Path) -> None:
    outputs = generate_outputs(ROOT, output_root=tmp_path)
    assert outputs["weight_grid_metrics"]["n_weight_settings"] == 203
    assert round(outputs["lomco_metrics"]["spearman_rho"], 3) == -0.111
    for path in outputs["paths"]:
        assert path.is_file()
        assert path.stat().st_size > 0
