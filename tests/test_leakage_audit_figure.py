from pathlib import Path

from analysis.leakage_audit_figure import build_figure


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_build_leakage_audit_figure() -> None:
    result = build_figure(PROJECT_ROOT)
    for path in result["outputs"].values():
        assert path.is_file()
        assert path.stat().st_size > 1000
    assert result["manifest"].is_file()
    for path in result["figure_data"].values():
        assert path.is_file()

