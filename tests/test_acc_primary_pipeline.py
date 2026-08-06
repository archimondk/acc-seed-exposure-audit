from pathlib import Path

import pytest

from analysis.acc_primary_pipeline import (
    compute_primary_analysis,
    load_inputs,
    write_outputs,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_locked_primary_model_contract() -> None:
    inputs = load_inputs(PROJECT_ROOT)
    result = compute_primary_analysis(inputs)

    assert result.metrics["n_edges"] == 1304
    assert result.metrics["n_genes"] == 399
    assert result.metrics["n_drugs_all"] == 124
    assert result.metrics["n_primary"] == 108
    assert len(result.primary_rows) == 108
    assert {row["rank_comp"] for row in result.primary_rows} == set(range(1, 109))
    assert all(row["residual_pct"] is not None for row in result.primary_rows)

    for row in result.primary_rows:
        expected = 0.5 * row["C_ACC_pct"] + 0.5 * row["residual_pct"]
        assert row["ADRS_comp"] == pytest.approx(expected, abs=1e-12)


def test_primary_snapshot_from_locked_formula() -> None:
    result = compute_primary_analysis(load_inputs(PROJECT_ROOT))
    rows = {row["drug"]: row for row in result.primary_rows}

    assert "legacy_mixed_evidence_auc" not in result.metrics
    assert result.metrics["ols_slope"] == pytest.approx(0.4201796755)
    assert rows["Mitotane"]["rank_comp"] == 17
    assert rows["Doxorubicin"]["rank_comp"] == 23
    assert rows["Abemaciclib"]["rank_comp"] == 8
    assert rows["Ixazomib"]["rank_comp"] == 3


def test_outputs_are_complete_and_portable(tmp_path: Path) -> None:
    inputs = load_inputs(PROJECT_ROOT)
    result = compute_primary_analysis(inputs)
    paths = write_outputs(inputs, result, tmp_path)

    assert all(path.is_file() for path in paths.values())
    assert len(result.context_only_rows) == 16
    primary_drugs = {row["drug"] for row in result.primary_rows}
    context_only_drugs = {row["drug"] for row in result.context_only_rows}
    assert primary_drugs.isdisjoint(context_only_drugs)
    assert primary_drugs | context_only_drugs == set(inputs.associations)

    manifest = paths["manifest"].read_text(encoding="utf-8")
    assert "primary-108-v2" in manifest
    assert "/sessions/" not in manifest
