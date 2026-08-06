from pathlib import Path

import pytest

from analysis.leave_one_seed_out import run_analysis


PROJECT_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def result():
    return run_analysis(PROJECT_ROOT)


def test_leave_one_seed_out_reproduces_registered_design(result):
    assert len(result.rows) == 45 * 4
    assert result.metrics["seed_count"] == 45
    assert result.metrics["drug_count"] == 108
    assert result.metrics["variant_count"] == 4
    assert result.metrics["run_count"] == 180


def test_leave_one_seed_out_reproduces_frozen_headline(result):
    assert result.metrics["minimum_rho"] == pytest.approx(0.9167, abs=5e-5)
    assert result.metrics["minimum_top20_jaccard"] == pytest.approx(
        0.7391304348, abs=1e-10
    )
    assert result.metrics["zero_exposure_seed_count"] == 28
    assert result.metrics["rb1_max_abs_shift"] == 57
    assert result.metrics["rb1_worst_shift_rank"] == 1
    assert result.metrics["rb1_exposed_shift_rank"] == 1


def test_leave_one_seed_out_reproduces_variant_extrema(result):
    by_variant = result.variant_summary
    assert by_variant["column_minmax"]["largest_abs_shift"] == 57
    assert by_variant["column_minmax"]["minimum_rho"] == pytest.approx(
        0.9167, abs=5e-5
    )
    assert by_variant["column_gene_rank"]["largest_abs_shift"] == 19
    assert by_variant["uniform_ratio_gene_rank"]["largest_abs_shift"] == 25
    assert by_variant["symmetric_gene_rank"]["largest_abs_shift"] == 20

