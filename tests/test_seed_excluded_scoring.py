from pathlib import Path

import pytest

from analysis.seed_excluded_scoring import run_analysis


PROJECT_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def result():
    return run_analysis(PROJECT_ROOT)


def test_seed_excluded_scoring_reproduces_primary_metrics(result):
    metrics = result.metrics
    assert metrics["drug_count"] == 108
    assert metrics["seed_count"] == 45
    assert metrics["mu_0"] == pytest.approx(0.0525027124, abs=5e-11)
    assert metrics["context_spearman"] == pytest.approx(0.2955, abs=5e-5)
    assert metrics["context_top20_intersection"] == 3
    assert metrics["context_top20_jaccard"] == pytest.approx(
        0.0810810811, abs=1e-10
    )
    assert metrics["composite_spearman"] == pytest.approx(
        0.6459001256, abs=1e-10
    )
    assert metrics["composite_top20_intersection"] == 7
    assert metrics["composite_top20_jaccard"] == pytest.approx(
        0.2121212121, abs=1e-10
    )


def test_seed_excluded_scoring_reproduces_coverage_and_focal_ranks(result):
    metrics = result.metrics
    assert metrics["directly_exposed_drug_count"] == 46
    assert metrics["unexposed_drug_count"] == 62
    assert metrics["zero_nonseed_gene_count"] == 2
    assert metrics["one_nonseed_gene_count"] == 17
    assert metrics["at_least_two_nonseed_gene_count"] == 89
    assert metrics["maximum_composite_rank_shift"] == 60

    rows = {row["drug"]: row for row in result.rows}
    assert rows["Abemaciclib"]["rank_nonseed_composite"] == 24
    assert rows["Palbociclib"]["rank_nonseed_composite"] == 69
    assert rows["Ribociclib"]["rank_nonseed_composite"] == 20
    assert rows["Doxorubicin"]["nonseed_gene_count"] == 0
    assert rows["Pralatrexate"]["nonseed_gene_count"] == 0

