import json
from pathlib import Path

from scripts.reproduce_rev13 import _prepare_manifest, build_stages


def test_rev13_orchestrator_covers_all_central_analysis_families() -> None:
    stages = build_stages(Path("."), python="python")
    names = {stage.name for stage in stages}
    assert {
        "core_and_amendments_4_5",
        "positive_control",
        "leakage_make_arms",
        "leakage_evaluate",
        "leakage_figure",
        "dirichlet_weights",
        "seed_weight_assignment",
        "leave_one_seed_out",
        "seed_excluded_scoring",
        "supplementary_figures_s1_s2",
        "scientific_regression_tests",
    }.issubset(names)
    leakage_runs = [name for name in names if name.startswith("leakage_run_")]
    assert len(leakage_runs) == 6


def test_resume_preserves_all_preceding_successful_stages(tmp_path: Path) -> None:
    stages = build_stages(tmp_path, python="python")
    start_name = "leakage_run_B2"
    start_index = [stage.name for stage in stages].index(start_name)
    manifest_path = tmp_path / "manifest.json"
    prior_records = [
        {"name": stage.name, "command": list(stage.command), "returncode": 0}
        for stage in stages[:start_index]
    ]
    prior_records.append(
        {"name": start_name, "command": [], "returncode": 1073807364}
    )
    manifest_path.write_text(
        json.dumps(
            {
                "started_utc": "2026-08-03T00:00:00+00:00",
                "stages": prior_records,
            }
        ),
        encoding="utf-8",
    )
    resumed = _prepare_manifest(tmp_path, manifest_path, stages, start_name)
    assert resumed["schema_version"] == 2
    assert [record["name"] for record in resumed["stages"]] == [
        stage.name for stage in stages[:start_index]
    ]
