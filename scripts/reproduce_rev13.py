"""Run every central rev13 analysis family from one documented command."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence


ARM_IDS = ("A1", "A2", "B1", "B2", "B2_lo", "B2_hi")


@dataclass(frozen=True)
class Stage:
    name: str
    command: tuple[str, ...]


def build_stages(
    project_root: Path,
    python: str | None = None,
    core_run_dir: Path | None = None,
) -> list[Stage]:
    root = project_root.resolve()
    executable = python or sys.executable
    root_arg = str(root)
    core_command = [
        executable,
        "-m",
        "scripts.reproduce",
        "--project-root",
        root_arg,
    ]
    if core_run_dir is not None:
        core_command.extend(["--run-dir", str(core_run_dir.resolve())])

    stages = [
        Stage("fetch_string_inputs", (executable, "-m", "scripts.fetch_external_inputs", "--project-root", root_arg)),
        Stage("core_and_amendments_4_5", tuple(core_command)),
        Stage("positive_control", (executable, "-m", "analysis.positive_control", "--project-root", root_arg)),
        Stage("leakage_make_arms", (executable, "-m", "analysis.leakage_audit", "--project-root", root_arg, "make-arms")),
    ]
    stages.extend(
        Stage(
            f"leakage_run_{arm}",
            (
                executable,
                "-m",
                "analysis.leakage_audit",
                "--project-root",
                root_arg,
                "run-arm",
                "--arm",
                arm,
            ),
        )
        for arm in ARM_IDS
    )
    stages.extend(
        [
            Stage("leakage_evaluate", (executable, "-m", "analysis.leakage_audit", "--project-root", root_arg, "evaluate")),
            Stage("leakage_figure", (executable, "-m", "analysis.leakage_audit_figure", "--project-root", root_arg)),
            Stage("amendment1_scale_free", (executable, str(root / "variant_scale_free_effect_v2.py"), "--project-root", root_arg)),
            Stage("network_resolution", (executable, "-m", "analysis.network_resolution_audit", "--project-root", root_arg)),
            Stage("dirichlet_weights", (executable, "-m", "analysis.dirichlet_weight_sensitivity", "--project-root", root_arg)),
            Stage("seed_weight_assignment", (executable, "-m", "analysis.seed_weight_assignment_sensitivity", "--project-root", root_arg)),
            Stage("leave_one_seed_out", (executable, "-m", "analysis.leave_one_seed_out", "--project-root", root_arg)),
            Stage("seed_excluded_scoring", (executable, "-m", "analysis.seed_excluded_scoring", "--project-root", root_arg)),
            Stage("evidence_audit", (executable, "-m", "analysis.evidence_label_audit", "--project-root", root_arg)),
            Stage("supplementary_figures_s1_s2", (executable, "-m", "analysis.supplementary_legacy_figures", "--project-root", root_arg)),
            Stage("scientific_regression_tests", (executable, "-m", "pytest", "-q")),
        ]
    )
    return stages


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--core-run-dir", type=Path)
    parser.add_argument("--list", action="store_true", help="List stages without running them")
    parser.add_argument("--start-at", help="Resume at the named stage")
    parser.add_argument("--stop-after", help="Stop after the named stage")
    return parser.parse_args(argv)


def _select_stages(stages: list[Stage], start_at: str | None, stop_after: str | None) -> list[Stage]:
    names = [stage.name for stage in stages]
    start = names.index(start_at) if start_at else 0
    stop = names.index(stop_after) + 1 if stop_after else len(stages)
    if start >= stop:
        raise ValueError("--start-at must precede or equal --stop-after")
    return stages[start:stop]


def _prepare_manifest(
    root: Path,
    manifest_path: Path,
    stages: list[Stage],
    start_at: str | None,
) -> dict[str, object]:
    now = datetime.now(timezone.utc).isoformat()
    if not start_at:
        return {
            "schema_version": 2,
            "status": "running",
            "started_utc": now,
            "project_root": str(root),
            "stages": [],
        }

    if not manifest_path.is_file():
        raise ValueError("Cannot resume: orchestration manifest does not exist")
    previous = json.loads(manifest_path.read_text(encoding="utf-8"))
    start_index = [stage.name for stage in stages].index(start_at)
    required_names = [stage.name for stage in stages[:start_index]]
    successful = {
        record.get("name"): record
        for record in previous.get("stages", [])
        if record.get("returncode") == 0
    }
    missing = [name for name in required_names if name not in successful]
    if missing:
        raise ValueError(
            "Cannot resume because preceding stages are not recorded as "
            f"successful: {missing}"
        )
    return {
        "schema_version": 2,
        "status": "running",
        "started_utc": previous.get("started_utc", now),
        "resumed_utc": now,
        "project_root": str(root),
        "stages": [successful[name] for name in required_names],
    }


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    root = args.project_root.resolve()
    stages = build_stages(root, core_run_dir=args.core_run_dir)
    if args.list:
        for index, stage in enumerate(stages, start=1):
            print(f"{index:02d}. {stage.name}: {' '.join(stage.command)}")
        return 0
    selected = _select_stages(stages, args.start_at, args.stop_after)
    manifest_path = root / "results/reproducibility/rev13_orchestration_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest = _prepare_manifest(root, manifest_path, stages, args.start_at)
    for index, stage in enumerate(selected, start=1):
        print(f"\n[{index}/{len(selected)}] {stage.name}", flush=True)
        started = time.perf_counter()
        completed = subprocess.run(stage.command, cwd=root, check=False)
        record = {
            **asdict(stage),
            "returncode": completed.returncode,
            "duration_seconds": time.perf_counter() - started,
        }
        manifest["stages"].append(record)  # type: ignore[union-attr]
        if completed.returncode != 0:
            manifest["status"] = "failed"
            manifest["failed_stage"] = stage.name
            manifest["completed_utc"] = datetime.now(timezone.utc).isoformat()
            manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
            return completed.returncode
        manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    manifest["status"] = "passed"
    manifest["completed_utc"] = datetime.now(timezone.utc).isoformat()
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"\nRev13 reproduction passed. Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
