"""Frozen seed-target leakage audit and interventional-arm runner.

The observational O1-O4 analyses generated the leakage hypothesis. This
module implements only the subsequently frozen A/B arm construction,
same-engine propagation runs, hard integrity checks and prespecified verdict.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import re
import sys
import time
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import scipy

from analysis.acc_primary_pipeline import load_inputs
from analysis.method_strengthening import (
    NULL_DRAWS,
    RNG_SEED,
    load_disease_seed_weights,
)
from analysis.normalization_sensitivity import VARIANT_DESCRIPTIONS
from analysis.positive_control import (
    load_positive_control_seed_weights,
    run_positive_control_analysis,
)


PROTOCOL_ID = "leakage_audit_v1"
PROTOCOL_RELATIVE_PATH = Path("leakage_audit/leakage_audit_protocol_v1.md")
FREEZE_RELATIVE_PATH = Path("leakage_audit/FREEZE.txt")
ARMS_RELATIVE_PATH = Path("leakage_audit/arms")
PRIMARY_VARIANT = "column_minmax"
PRIMARY_DRUG = "Abemaciclib"
NEGATIVE_CONTROL_DRUG = "Ribociclib"
ARM_IDS = ("A1", "A2", "B1", "B2", "B2_lo", "B2_hi")
L1_THRESHOLD = 2.0
L2_THRESHOLD = 1.5
L3_THRESHOLD = 2.0
NC_THRESHOLD = 0.5
F1_THRESHOLD = 2.0


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"CSV has no rows: {path}")
    return rows


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: f"{value:.12g}" if isinstance(value, float) else value
                    for key, value in row.items()
                }
            )


def verify_frozen_protocol(project_root: Path) -> dict[str, Any]:
    """Verify the external, non-self-referential freeze record."""

    root = project_root.resolve()
    protocol_path = root / PROTOCOL_RELATIVE_PATH
    freeze_path = root / FREEZE_RELATIVE_PATH
    protocol = protocol_path.read_text(encoding="utf-8")
    freeze = freeze_path.read_text(encoding="utf-8")
    actual_hash = _sha256(protocol_path)
    if actual_hash not in freeze:
        raise ValueError("Frozen protocol SHA-256 does not match FREEZE.txt")
    required_protocol_tokens = (
        "**L2** | Δz_Breast | ≥ 1.5",
        r"**L3** | Effect symmetry: \|Δz_ACC − Δz_Breast\| | ≤ 2.0",
        "O4 is",
        "descriptive: no threshold",
        "No A2, B2 or sensitivity-arm result has been evaluated",
    )
    missing = [token for token in required_protocol_tokens if token not in protocol]
    if missing:
        raise ValueError(f"Frozen protocol is missing required text: {missing}")
    if "<insert" in protocol or "<hash computed" in protocol:
        raise ValueError("Frozen protocol still contains a freeze placeholder")
    if "No A2, B2, B2_lo, B2_hi" not in freeze:
        raise ValueError("FREEZE.txt does not state the interventional blind state")
    return {
        "protocol_id": PROTOCOL_ID,
        "protocol_path": PROTOCOL_RELATIVE_PATH.as_posix(),
        "protocol_sha256": actual_hash,
        "freeze_path": FREEZE_RELATIVE_PATH.as_posix(),
        "freeze_sha256": _sha256(freeze_path),
        "l1_threshold": L1_THRESHOLD,
        "l2_threshold": L2_THRESHOLD,
        "l3_threshold": L3_THRESHOLD,
        "interventional_results_observed_at_freeze": False,
    }


def build_arm_seed_sets(project_root: Path) -> dict[str, dict[str, float]]:
    """Return the six frozen reference, intervention and sensitivity seed sets."""

    root = project_root.resolve()
    a1 = load_disease_seed_weights(
        root / "data" / "ACC_P0.5C_gene_weights_v1.csv"
    )
    if len(a1) != 45 or "RB1" not in a1:
        raise ValueError("A1 must contain 45 disease-only seeds including RB1")
    a2 = {gene: weight for gene, weight in a1.items() if gene != "RB1"}

    b1 = load_positive_control_seed_weights(
        root / "results" / "positive_control" / "positive_control_seed_frozen.csv"
    )
    if len(b1) != 24 or "RB1" in b1:
        raise ValueError("B1 must contain 24 frozen seeds and exclude RB1")
    median_weight = float(np.median(np.asarray(list(b1.values()), dtype=float)))

    def add_rb1(multiplier: float) -> dict[str, float]:
        weights = dict(b1)
        weights["RB1"] = median_weight * multiplier
        return weights

    arms = {
        "A1": dict(a1),
        "A2": a2,
        "B1": dict(b1),
        "B2": add_rb1(1.0),
        "B2_lo": add_rb1(0.5),
        "B2_hi": add_rb1(1.5),
    }
    expected_sizes = {
        "A1": 45,
        "A2": 44,
        "B1": 24,
        "B2": 25,
        "B2_lo": 25,
        "B2_hi": 25,
    }
    for arm_id, weights in arms.items():
        if len(weights) != expected_sizes[arm_id]:
            raise ValueError(f"Unexpected seed count for {arm_id}")
        if any(not math.isfinite(value) or value <= 0 for value in weights.values()):
            raise ValueError(f"{arm_id} contains a non-positive seed weight")
    return arms


def _seed_rows(
    arm_id: str,
    weights: Mapping[str, float],
) -> list[dict[str, Any]]:
    total = float(sum(weights.values()))
    if total <= 0:
        raise ValueError(f"{arm_id} has no restart mass")
    return [
        {
            "arm_id": arm_id,
            "gene": gene,
            "raw_weight": float(weight),
            "normalized_weight": float(weight / total),
            "manipulated_gene": "yes" if gene == "RB1" else "no",
        }
        for gene, weight in sorted(weights.items())
    ]


def write_frozen_arm_inputs(project_root: Path) -> dict[str, Any]:
    """Write seed files and a manifest before any arm computation."""

    root = project_root.resolve()
    freeze_record = verify_frozen_protocol(root)
    arms = build_arm_seed_sets(root)
    arms_root = root / ARMS_RELATIVE_PATH
    arms_root.mkdir(parents=True, exist_ok=True)
    arm_manifest: dict[str, Any] = {
        "protocol": freeze_record,
        "rng_seed": RNG_SEED,
        "null_draws": NULL_DRAWS,
        "results_inspected_before_seed_construction": False,
        "arms": {},
    }
    for arm_id, weights in arms.items():
        arm_dir = arms_root / arm_id
        arm_dir.mkdir(parents=True, exist_ok=True)
        seed_path = arm_dir / "seeds_frozen.csv"
        _write_csv(seed_path, _seed_rows(arm_id, weights))
        arm_manifest["arms"][arm_id] = {
            "seed_n": len(weights),
            "seed_path": seed_path.relative_to(root).as_posix(),
            "seed_sha256": _sha256(seed_path),
            "raw_weight_sum": float(sum(weights.values())),
            "rb1_raw_weight": weights.get("RB1"),
            "status": (
                "reference_complete_before_freeze"
                if arm_id in {"A1", "B1"}
                else "not_executed_at_seed_freeze"
            ),
        }
    manifest_path = arms_root / "arms_manifest.json"
    manifest_path.write_text(
        json.dumps(arm_manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return arm_manifest


def _load_arm_weights(seed_path: Path) -> dict[str, float]:
    rows = _read_csv(seed_path)
    required = {"gene", "raw_weight"}
    if required - set(rows[0]):
        raise ValueError(f"{seed_path} is missing {sorted(required - set(rows[0]))}")
    weights = {row["gene"].strip(): float(row["raw_weight"]) for row in rows}
    if len(weights) != len(rows):
        raise ValueError(f"{seed_path} contains duplicate genes")
    return weights


def run_arm(
    project_root: Path,
    arm_id: str,
    n_null: int = NULL_DRAWS,
    null_batch_size: int = 64,
    output_dir: Path | None = None,
) -> dict[str, Any]:
    """Run one arm through the identical four-variant propagation engine."""

    if arm_id not in ARM_IDS:
        raise ValueError(f"Unknown arm: {arm_id}")
    root = project_root.resolve()
    freeze_record = verify_frozen_protocol(root)
    canonical_arm_dir = root / ARMS_RELATIVE_PATH / arm_id
    seed_path = canonical_arm_dir / "seeds_frozen.csv"
    arm_dir = (
        output_dir.resolve()
        if output_dir is not None
        else canonical_arm_dir
    )
    arm_dir.mkdir(parents=True, exist_ok=True)
    weights = _load_arm_weights(seed_path)
    started = time.perf_counter()
    (
        summary_rows,
        drug_rows,
        group_null_rows,
        null_matrices,
        matched_rows,
        metrics,
    ) = run_positive_control_analysis(
        root,
        seed_path,
        n_null=n_null,
        null_batch_size=null_batch_size,
        seed_weights_override=weights,
    )
    wall_seconds = time.perf_counter() - started

    output_paths = {
        "variant_summary": arm_dir / "variant_summary.csv",
        "drug_ranks": arm_dir / "drug_ranks_108.csv",
        "group_null": arm_dir / "CDK46_group_null.csv",
        "degree_matched_seeds": arm_dir / "degree_matched_seeds.csv",
        "drug_null": arm_dir / f"drug_null_{n_null}.npz",
        "metrics": arm_dir / "arm_metrics.json",
    }
    _write_csv(output_paths["variant_summary"], summary_rows)
    _write_csv(output_paths["drug_ranks"], drug_rows)
    _write_csv(output_paths["group_null"], group_null_rows)
    _write_csv(output_paths["degree_matched_seeds"], matched_rows)
    np.savez_compressed(
        output_paths["drug_null"],
        drug_names=np.asarray(metrics["primary_drug_order"], dtype="U"),
        **null_matrices,
    )
    arm_metrics = dict(metrics)
    arm_metrics.update(
        {
            "analysis_version": "seed-target-leakage-audit-v1",
            "arm_id": arm_id,
            "seed_n": len(weights),
            "seed_raw_weight_sum": float(sum(weights.values())),
            "rb1_raw_weight": weights.get("RB1"),
            "wall_clock_seconds": wall_seconds,
            "protocol_sha256": freeze_record["protocol_sha256"],
        }
    )
    output_paths["metrics"].write_text(
        json.dumps(arm_metrics, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    input_paths = (
        seed_path,
        root / PROTOCOL_RELATIVE_PATH,
        root / "data" / "bindex_network" / "bindex_edges_1304.csv",
        root / "9606.protein.info.v12.0.txt.gz",
        root / "9606.protein.links.v12.0.txt.gz",
    )
    manifest = {
        "analysis_version": "seed-target-leakage-audit-v1",
        "arm_id": arm_id,
        "command": (
            f"python -m analysis.leakage_audit run-arm --arm {arm_id} "
            f"--n-null {n_null} --batch-size {null_batch_size}"
        ),
        "rng_seed": RNG_SEED,
        "null_draws": n_null,
        "protocol_sha256": freeze_record["protocol_sha256"],
        "python": sys.version,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "logical_cpu_count": __import__("os").cpu_count(),
        "wall_clock_seconds": wall_seconds,
        "inputs": {
            path.relative_to(root).as_posix(): _sha256(path)
            for path in input_paths
        },
        "outputs": {
            path.relative_to(root).as_posix(): _sha256(path)
            for path in output_paths.values()
        },
    }
    manifest_path = arm_dir / "run_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return manifest


def evaluate_leakage_verdict(
    variant_values: Mapping[str, Mapping[str, float]],
) -> dict[str, Any]:
    """Apply the frozen L1-L4, NC1 and F1 rules without tuning."""

    missing = set(VARIANT_DESCRIPTIONS) - set(variant_values)
    if missing:
        raise ValueError(f"Missing normalization variants: {sorted(missing)}")
    per_variant: dict[str, dict[str, Any]] = {}
    for variant in VARIANT_DESCRIPTIONS:
        values = variant_values[variant]
        required = {"z_A1", "z_A2", "z_B1", "z_B2"}
        absent = required - set(values)
        if absent:
            raise ValueError(f"{variant} is missing {sorted(absent)}")
        delta_acc = float(values["z_A1"] - values["z_A2"])
        delta_breast = float(values["z_B2"] - values["z_B1"])
        l1 = delta_acc >= L1_THRESHOLD
        l2 = delta_breast >= L2_THRESHOLD
        l3 = abs(delta_acc - delta_breast) <= L3_THRESHOLD
        same_direction = delta_acc > 0 and delta_breast > 0
        per_variant[variant] = {
            "delta_z_acc": delta_acc,
            "delta_z_breast": delta_breast,
            "L1": l1,
            "L2": l2,
            "L3": l3,
            "same_direction": same_direction,
            "all_primary_rules": l1 and l2 and l3 and same_direction,
        }
    n_passing = sum(
        bool(values["all_primary_rules"]) for values in per_variant.values()
    )
    primary = per_variant[PRIMARY_VARIANT]
    l4 = n_passing >= 3
    f1 = float(variant_values[PRIMARY_VARIANT]["z_A2"]) >= F1_THRESHOLD
    nc1_acc = abs(
        float(
            variant_values[PRIMARY_VARIANT].get(
                "ribociclib_delta_acc", math.nan
            )
        )
    )
    nc1_breast = abs(
        float(
            variant_values[PRIMARY_VARIANT].get(
                "ribociclib_delta_breast", math.nan
            )
        )
    )
    nc1_available = math.isfinite(nc1_acc) and math.isfinite(nc1_breast)
    nc1_pass = (
        nc1_available
        and nc1_acc <= NC_THRESHOLD
        and nc1_breast <= NC_THRESHOLD
    )
    if f1:
        status = "FALSIFIED_ACC_SIGNAL_SURVIVES"
    elif primary["L1"] and primary["L2"] and primary["L3"] and l4:
        status = "LEAKAGE_SUPPORTED"
    else:
        status = "PARTIAL_OR_NOT_SUPPORTED"
    return {
        "status": status,
        "primary_variant": PRIMARY_VARIANT,
        "criteria": {
            "L1": {
                "passed": bool(primary["L1"]),
                "value": primary["delta_z_acc"],
                "threshold": L1_THRESHOLD,
            },
            "L2": {
                "passed": bool(primary["L2"]),
                "value": primary["delta_z_breast"],
                "threshold": L2_THRESHOLD,
            },
            "L3": {
                "passed": bool(primary["L3"]),
                "value": abs(
                    primary["delta_z_acc"] - primary["delta_z_breast"]
                ),
                "threshold": L3_THRESHOLD,
            },
            "L4": {
                "passed": l4,
                "n_variants_passing": n_passing,
                "threshold": 3,
            },
            "NC1": {
                "passed": nc1_pass if nc1_available else None,
                "delta_acc_absolute": nc1_acc if nc1_available else None,
                "delta_breast_absolute": (
                    nc1_breast if nc1_available else None
                ),
                "threshold": NC_THRESHOLD,
            },
            "F1": {
                "passed": f1,
                "z_A2": float(variant_values[PRIMARY_VARIANT]["z_A2"]),
                "threshold": F1_THRESHOLD,
            },
        },
        "per_variant": per_variant,
    }


def _drug_rows_by_key(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    rows = _read_csv(path)
    keyed = {(row["variant"], row["drug"]): row for row in rows}
    if len(keyed) != len(rows):
        raise ValueError(f"Duplicate variant-drug rows in {path}")
    return keyed


def evaluate_formal_outputs(project_root: Path) -> dict[str, Any]:
    """Validate all formal arms, calculate contrasts and write the verdict."""

    root = project_root.resolve()
    freeze_record = verify_frozen_protocol(root)
    arms_root = root / ARMS_RELATIVE_PATH
    required_arms = ("A1", "A2", "B1", "B2", "B2_lo", "B2_hi")
    arm_rows: dict[str, dict[tuple[str, str], dict[str, str]]] = {}
    checks: list[dict[str, Any]] = []
    for arm_id in required_arms:
        arm_dir = arms_root / arm_id
        seed_rows = _read_csv(arm_dir / "seeds_frozen.csv")
        matched_path = arm_dir / "degree_matched_seeds.csv"
        matched_rows = _read_csv(matched_path)
        drug_path = arm_dir / "drug_ranks_108.csv"
        drug_rows = _read_csv(drug_path)
        expected_matched = 10_000 * len(seed_rows)
        checks.extend(
            (
                {
                    "check": f"{arm_id}_matched_row_count",
                    "passed": len(matched_rows) == expected_matched,
                    "observed": len(matched_rows),
                    "expected": expected_matched,
                },
                {
                    "check": f"{arm_id}_drug_row_count",
                    "passed": len(drug_rows) == 4 * 108,
                    "observed": len(drug_rows),
                    "expected": 4 * 108,
                },
                {
                    "check": f"{arm_id}_replicate_count",
                    "passed": len(
                        {int(row["replicate"]) for row in matched_rows}
                    )
                    == 10_000,
                    "observed": len(
                        {int(row["replicate"]) for row in matched_rows}
                    ),
                    "expected": 10_000,
                },
            )
        )
        arm_rows[arm_id] = _drug_rows_by_key(drug_path)

    a2_replacements = {
        row["null_seed"]
        for row in _read_csv(arms_root / "A2" / "degree_matched_seeds.csv")
    }
    checks.append(
        {
            "check": "A2_RB1_is_replacement_candidate",
            "passed": "RB1" in a2_replacements,
            "observed": "RB1" in a2_replacements,
            "expected": True,
        }
    )
    a1_hash = _sha256(arms_root / "A1" / "degree_matched_seeds.csv")
    a2_hash = _sha256(arms_root / "A2" / "degree_matched_seeds.csv")
    checks.append(
        {
            "check": "A1_A2_matched_seed_hashes_differ",
            "passed": a1_hash != a2_hash,
            "observed": f"{a1_hash};{a2_hash}",
            "expected": "different",
        }
    )
    if not all(bool(row["passed"]) for row in checks):
        _write_csv(arms_root.parent / "verification_checks.csv", checks)
        raise ValueError("One or more formal leakage-audit checks failed")

    curve_rows: list[dict[str, Any]] = []
    primary_values: dict[str, dict[str, float]] = {}
    for variant in VARIANT_DESCRIPTIONS:
        for drug in sorted(
            {
                key[1]
                for key in arm_rows["A1"]
                if key[0] == variant
            }
        ):
            z_values = {
                arm_id: float(
                    arm_rows[arm_id][(variant, drug)]["z_degree_matched"]
                )
                for arm_id in required_arms
            }
            curve_rows.append(
                {
                    "variant": variant,
                    "drug": drug,
                    **{f"z_{arm_id}": value for arm_id, value in z_values.items()},
                    "delta_z_acc_A1_minus_A2": z_values["A1"] - z_values["A2"],
                    "delta_z_breast_B2_minus_B1": z_values["B2"] - z_values["B1"],
                    "delta_z_breast_B2_lo_minus_B1": (
                        z_values["B2_lo"] - z_values["B1"]
                    ),
                    "delta_z_breast_B2_hi_minus_B1": (
                        z_values["B2_hi"] - z_values["B1"]
                    ),
                }
            )
        primary = {
            arm_id: float(
                arm_rows[arm_id][(variant, PRIMARY_DRUG)]["z_degree_matched"]
            )
            for arm_id in required_arms
        }
        ribo = {
            arm_id: float(
                arm_rows[arm_id][
                    (variant, NEGATIVE_CONTROL_DRUG)
                ]["z_degree_matched"]
            )
            for arm_id in ("A1", "A2", "B1", "B2")
        }
        primary_values[variant] = {
            "z_A1": primary["A1"],
            "z_A2": primary["A2"],
            "z_B1": primary["B1"],
            "z_B2": primary["B2"],
            "z_B2_lo": primary["B2_lo"],
            "z_B2_hi": primary["B2_hi"],
            "ribociclib_delta_acc": ribo["A1"] - ribo["A2"],
            "ribociclib_delta_breast": ribo["B2"] - ribo["B1"],
        }
    verdict = evaluate_leakage_verdict(primary_values)
    verdict["protocol"] = freeze_record

    associations = load_inputs(root).associations
    rb1_exposed = {
        drug for drug, genes in associations.items() if "RB1" in genes
    }
    unexposed_rows = [
        row
        for row in curve_rows
        if row["variant"] == PRIMARY_VARIANT and row["drug"] not in rb1_exposed
    ]
    unexposed_acc = np.abs(
        [float(row["delta_z_acc_A1_minus_A2"]) for row in unexposed_rows]
    )
    unexposed_breast = np.abs(
        [float(row["delta_z_breast_B2_minus_B1"]) for row in unexposed_rows]
    )
    verdict["criteria"]["NC2"] = {
        "threshold": NC_THRESHOLD,
        "unexposed_drug_n": len(unexposed_rows),
        "acc_n_within_threshold": int(np.count_nonzero(unexposed_acc <= NC_THRESHOLD)),
        "breast_n_within_threshold": int(
            np.count_nonzero(unexposed_breast <= NC_THRESHOLD)
        ),
        "acc_max_absolute_delta": float(unexposed_acc.max()),
        "breast_max_absolute_delta": float(unexposed_breast.max()),
        "passed": bool(
            np.all(unexposed_acc <= NC_THRESHOLD)
            and np.all(unexposed_breast <= NC_THRESHOLD)
        ),
    }

    summary_rows: list[dict[str, Any]] = []
    for arm_id in required_arms:
        for row in _read_csv(arms_root / arm_id / "variant_summary.csv"):
            summary_rows.append({"arm_id": arm_id, **row})
    _write_csv(arms_root.parent / "leakage_curve_108.csv", curve_rows)
    _write_csv(arms_root.parent / "arm_summary.csv", summary_rows)
    _write_csv(arms_root.parent / "verification_checks.csv", checks)
    verdict_path = arms_root.parent / "verdict.json"
    verdict_path.write_text(
        json.dumps(verdict, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return verdict


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("make-arms")
    run_parser = subparsers.add_parser("run-arm")
    run_parser.add_argument("--arm", choices=ARM_IDS, required=True)
    run_parser.add_argument("--n-null", type=int, default=NULL_DRAWS)
    run_parser.add_argument("--batch-size", type=int, default=64)
    subparsers.add_parser("evaluate")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    if args.command == "make-arms":
        payload = write_frozen_arm_inputs(args.project_root)
    elif args.command == "run-arm":
        payload = run_arm(
            args.project_root,
            args.arm,
            n_null=args.n_null,
            null_batch_size=args.batch_size,
        )
    else:
        payload = evaluate_formal_outputs(args.project_root)
    print(json.dumps(payload, indent=2, ensure_ascii=False, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
