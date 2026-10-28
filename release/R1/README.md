# R1 reproducibility and provenance

Prepared 28 October 2026 for PLOS ONE manuscript PONE-D-26-43321, release v1.1.0.
The public base is commit bdd0c7030d0cb024427dfa4d7e5faa2a5e92ba23.

## Contents

- `analysis/melanoma_*.py` and `tests/test_melanoma*.py`: implementation and checks.
- `04_melanoma_protocol/`: unchanged publicly frozen stage-A/stage-B files.
- `data/melanoma_feasibility_20261018/*source_manifest.json`: provider routes and hashes.
- `results/melanoma_stageA_seeds_python314_20261019/`: local seed-run aggregate outputs.
- `results/melanoma_stageB_20261019/`: local stage-B rerun and validation outputs.
- `results/melanoma_S4_S5_propagation/`, `results/melanoma_S4_S6_local/`: additional branches.
- `results/melanoma_h5_124_supplement_20261019/`: expanded cross-disease comparison.
- `results/acc_parameter_sensitivity_R1_20261019_rerun/` and `results/acc_task5/`: ACC diagnostics.
- `independent_review/human_R1_20261026/`: locked human labels and final adjudicated R1 labels.
- `release/R1/supplied_packages/`: supplementary results supplied to the project,
  kept distinct from locally rerun outputs; their presence does not imply a new rerun.

## Reproduction

Use the repository-root environment and core reproduction instructions first.
Obtain provider raw inputs using the source-manifest URLs and verify SHA-256.
The archived run manifests record the executed arguments and frozen input hashes.
They retain original local path names as historical provenance; map those paths
to the extracted repository and your raw-data directory. The historical runtime
folder must contain `04_melanoma_protocol/` and the seed outputs specified in
the run manifest. `RUN_STAGE_B_original.md` and `RUN_S4_S6.md` retain the original
execution guidance. Use `python -m analysis.melanoma_audit --help` and the
corresponding module help for current arguments. Raw provider files and
patient/sample extract tables are regenerated locally, not redistributed here.

The two large `null_C_raw_*.npz` arrays are distributed in the separate
`melanoma_null_arrays_v1.1.0.zip` release asset. Extract it at the repository
root. `release/R1_NULL_ARRAYS.csv` records their paths, byte sizes and hashes.
`release/R1_MANIFEST.csv` checks all source-archive files except itself.

## Chronology limits

The protocol records local freezes on 18 October 2026. GitHub releases were
published at 21:56:55 EDT (stage A) and 22:44:36 EDT (stage B), equivalent to
19 October 01:56:55 and 02:44:36 UTC. The local archived rerun followed on
19 October. Public timestamps establish publication, not that all prepared
outputs were first generated afterward. Frozen historical records are preserved,
and no new prospective-validation claim is made by this release.

## Rights and review boundaries

MIT applies to author software. CC BY 4.0 applies to author-generated material
only within DATA_LICENSE.md. Third-party source terms remain separate; no
provider redistribution clearance is asserted. No publisher full texts or
signed human declaration are added. The completed declaration is held privately.
The reviewer distribution whitelist was distinct from this public code release.
