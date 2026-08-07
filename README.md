# Seed exposure and normalization sensitivity in network-propagation drug ranking by gene-set averaging: an adrenocortical carcinoma audit

Reproducible analysis package for a methodological audit of gene-set-averaged
network-propagation drug ranking in adrenocortical carcinoma (ACC).

## Scope and scientific status

This repository does not present a validated ACC efficacy predictor. It tests
how direct disease-seed membership, weight assignment, network topology,
normalization and drug gene-set size structure a fixed pharmacogenomic drug
ranking.

The locked primary universe contains 108 drugs with complete NCATS MIPE and
NCI-60 measurements. The primary score is:

```text
ADRS_comp = 0.50 * percentile(C_ACC) + 0.50 * percentile(ACC-relative residual)
```

The final audit reports negative and normalization-sensitive findings rather
than selecting a favourable branch. The frozen RB1 intervention met its full
rule set under one of four normalizations; the authoritative verdict is
`PARTIAL_OR_NOT_SUPPORTED`. A complete 45-seed influence scan and a
seed-excluded scoring arm are reported as result-known post-hoc sensitivity
analyses.

## Evidence-review provenance

Anthropic Claude performed a protocol-locked, blinded, **model-based** second
classification of 19 literature-evidence records. The exact Claude
model/version was not recorded and could not be recovered from the saved
prompt, output or file metadata. No independent human rereview was performed,
and this repository makes no human inter-rater reliability claim. The
continuous `S_external` scores were not independently rescored and are used
only for supplementary reprioritization, not as an independent validation
endpoint.

## Environment

- Python: 3.11-3.14; verified reference environment: Python 3.14.3
- Exact dependencies: `requirements-lock.txt`
- Package declaration: `pyproject.toml`
- Frozen input hashes: `reproducibility/input_manifest.json`

Windows PowerShell example:

```powershell
py -3.14 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements-lock.txt
```

The two large STRING v12.0 inputs are intentionally excluded from Git. Fetch
them from the official STRING download server and verify the frozen SHA-256
values with:

```powershell
python -m scripts.fetch_external_inputs --project-root .
```

## Reproduce the complete rev13 package

Run all 21 central stages, including the isolated core, positive control, six
seed-exposure intervention arms, post-hoc amendments, figures and scientific
regression tests, with:

```powershell
python -m scripts.reproduce_rev13 --project-root .
```

The final orchestration manifest is written to
`results/reproducibility/rev13_orchestration_manifest.json`. On the reference
Windows workstation, the verified run required 100.8 minutes in summed stage
time; most of that time was spent in the seven 10,000-draw null analyses.
After a host interruption, resume at a named stage recorded by `--list`, for
example:

```powershell
python -m scripts.reproduce_rev13 --project-root . --start-at leakage_run_B2
```

Resume is allowed only when every preceding stage is already recorded as
successful. The completed verification report is
`results/reproducibility/REV13_VERIFICATION_REPORT_20260805.md`.

## Reproduce the locked core and Amendments 4-5

From the repository root:

```powershell
python -m scripts.reproduce --project-root .
```

The isolated runner:

1. verifies SHA-256 for 14 frozen inputs;
2. creates a new timestamped directory under `repro_outputs/`;
3. copies only the declared inputs into the isolated workspace;
4. regenerates the locked primary analysis, exact mechanism tests, evidence
   audit, method-strengthening baselines/nulls, four normalization variants,
   shrinkage/missingness analyses, the full 45-seed leave-one-out scan and
   seed-excluded scoring;
5. compares 38 scientific outputs with the authoritative files; and
6. writes machine-readable and human-readable run manifests.

The degree-matched null uses 10,000 draws; the verified isolated core required
18.8 minutes on the reference Windows workstation. The runner refuses to overwrite
a non-empty run directory. To choose the destination:

```powershell
python -m scripts.reproduce --project-root . --run-dir repro_outputs/my_run
```

## Additional frozen analyses

The following analyses have independent, deterministic entry points and
retained protocol/freeze records:

```powershell
python -m analysis.positive_control
python -m analysis.dirichlet_weight_sensitivity
python -m analysis.seed_weight_assignment_sensitivity
python -m analysis.leave_one_seed_out
python -m analysis.seed_excluded_scoring
python -m analysis.supplementary_legacy_figures
```

The six-arm RB1 intervention uses:

```powershell
python make_arms.py
python -m analysis.leakage_audit run-arm --arm A1
python -m analysis.leakage_audit run-arm --arm A2
python -m analysis.leakage_audit run-arm --arm B1
python -m analysis.leakage_audit run-arm --arm B2
python -m analysis.leakage_audit run-arm --arm B2_lo
python -m analysis.leakage_audit run-arm --arm B2_hi
python evaluate_arms.py
```

All formal runs use the frozen random seeds, inputs and decision rules recorded
under `experiments/` and `leakage_audit/`.

## Tests

```powershell
python -m pytest -q
```

Tests cover formulas, complete-case membership, exact mechanism nulls, figure
source data, evidence-label gates, frozen-input integrity, positive-control
rules, intervention rules, normalization, component/assignment sensitivity,
the 45-seed influence scan and seed-excluded scoring.

## Project layout

- `analysis/`: deterministic analysis and figure modules.
- `scripts/`: isolated reproduction entry point.
- `tests/`: scientific regression and integrity tests.
- `experiments/`: protocols, freeze records and reproducibility checklist.
- `data/`: small frozen analytical inputs and derived data.
- `results/`: authoritative tables, metrics, manifests and audit reports.
- `leakage_audit/`: observational and six-arm intervention artifacts.
- `figure_data/`: plotted source data.
- `figures/revision/`: main and supplementary figure assets.
- `reproducibility/`: frozen input inventory and release manifests.

## Third-party data

Third-party raw files are not automatically licensed for redistribution by
this repository. Where redistribution is not explicitly permitted, the public
release will provide the official source, version, retrieval date, byte count
and SHA-256 instead of the raw file. See `DATA_SOURCES.md` and the positive-
control source manifest for provenance.

## Licences

Project-authored analysis code is released under the MIT License; see
`LICENSE`. Project-authored derived data, figure-source data, protocols and
documentation are released under CC BY 4.0 within the scope defined in
`DATA_LICENSE.md`. Neither licence applies to third-party raw data, publisher
content or upstream materials identified in `DATA_SOURCES.md`.

## Interpretation guardrails

- Pharmacogenomic associations are not confirmed molecular targets.
- Direct seed overlap is a diagnostic baseline, not biological validation.
- The CDK4/6 result is not efficacy evidence.
- Same-resource NCI-CCR/NCATS comparisons are cross-platform concordance.
- High global rank correlation does not imply exact-rank or candidate-level
  stability.
- Seed exclusion is a computational mitigation and sensitivity analysis, not
  proof that the remaining ranking is biologically superior.
- No wet-lab, prospective or independent clinical validation is claimed.

## Citation and release status

Citation metadata are provided in `CITATION.cff`. The current archival release
is [`v1.0.3`](https://github.com/archimondk/acc-seed-exposure-audit/releases/tag/v1.0.3),
and its Zenodo concept DOI is
[`10.5281/zenodo.21824511`](https://doi.org/10.5281/zenodo.21824511). MIT and
CC BY 4.0 apply to author-owned code and author-generated data/documentation,
respectively, with the third-party data restrictions described above.
