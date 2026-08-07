# Rev13 full-reproduction verification report

Verification date: 2026-08-05
Target: Pharmaceutics rev13 manuscript and supplementary package
Overall computational status: **PASSED**

## Complete orchestration

- Entry point: `python -m scripts.reproduce_rev13 --project-root .`
- Final manifest: `results/reproducibility/rev13_orchestration_manifest.json`
- Manifest schema: 2
- Successful stages: 21/21
- Non-zero return codes: 0
- Summed stage time: 6046.3 s (100.8 min)
- Initial start: 2026-08-03T15:11:42Z
- Resume start: 2026-08-05T12:12:03Z
- Completion: 2026-08-05T12:47:33Z

The first orchestration process was interrupted by the Windows host while the
B2 arm was running. This was not a scientific-analysis failure. The workflow
was resumed at `leakage_run_B2`; schema-2 resume logic retained only preceding
stages with return code 0. B2, B2_lo and B2_hi then completed successfully.

## Isolated core verification

- Runner: `rev13-reproduce-v5-amendments-4-5`
- Isolated run: `repro_outputs/run_20260803T151144Z`
- Runtime: 1130.6 s (18.8 min)
- Frozen input hashes: 14/14 passed
- Authoritative scientific outputs: 38/38 byte-identical
- Expected core figure files: 9/9 present
- C4 evidence, method-strengthening, normalization, reviewer-minor-audit and
  Amendments 4-5 gates: passed
- Scientific-output mismatches: none

## Final scientific regression tests

- Passed: 87
- Skipped: 1
- Failed: 0

The single skip is the obsolete rev5 Word-layout asset test. That asset is
intentionally excluded from the rev13 public package; its skip does not remove
any scientific regression test.

## Key reproduced audit outcomes

- The outcome-blind positive control completed as specified but did not meet
  its four prespecified recovery criteria. This is the scientific result, not
  an execution error.
- The frozen six-arm RB1 intervention verdict remains
  `PARTIAL_OR_NOT_SUPPORTED`.
- Dirichlet sensitivity reproduced 1,000 draws.
- W1/W2 weight-assignment sensitivity reproduced 1,000 permutations.
- The 45-seed leave-one-out scan reproduced 180 seed-by-normalization runs;
  minimum rank rho was 0.9167 and RB1 produced the largest focal shift.
- Seed-excluded scoring reproduced rank rho 0.6459 and retained 7/20 locked
  Top-20 drugs.
- The evidence audit reproduced 14 legacy labels, two strict positive
  candidates, no unambiguous strict negative comparator and non-estimable AUC.
- Supplementary Figure S1 reproduced 203 weight settings; Figure S2 reproduced
  rho = -0.111 with n = 103.

## Figure and supplement integrity

- Main Figures 1-5 are present; Figures 2-5 have PDF, SVG and PNG exports.
- Supplementary Figures S1-S6 each have PDF, SVG and PNG exports.
- All actual file references in the S1-S29/S1-S6 supplementary index resolve.
- Figure 4 and Figure 5 were regenerated during the completed run.

## Remaining non-computational release gates

1. A GitHub repository URL and release tag must be created.
2. The tagged release must be archived and assigned a real DOI.
3. All five authors must approve the final CRediT statement and manuscript.
4. The redistribution status of upstream B-index-derived tables must be
   confirmed before public upload.
5. The manuscript Data Availability statement must then be replaced with the
   real repository URL, tag and DOI.

No URL, DOI, upstream redistribution approval or all-author manuscript approval
is inferred by this report. MIT and CC BY 4.0 were selected on 2026-08-05 for
the author-owned materials within their recorded scopes.
