# Amendment 1 to `leakage_audit_v1`

**Date:** 2026-07-29  
**Results observed at registration:** YES — all six arms had been executed and
the frozen verdict (`PARTIAL_OR_NOT_SUPPORTED`) had been recorded.

## Rationale

Frozen rules L1 and L2 applied absolute delta-z thresholds (at least 2.0 and
1.5, respectively). Those thresholds were calibrated on the primary
column-min-max scale. Rank-percentile gene scaling changes both the magnitude
and dispersion of intervention responses, so the same absolute threshold is
not directly comparable across normalization variants. L4, which required the
magnitude rules to hold in at least three of four variants, therefore tested
recurrence of an identical absolute effect size on scales with different
dispersion.

## Post-hoc additions

Within each disease × normalization cell, we additionally report:

1. the rank of the abemaciclib response among all 108 drugs by absolute
   delta-z;
2. the one-sided empirical P value of that response against the 106 drugs
   unexposed to the manipulated RB1 association; and
3. the prespecified stability control NC2 evaluated separately in every
   normalization variant.

The outputs are `amendment1_scale_free_effect.csv`,
`amendment1_delta_z_all_drugs.csv`, and `amendment1_summary.json`.

## Decision status

The frozen rules L1–L4 and their execution are unchanged. The frozen verdict
remains `PARTIAL_OR_NOT_SUPPORTED`. Every analysis added by this amendment is
labelled post hoc. The SHA-256 in `FREEZE.txt` continues to identify the
unchanged frozen protocol; this sidecar preserves that protocol file
byte-for-byte.
