# Protocol Amendment 5: seed-excluded drug scoring

Status: result-known post-hoc; frozen before Amendment 5 outputs  
Date frozen: 30 July 2026  
Purpose: evaluate a direct mitigation for seed-exposure sensitivity requested
in the Pharmaceutics rev12 referee report

## Known information at registration

The locked primary ranking, replacement-baseline audit, control-disease
comparison, RB1 intervention, Amendments 1–4 and all associated manuscript
interpretations were known. In particular, direct weighted seed overlap
reproduced the locked primary Top 20, and the ACC and breast-cancer arms had
different direct seed-association coverage. No seed-excluded score or rank was
inspected before this protocol was frozen.

## Locked inputs

- `data/bindex_network/bindex_edges_1304.csv`
- `data/bindex_network/rACC_399_fullSTRING.csv`
- `data/ACC_P0.5C_gene_weights_v1.csv`
- `results/primary_analysis/ADRS_comp_primary_108.csv`

The disease-biology-only seed set is the same 45-gene set used by the primary
analysis: genes with a positive disease-only weight after excluding the
therapeutic component and MGMT, SLFN11, ABCB1, SOAT1 and UBA1.

## Estimand and scoring rule

For drug d with locked associated-gene set A_d and disease seed set S, define
the non-seed association set:

`A_d,nonseed = A_d \ S`.

The primary C_ACC pseudo-count and reference mean are retained. Let k = 3 and
let mu_0 be the locked association-weighted mean of r_ACC across all 1304
drug–gene edges. The seed-excluded context score is:

`C_ACC,nonseed(d) = [sum_{g in A_d,nonseed} r_ACC(g) + k*mu_0] /
                    [|A_d,nonseed| + k]`.

This keeps every drug in the locked universe. If a drug has no non-seed
associated gene, the formula returns the neutral prior mu_0 rather than
deleting the drug or assigning an arbitrary extreme value.

Rank `C_ACC,nonseed` within the same 108 complete-case drugs using
average-rank percentiles. Retain the locked ACC-relative residual percentile
and define the remedial composite:

`ADRS_nonseed = 0.50*P(C_ACC,nonseed) + 0.50*P(locked residual)`.

Deterministic ordinal ranks are descending by score with drug name as the
tie-breaker, matching the primary pipeline.

## Prespecified outputs

Primary:

1. Spearman correlation between `ADRS_nonseed` and locked `ADRS_comp`;
2. Top-20 Jaccard overlap and the entering/leaving drug sets.

Secondary:

3. Spearman correlation and Top-20 Jaccard overlap between
   `C_ACC,nonseed` and locked `C_ACC`;
4. number of drugs with zero, one and at least two retained non-seed genes;
5. maximum absolute rank shift and selected shifts for abemaciclib,
   palbociclib and ribociclib;
6. seed-exposed versus unexposed summaries.

## Interpretation boundary

The analysis is a deterministic sensitivity/remedy arm, not an independent
validation. A lower overlap would show that direct seed contributions can be
removed computationally, but would not establish that the remaining score is
more biologically valid. A high overlap would show that non-seed propagation
and the locked activity component retain much of the ordering. Either result
must be reported.
