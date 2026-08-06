# Protocol Amendment 4: full leave-one-seed-out influence scan

## Status and provenance

- Protocol ID: `amendment4_leave_one_seed_out_v1`
- Date frozen: 2026-07-30
- Status: reviewer-requested, post-hoc descriptive analysis.
- Timing: the primary analysis, RB1 add/remove intervention, and Amendments 2–3 were already known. This protocol was written and frozen before any full 45-seed leave-one-out result was generated.
- Purpose: determine whether the influence pattern documented for RB1 is exceptional or part of a broader dependence on individual seed membership.

## Fixed inputs

- The locked 45 disease-biology-only ACC seeds and their baseline component-derived weights.
- The frozen STRING v12.0 graph and 399-gene pharmacogenomic association universe.
- The locked 108-drug complete-case universe.
- Restart probability `alpha = 0.40`, shrinkage pseudo-count `k = 3`, and locked ACC-relative potency residual.
- The four already defined propagation/normalization variants:
  1. `column_minmax`;
  2. `column_gene_rank`;
  3. `uniform_ratio_gene_rank`;
  4. `symmetric_gene_rank`.

## Intervention

For each of the 45 seeds in turn:

1. remove that seed from the restart vector;
2. renormalize the remaining 44 baseline weights to sum to one;
3. rerun each of the four propagation/normalization variants;
4. recompute `C_ACC`, its within-universe percentile, `ADRS_comp = 0.50·P(C_ACC) + 0.50·P(residual)`, and deterministic ordinal ranks (descending score, drug name as the tie breaker).

No seed is added, no coefficient is re-estimated, no outcome is used to choose a branch, and no degree-matched null is regenerated.

## Prespecified descriptive outcomes

For every omitted seed and normalization variant:

- Spearman correlation of leave-one-out versus locked `ADRS_comp` ranks;
- Top-20 Jaccard overlap with the locked Top 20;
- median and maximum absolute rank change across 108 drugs;
- number of drugs changing by at least 5 and at least 10 ranks;
- number of drugs directly associated with the omitted seed in the frozen drug–gene map;
- median and maximum absolute rank change among directly associated drugs;
- median and maximum absolute rank change among drugs not directly associated with the omitted seed.

Across seeds:

- rank RB1 by maximum absolute rank change and by the exposed-versus-unexposed contrast;
- report how many seeds have no direct drug association in the locked 108-drug universe;
- identify the most influential seeds, but make no confirmatory significance claim.

## Interpretation rules

- A large effect for an omitted seed supports sensitivity to seed membership; it does not establish that the seed is biologically invalid.
- A larger change among directly associated drugs than among unexposed drugs is consistent with direct seed exposure, but remains descriptive because seed degree, baseline weight, and network position are not randomized.
- If multiple seeds produce effects comparable with or larger than RB1, the manuscript must generalize the RB1 case to a seed-influence phenomenon.
- If RB1 is uniquely influential, the manuscript must retain the narrower RB1-specific boundary.
- Stability of the whole ranking must not be used as evidence of efficacy or external validity.

## Output and audit requirements

- Only a Markdown audit report is created in this revision cycle.
- The report must include the protocol SHA-256, software versions, input SHA-256 values, convergence checks, aggregate results, and the complete 45-seed-by-variant summary as a Markdown table.
- Any manuscript sentence using these results must label the scan as reviewer-requested and post-hoc.
