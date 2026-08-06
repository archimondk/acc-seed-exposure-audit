# ACC-PHARMA-NET Protocol Amendment 3

## Post-hoc sensitivity to seed-weight magnitude and gene assignment

**Protocol ID:** `seed_weight_assignment_sensitivity_v1`  
**Registration date:** 2026-07-29  
**Status at registration:** Post hoc. The primary ACC ranking, all four
normalization analyses, the positive-control analysis, the six-arm RB1
intervention, Protocol Amendments 1 and 2, the frozen
`PARTIAL_OR_NOT_SUPPORTED` verdict, and the hypothesis that the curated
weights may add little beyond seed membership had already been observed or
formulated. No W1 or W2 output had been generated or inspected.

## 1. Objective

Test whether the magnitude and gene assignment of the 45 fixed ACC restart
weights carry appreciable information for the locked primary 108-drug
ranking beyond binary seed membership.

This is a descriptive, result-known stress test. It is not an efficacy
validation, a prospective hypothesis test, or a basis for revising the frozen
leakage-audit verdict.

## 2. Weight definition and fixed inputs

The baseline restart weights are the weights actually used by the locked
disease-only primary model:

`w_g = 0.30 G_g + 0.20 R_g + 0.20 P_g + 0.10 L_g + 0.05 S_g`.

They exclude the therapeutic-relevance component `T` and exclude `MGMT`,
`SLFN11`, `ABCB1`, `SOAT1`, and `UBA1`. Among the 45 retained seeds, the
unnormalized five-component weights range from 0.260 to 0.770. After
normalization to a restart simplex, they range from 0.0114085 to 0.0337867
(1.14%–3.38% of restart mass), with a maximum/minimum ratio of 2.9615.

The drug–gene association network, STRING v12 graph, 45 seed genes, restart
probability (`alpha = 0.40`), shrinkage (`k = 3`), locked MIPE/NCI-60 residual,
ADRS 0.50/0.50 mixing weights, and locked 108-drug complete-case universe are
unchanged. Only the restart weights assigned to the 45 fixed seeds are
altered.

The CSV field `ACC_weight` is not used in this amendment because it includes
the excluded therapeutic component and is not the weight vector used by the
locked disease-only primary analysis.

## 3. Descriptive arms

### W1: uniform membership-only restart

Every retained seed receives restart weight `1/45`.

This arm asks whether retaining seed membership while discarding all curated
weight magnitudes materially changes the locked ranking.

### W2: permuted weight-to-gene assignment

- Permutation count: `B = 1000`.
- Random-number generator: NumPy `default_rng`.
- RNG seed: `20260729`.
- For each draw, the 45 normalized baseline restart weights are randomly
  permuted across the same 45 genes.
- Each draw preserves the weight distribution exactly and breaks only the
  curated pairing between gene identity and weight.
- Permutations are sampled independently; no draw is removed because it is
  close to the identity, extreme, or unfavorable to a named drug.

This arm asks whether the curated assignment of the observed weight
distribution to specific genes materially changes the locked ranking.

## 4. Locked computational pathway

Both arms follow the same primary pathway used in Protocol Amendment 2:

`45 fixed seed genes and arm-specific restart weights
-> column-normalized STRING RWR
-> min–max r_ACC over the 399 associated genes
-> shrinkage C_ACC
-> C_ACC percentile in the locked 108 drugs
-> ADRS_comp = 0.50*C_ACC percentile + 0.50*locked residual percentile`.

Ties are broken deterministically by ascending drug name after descending
ADRS score. Degree-matched nulls are not regenerated because the estimand is
ranking sensitivity to weight magnitude and assignment; no empirical P value
or FDR claim is made.

## 5. Prespecified summaries

### W1

- Spearman correlation between W1 and locked ADRS ranks.
- Top-20 Jaccard overlap between W1 and the locked Top 20.
- Mean and maximum absolute drug-rank change.
- Locked and W1 ranks for all 108 drugs and for the three CDK4/6 inhibitors.

### W2

- Distribution of Spearman correlation with the locked ADRS rank:
  minimum, 5th percentile, median, 95th percentile, and maximum.
- Distribution of Top-20 Jaccard overlap with the locked Top 20 using the same
  five summaries.
- The same two distributions relative to W1.
- Distribution of the mean and maximum absolute drug-rank change from the
  locked ranking.
- For every drug: locked rank, W1 rank, W2 median and 5th/25th/75th/95th rank
  percentiles, rank standard deviation, and probabilities of Top-10 and
  Top-20 membership.
- The same per-drug summaries for abemaciclib, palbociclib, and ribociclib.
- Descriptive comparison of the W2 Spearman and Top-20 Jaccard 5th–95th
  percentile intervals with the corresponding Protocol Amendment 2
  Dirichlet intervals.

No binary equivalence margin, null-hypothesis P value, or FDR procedure is
defined.

## 6. Interpretation rule

If W1 remains highly concordant with the locked ranking and the W2
concordance distributions are concentrated near one, the result will support
the bounded claim that, within this fixed seed set and fixed primary
estimator, curated weight magnitude and weight-to-gene assignment add little
ranking information beyond seed membership.

The phrase "mathematically equivalent" will not be used unless all 108 ranks
are exactly identical. The phrase "statistically undetectable" will not be
used because no equivalence margin or inferential test was prespecified.
Regardless of direction, the result cannot validate the seed set or drug
efficacy and cannot revise the frozen leakage verdict.

## 7. Quality-control gates

The run is valid only if all of the following hold:

1. W1 assigns exactly `1/45` to each of the 45 retained seeds.
2. Exactly 1000 W2 permutations are generated.
3. Every W2 column is an exact permutation of the normalized baseline weight
   vector and sums to one within `1e-12`.
4. Every arm retains exactly the same 45 seed genes.
5. Every RWR column conserves probability and converges at L1 tolerance
   `< 1e-10` within 500 iterations.
6. The unpermuted five-component reconstruction reproduces frozen r_ACC with
   maximum absolute difference `<= 1.1e-6` and Spearman correlation
   `>= 0.999999`.
7. Every arm yields a complete permutation of ranks 1–108.
8. Output row counts, input SHA-256 hashes, software versions, logical CPU
   count, and wall-clock time are recorded.

## 8. Frozen outputs

- `results/seed_weight_assignment_sensitivity/uniform_drug_ranks.csv`
- `results/seed_weight_assignment_sensitivity/permutation_draw_summary.csv`
- `results/seed_weight_assignment_sensitivity/permutation_drug_rank_draws.csv`
- `results/seed_weight_assignment_sensitivity/permutation_drug_rank_summary.csv`
- `results/seed_weight_assignment_sensitivity/seed_weight_assignment_summary.json`
- `results/seed_weight_assignment_sensitivity/run_manifest.md`
- `figures/revision/FigS4_seed_weight_assignment_sensitivity.{png,pdf,svg}`

The manuscript will report all prespecified summaries regardless of direction.
