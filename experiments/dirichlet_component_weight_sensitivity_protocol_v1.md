# ACC-PHARMA-NET Protocol Amendment 2

## Post-hoc Dirichlet sensitivity of disease-component weights

**Protocol ID:** `dirichlet_component_weight_sensitivity_v1`  
**Registration date:** 2026-07-29  
**Status at registration:** Post hoc. The primary ACC ranking, all four
normalization analyses, the positive-control analysis, the six-arm RB1
intervention, Protocol Amendment 1 and the frozen
`PARTIAL_OR_NOT_SUPPORTED` verdict had already been observed.

## 1. Objective

Quantify how uncertainty in the expert-assigned disease-component coefficients
propagates through the primary ACC-PHARMA-NET estimator and the locked
108-drug external-evidence-free ranking.

This is a descriptive stress test. It is not an efficacy validation, a
prospective hypothesis test or a basis for revising the frozen leakage-audit
verdict.

## 2. Fixed inputs and analysis universe

- The drug–gene association network, STRING v12 files, MIPE/NCI-60 residual,
  restart probability (`alpha = 0.40`), shrinkage (`k = 3`) and locked
  108-drug complete-case universe are unchanged.
- The 45 disease-biology-only ACC seeds are fixed.
- `MGMT`, `SLFN11`, `ABCB1`, `SOAT1` and `UBA1` remain excluded.
- The therapeutic-relevance component `T` remains excluded. Introducing `T`
  would alter both the estimand and, under the therapy-informed construction,
  seed membership; it is therefore outside this weight-only sensitivity.
- Only the five active component coefficients are varied:
  genomic driver (`G`), recurrence (`R`), core pathway (`P`),
  lineage/biomarker (`L`) and prognostic/subtype (`S`).

The locked primary coefficient proportions, after normalization to the
five-component simplex, are:

| Component | Proportion |
|---|---:|
| G | 0.3529411764705882 |
| R | 0.2352941176470588 |
| P | 0.2352941176470588 |
| L | 0.1176470588235294 |
| S | 0.0588235294117647 |

## 3. Perturbation distribution

- Draw count: `B = 1000`.
- Random-number generator: NumPy `default_rng`.
- RNG seed: `20260729`.
- Distribution: `Dirichlet(1, 1, 1, 1, 1)`.
- Each draw is an unconstrained global simplex stress test, not a probability
  model for clinically plausible expert uncertainty.
- No draw will be removed because it is extreme, produces an unfavorable
  ranking or changes a named candidate.

For draw `b`, the raw weight of retained seed gene `g` is

`w_g,b = theta_G,b G_g + theta_R,b R_g + theta_P,b P_g
         + theta_L,b L_g + theta_S,b S_g`.

The 45 positive gene weights are normalized to sum to one before forming the
restart vector.

## 4. Locked computational pathway

Each draw follows the primary pathway only:

`component coefficients -> 45 seed weights -> column-normalized STRING RWR
-> min-max r_ACC over the 399 associated genes -> shrinkage C_ACC
-> C_ACC percentile in the locked 108 drugs
-> ADRS_comp = 0.50*C_ACC percentile + 0.50*locked residual percentile`.

The MIPE/NCI-60 residual percentile and the 0.50/0.50 ADRS component weights
are held fixed. Normalization sensitivity is not crossed with this experiment
because the four-estimator multiverse has already been reported separately.
Degree-matched nulls are not regenerated: this experiment estimates ranking
sensitivity to component weights and makes no empirical-P or FDR claim.

Ties are broken deterministically by ascending drug name after descending
ADRS score, matching the locked primary-ranking convention.

## 5. Prespecified summaries

### Primary summary

- Distribution across the 1000 draws of Spearman correlation between each
  ADRS rank vector and the locked primary ADRS rank vector:
  minimum, 5th percentile, median, 95th percentile and maximum.

### Secondary summaries

- Top-20 Jaccard overlap with the locked primary Top 20:
  minimum, 5th percentile, median, 95th percentile and maximum.
- For every drug: locked rank, median draw rank, 5th/25th/75th/95th rank
  percentiles, rank standard deviation, and probabilities of Top 10 and
  Top 20 membership.
- For abemaciclib, palbociclib and ribociclib: the same per-drug summaries.
- Distribution of the mean rank of the three CDK4/6 drugs.
- Counts of drugs with Top-20 membership probability at least 0.80 and at
  most 0.20.
- No binary robustness threshold and no null-hypothesis P value are defined.

## 6. Quality-control gates

The run is valid only if all of the following hold:

1. Exactly 1000 coefficient vectors are generated.
2. Every coefficient vector is finite, strictly positive and sums to one
   within `1e-12`.
3. All draws retain exactly the same 45 seed genes.
4. Every restart vector is finite, strictly positive on those seeds and sums
   to one within `1e-12`.
5. Every RWR column conserves probability and converges at L1 tolerance
   `< 1e-10` within 500 iterations.
6. The baseline five-component reconstruction reproduces the frozen r_ACC
   vector with maximum absolute difference `<= 1.1e-6` and Spearman
   correlation `>= 0.999999`.
7. Every draw yields a complete permutation of ranks 1–108.
8. Output row counts, input SHA-256 hashes, software versions, logical CPU
   count and wall-clock time are recorded.

## 7. Frozen outputs

- `results/dirichlet_weight_sensitivity/component_weight_draws.csv`
- `results/dirichlet_weight_sensitivity/drug_rank_draws.csv`
- `results/dirichlet_weight_sensitivity/drug_rank_summary.csv`
- `results/dirichlet_weight_sensitivity/draw_summary.csv`
- `results/dirichlet_weight_sensitivity/dirichlet_weight_sensitivity_summary.json`
- `results/dirichlet_weight_sensitivity/run_manifest.md`
- `figures/revision/FigS3_dirichlet_weight_sensitivity.{png,pdf,svg}`

The manuscript will report all prespecified summaries regardless of direction.

