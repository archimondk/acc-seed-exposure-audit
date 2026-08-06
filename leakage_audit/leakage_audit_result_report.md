# Seed–target leakage audit: frozen-result interpretation

**Protocol:** `leakage_audit_v1`  
**Frozen protocol SHA-256:** `36c9638ade80bb761f6e8481889575b1b80feb460bbd51c9d5d68617a4155e85`  
**RNG seed:** 20260727  
**Matched-null draws:** 10,000 per arm  
**Decision status:** `PARTIAL_OR_NOT_SUPPORTED`

## 1. Primary finding

The prespecified strong leakage hypothesis was **not supported across the normalization multiverse**. In the locked primary `column_minmax` construction, removing RB1 from the ACC restart vector reduced the abemaciclib matched-null z score by 3.907, whereas adding RB1 to the breast-cancer restart vector increased it by 1.961. These effects were directionally concordant and satisfied L1–L3. However, only one of four normalization variants satisfied all three rules, below the frozen L4 requirement of at least three variants.

This is a mixed result, not a null result. It shows that direct RB1 exposure can dominate abemaciclib's score under the primary min–max scaling, but it does not establish seed–target leakage as a normalization-robust explanation of the complete propagation output.

## 2. Prespecified criteria

| Criterion | Frozen rule | Observed | Result |
|---|---:|---:|---|
| L1 | Δz_ACC ≥ 2.0 | 3.907 | Pass |
| L2 | Δz_Breast ≥ 1.5 | 1.961 | Pass |
| L3 | \|Δz_ACC − Δz_Breast\| ≤ 2.0 | 1.947 | Pass |
| L4 | L1–L3 pass in ≥3/4 variants | 1/4 | **Fail** |
| F1 | z_A2(abemaciclib) ≥ 2.0 | −1.173 | Not triggered |
| NC1 | \|Δz(ribociclib)\| ≤ 0.5 in both diseases | 0.014/0.015 | Pass |
| NC2 | \|Δz\| ≤ 0.5 for every unexposed drug | maxima 0.319/0.367 | Pass |

## 3. Normalization-specific effects

| Variant | Δz_ACC | Δz_Breast | L1 | L2 | L3 | All primary rules |
|---|---:|---:|---|---|---|---|
| column_minmax | 3.907 | 1.961 | Pass | Pass | Pass | **Pass** |
| column_gene_rank | 0.713 | 0.425 | Fail | Fail | Pass | Fail |
| uniform_ratio_gene_rank | 0.731 | 1.144 | Fail | Fail | Pass | Fail |
| symmetric_gene_rank | 0.680 | 0.405 | Fail | Fail | Pass | Fail |

The result localizes the large RB1 effect to the primary min–max score construction. Rank-transforming the propagated gene scores compressed the intervention effect below the frozen L1 and L2 thresholds. The audit therefore identifies an interaction between seed exposure and score normalization rather than a universally dominant leakage mechanism.

## 4. Weight-perturbation sensitivity

In the breast-cancer B2 arm, the frozen RB1 raw weight was the median B1 raw seed weight. The abemaciclib Δz was 1.772 at 0.5× weight, 1.961 at 1.0× and 1.965 at 1.5× under `column_minmax`. Thus the primary direction is not created by a single finely tuned RB1 weight. The response largely saturates by the median anchor, although this sensitivity does not rescue the failed cross-normalization L4 criterion.

## 5. Controls and diagnostic meaning

- F1 was not triggered: after removing RB1, abemaciclib had z = −1.173 under the primary variant. The primary signal therefore did not persist independently of RB1 exposure.
- The unexposed-drug controls passed. Across 106 drugs without the relevant RB1 exposure, the largest absolute intervention shift was 0.319 in ACC and 0.367 in breast cancer, both below 0.5.
- Ribociclib, which lacks an RB1 association in the recovered pharmacogenomic network, changed by only 0.014 in ACC and 0.015 in breast cancer.

These controls argue against a global renormalization artifact moving all drug z scores. The large primary-variant abemaciclib change is specific to direct exposure to the manipulated seed, but its magnitude depends on normalization.

## 6. Relationship to the observational audit

The interventional result must be read together with, but not substituted for, the exploratory O1–O4 analyses:

- shrinkage-weighted seed overlap correlated with C_ACC at Spearman ρ = 0.828 and exactly reproduced its Top-20 set;
- 0/62 drugs without a direct ACC-seed association reached the top quartile;
- after controlling ranked seed-overlap count, the partial Spearman correlation between C_ACC and the matched-null z score was 0.042 (P = 0.667).

O1–O3 show that the primary ACC ranking is heavily structured by seed adjacency. The frozen 2×2 intervention adds causal-design evidence that RB1 drives the primary min–max abemaciclib signal, while L4 shows that this conclusion cannot be generalized across score normalizations.

## 7. Manuscript-safe interpretation

Allowed claim:

> Direct RB1 exposure produced a large, symmetric abemaciclib shift under the primary min–max construction, while unexposed-drug controls remained stable; however, the effect failed the prespecified cross-normalization robustness criterion (1/4 variants), so seed–target leakage was not established as a normalization-invariant explanation.

Claims not supported:

- “Seed–target leakage was confirmed across the framework.”
- “The positive-control failure proves that network propagation lacks disease transportability.”
- “RWR contributes no information under every normalization.”
- “The primary ACC candidates are pharmacologically validated after leakage correction.”

## 8. Required paper changes

1. Recast the positive-control section as a design audit rather than a simple failed-control narrative.
2. Disclose the ACC/breast seed asymmetry and the presence of RB1 in the ACC restart vector.
3. Correct the rev6 statement that the breast CDK4/6 group mean was below its null under every variant; it was above the null under the primary variant only.
4. Report both layers of evidence: strong observational seed-adjacency dependence and partial, normalization-specific interventional support.
5. Preserve the frozen L4 failure in the Abstract, Results, Discussion and Conclusions.
6. Retain the control-disease analysis as a diagnostic stress test, not as evidence against established CDK4/6 pharmacology.

