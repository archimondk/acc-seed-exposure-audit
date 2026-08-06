# Positive-control protocol for the reviewer-requested transportability test

## Status

Protocol and machine-readable selection rules frozen on 2026-07-28 before
inspection of any control-disease drug rank. Execution pending.

## Control disease and accepted mechanism

- Control disease: estrogen-receptor-positive, HER2-negative breast cancer.
- Accepted mechanism: CDK4/6 inhibition.
- Rationale: this disease-mechanism pair was suggested in the external review
  and is clinically established. It was selected before inspecting any
  control-disease ranking.

## Outcome-blind cohort and seed construction

1. Use the TCGA 2012 breast-cancer study (PMID 23000897; PMCID PMC3465532)
   and a byte-locked cBioPortal DataHub `brca_tcga` snapshot at commit
   `58341090c8bf0368ebe03f7aa95ec5137a8def25`.
2. Define ER positivity as `ER_STATUS_BY_IHC == Positive`. Define HER2
   negativity as FISH negative, or as IHC negative when FISH is unavailable
   or not evaluated. Any positive IHC or FISH result overrides a negative
   result. Exclude equivocal/indeterminate cases without a definitive negative
   FISH result. Retain one lexicographically first primary sample per patient
   and assay.
3. Mutation candidates are the genes reported as significantly mutated in
   the TCGA paper's main text. Include a candidate when a nonsynonymous
   mutation occurs in at least 10 sequenced control patients and at least 5%
   of sequenced control patients.
4. Independently select RPPA lineage/subtype features by comparing the locked
   control cohort with other primary breast tumours having definitive ER and
   HER2 classifications. Require at least 50 samples per group, a two-sided
   Mann-Whitney BH q < 0.01, and an absolute median z-score difference >= 0.5.
5. Include genes only through predefined genomic-driver, recurrence and
   lineage/biomarker fields. Mutation recurrence is capped at 1.0 after
   division by 0.20; the RPPA lineage score is capped at 1.0 after division
   of the absolute median difference by 2.0.
6. Use the same disease-only raw-weight formula as the ACC analysis:
   `0.30*genomic_driver + 0.20*recurrence + 0.20*core_pathway +
   0.10*lineage_biomarker + 0.05*prognostic_subtype`. Unused components are
   zero, and retained weights are normalized to sum to one only when building
   the restart vector.
7. Force-exclude `CDK4` and `CDK6` even if they satisfy an objective input
   rule. This prevents the accepted therapeutic targets from being copied
   directly into the disease restart vector.
8. Do not use drug screens, treatment response, drug names, pharmacologic
   targets or the expected CDK4/6 rank to select or weight seeds.
9. Exclude the therapeutic component and any explicit
   therapy-response/vulnerability extension, matching the ACC primary
   analysis.
10. Record every candidate gene, inclusion decision, source, evidence field
    and raw component score.
11. Freeze the CSV and its SHA-256 hash before running network propagation.

Required columns:

`gene, genomic_driver, recurrence, core_pathway, lineage_biomarker,
prognostic_subtype, raw_weight, include_primary, exclusion_reason, source_id,
source_version`

## Locked computational procedure

- Use the same STRING v12 graph, edge threshold, restart probability,
  convergence criterion, dangling-mass handling, 399-gene projection and
  108-drug complete-case universe as the ACC primary analysis.
- Use the same `C_ACC` shrinkage rule with `k = 3` and the same
  pharmacogenomic association network.
- Do not refit ADRS weights or tune the seed list after seeing ranks.
- Run the primary column-normalized construction and the same three
  normalization sensitivity variants.
- Run 10,000 degree-matched random-seed draws with add-one empirical P values.
- Apply BH correction across the 108 drug-level tests and separately across
  the four CDK4/6 group tests.

## Prespecified success criteria

The positive control is considered recovered only if all primary criteria
hold:

1. The CDK4/6 group has empirical P < 0.05 under the primary
   column-normalized matched-seed null.
2. The primary variant's group q, obtained by BH adjustment across the four
   prespecified normalization-variant group tests, is < 0.05.
3. At least two of abemaciclib, palbociclib and ribociclib fall in the Top 25%
   of the control-disease `C_ACC` ranking.
4. The observed CDK4/6 mean percentile exceeds the corresponding
   degree-matched null mean in at least three of four normalization variants.

Failure, partial recovery and normalization dependence will be reported
without changing the seed definition or success thresholds.

## Required execution artifacts

- Frozen seed CSV and SHA-256 manifest.
- Drug ranks and `C_ACC` values for all four normalization variants.
- Ten-thousand-draw drug-level and CDK4/6 group null tables.
- One summary JSON containing every prespecified criterion.
- A report stating pass, partial recovery or fail.
- An isolated reproduction record with software versions, platform, logical
  CPU count and wall-clock time.

## Interpretation boundary

Recovery would show that the implementation can retrieve one established
disease-mechanism relationship under a separately defined disease context.
It would not validate ACC candidates or establish clinical efficacy. Failure
would weaken claims that the pipeline transports across cancers and must be
reported as such.
