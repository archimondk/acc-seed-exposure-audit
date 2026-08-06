# C4 evidence-label data card

## Intended use

The table supports evidence-aware interpretation and reprioritization. It is not an independent clinical validation dataset.

## Unit and scope

- Unit: one named drug.
- Scope: 19 drugs with curated ACC evidence; 14 appeared in the legacy binary benchmark.
- Primary score join: frozen 108-drug complete-case ADRS_comp ranking.

## Evidence dimensions

- Domain: clinical, preclinical, or class extrapolation.
- Design: guideline/standard, prospective phase II, randomized regimen, retrospective series, case report, in-vitro experiment, or no direct study.
- Specificity: direct drug evidence versus class-only extrapolation.
- Exposure: monotherapy/standard versus combination-regimen context.
- Direction: positive, positive preclinical, mixed/limited, regimen effect, contextual case, class extrapolation, or no direct evidence.

## Benchmark gate

- Legacy labels: 10 positive and 4 negative (n=14).
- Strict drug-specific clinical candidates: 2 positive and 0 negative.
- ROC-AUC requires both classes. Because the strict negative class is empty, ROC-AUC and PR-AUC are not estimable.

## Independence and review status

- The same literature informed S_external; the table cannot be presented as an independent validation set.
- Source/locator review was completed by the primary reviewer workflow.
- Anthropic Claude performed a protocol-locked, blinded model-based classification of all 19 records; this was not a human rereview. Strict eligibility and strict binary labels agreed in 19/19 records.
- Six residual field disagreements across five records were adjudicated as taxonomy-only; none changed strict inclusion.
- Continuous S_external scores were not independently rescored; Monte-Carlo perturbation remains a sensitivity analysis, not an independent rescoring of the continuous rubric.
