# C4 change report: evidence-label audit

## Decision

The legacy 14-drug ROC-AUC/PR-AUC analysis is retired. Its labels are not a coherent drug-specific clinical endpoint.

## Audit result

- Drug-specific clinical monotherapy: 3
- Clinical regimen/case (confounded): 6
- Direct ACC preclinical: 2
- Class extrapolation: 3
- Strict drug-specific clinical subset: 2 positive, 0 negative.
- Consequence: ROC-AUC and PR-AUC are not estimable; no replacement numeric discrimination claim is reported.

## Corrected interpretations

- Etoposide, doxorubicin and cisplatin inherit only regimen-level EDP-M evidence; individual effects are not identified by FIRM-ACT.
- Gemcitabine and erlotinib were evaluated largely as combinations.
- Palbociclib and ribociclib provide direct preclinical, not clinical, evidence.
- Afatinib and osimertinib were mislabeled from gefitinib class extrapolation.
- Sunitinib is reclassified as limited/mixed rather than an unambiguous negative.
- Carboplatin cannot inherit cisplatin evidence.

## Manuscript synchronization map

- Remove the AUC 0.40, bootstrap interval, permutation P value, PR-AUC, subset AUC and leave-one-label-out claims from the abstract, Methods 2.11, Results, Discussion, Conclusion and Limitations.
- Replace the old Figure 4 file and legend with the evidence-audit figure and caption in the figure manifest.
- Distinguish the blinded Claude–human categorical traceability comparison from the continuous S_external scores, which were not independently rescored.
- Report clinical regimen evidence, preclinical evidence and class extrapolation separately wherever individual drugs are discussed.

## Integrity status

- The same literature informed S_external, so this evidence set is not described as independent validation.
- All 19 source records have a primary-source/guideline locator review.
- The 12 unique DOIs resolve to matching title/year records in Crossref and OpenAlex; no scripted retraction signal was found. This signal check does not guarantee that every historical retraction database is covered.
- The blinded Claude model classification is complete: strict eligibility and strict binary labels matched the primary human curator in 19/19 records. No human rereview was performed.
- Six residual field disagreements across five records were taxonomy-only and did not change strict inclusion.

## Files

- `data/evidence/evidence_labels_v3_adjudicated.csv`
- `results/evidence_audit/second_reviewer_agreement.json`
- `figure_data/revision/Fig4_evidence_audit_primary108.csv`
- `results/evidence_audit/evidence_audit_metrics.json`
- `results/evidence_audit/source_locator_audit.md`
- `results/evidence_audit/data_card.md`
- `projects/ACC-PHARMA-NET/figures/F4_plan_card.md`
- `figures/revision/Fig4_evidence_audit_primary108.pdf`
- `figures/revision/Fig4_evidence_audit_primary108.svg`
- `figures/revision/Fig4_evidence_audit_primary108.png`
- `results/evidence_audit/C4_figure_QA.json`
- `projects/ACC-PHARMA-NET/figures/manifest.md`
- `databases/db09-projects/projects/ACC-PHARMA-NET/version_history.md`
- `results/evidence_audit/citation_verification.json`
