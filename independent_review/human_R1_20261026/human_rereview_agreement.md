# Independent human classification of the 19 evidence records

Locked file SHA-256: `9254f2f2fe01361a8fb28e05cdf8a1f97d878d546499d971a8e01fc11cb234ed`. Vocabulary/consistency issues: 1. Rows with unquoted commas: B07, B16 (categorical fields unaffected; free-text columns after strict_exclusion_reason may be shifted).

| Field | vs primary curator | vs final adjudicated | vs Claude (model audit) |
|---|---|---|---|
| evidence_domain | 18/19 (0.740–0.999) | 18/19 (0.740–0.999) | 18/19 (0.740–0.999) |
| evidence_design | 17/19 (0.669–0.987) | 17/19 (0.669–0.987) | 17/19 (0.669–0.987) |
| drug_specificity | 17/19 (0.669–0.987) | 17/19 (0.669–0.987) | 18/19 (0.740–0.999) |
| exposure_context | 19/19 (0.824–1.000) | 19/19 (0.824–1.000) | 19/19 (0.824–1.000) |
| direction | 12/19 (0.384–0.837) | 16/19 (0.604–0.966) | 17/19 (0.669–0.987) |
| strict_candidate_eligible | 19/19 (0.824–1.000) | 19/19 (0.824–1.000) | 19/19 (0.824–1.000) |
| strict_binary_label | 19/19 (0.824–1.000) | 19/19 (0.824–1.000) | 19/19 (0.824–1.000) |

Human strict positives: Cabozantinib, Mitotane; strict negatives: none.

Intervals: exact binomial 95%. Field-level disagreements: human_rereview_disagreements.csv.
