# Archived independent human rereview protocol for the 19 external-evidence records

> **Status update, 1 August 2026:** This optional protocol was not executed.
> The author team decided not to pursue an independent human rereview. It is
> retained for provenance only; the manuscript makes no human inter-rater
> reliability claim.

## Purpose

This protocol closes the remaining reviewer concern that the current secondary classification was performed by Anthropic Claude rather than an independent human. It must be completed by a domain-qualified human who did not perform the primary curation.

## Blinding

Provide only:

- `independent_review/claude_blinded_evidence_review_v1/blind_review_pairs.csv`;
- `independent_review/claude_blinded_evidence_review_v1/review_protocol.md`;
- `independent_review/claude_blinded_evidence_review_v1/output_template.csv`;
- the legally shareable source package in `independent_review/claude_blinded_evidence_review_v1/sources/`;
- the B02 PubMed record in `independent_review/cabozantinib_b02_rerate/`.

Do not provide:

- primary-curator labels or continuous S_external scores;
- Claude outputs or adjudication notes;
- manuscript Results or Discussion;
- drug scores, ranks, figures or aggregate agreement counts.

## Reviewer eligibility and declaration

Before rating, the human reviewer must record:

- full name and affiliation;
- relevant expertise;
- confirmation of no involvement in the primary curation;
- conflicts of interest, if any;
- date the source package was received;
- date the completed file was locked.

Authorship is not automatic. Any authorship decision must follow the journal's authorship policy and the actual contribution.

## Required output

Return the completed 19-row template without changing row order or controlled vocabulary. Also provide a signed or emailed declaration containing:

> I classified the 19 prespecified drug–source records using only the supplied protocol and sources. I did not inspect the primary labels, continuous scores, Claude classifications, manuscript interpretation, drug ranks or aggregate results before locking my completed review.

Save the first completed file unchanged as `human_review_locked.csv`. Record its SHA-256 before any comparison.

## Comparison plan

After locking:

1. compare the human review with the primary human curator on strict eligibility and strict binary label;
2. report raw agreement and an exact binomial 95% confidence interval;
3. list field-level disagreements without collapsing taxonomy categories;
4. adjudicate with both human reviewers while preserving the locked originals;
5. compare the human result with Claude only as a secondary model-audit analysis;
6. update Methods, Results, Supplementary Tables S9/S17, AI disclosure and the release manifest.

No Cohen's kappa will be reported because the strict endpoint has sparse and structurally unbalanced categories.

## Current status

Never executed. On 1 August 2026, the author team decided not to pursue this optional rereview. The current manuscript describes the existing second classification as model-based and makes no claim of independent human inter-rater reliability.
