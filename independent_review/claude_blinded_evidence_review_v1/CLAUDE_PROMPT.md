# Task: independent blinded review of ACC external-evidence labels

You are the independent second reviewer for a scientific manuscript's
external-evidence audit.

Read these files first:

1. `README.md`
2. `review_protocol.md`
3. `blind_review_pairs.csv`
4. `sources/source_manifest.csv`
5. the available primary-source files in `sources/`
6. `output_template.csv`

Evaluate all 19 prespecified drug–source pairs independently. Do not infer or
try to reproduce the first reviewer's result. Do not use model ranks,
manuscript conclusions, legacy scores, or aggregate counts.

For each pair:

1. Verify whether the source actually administers, tests, or explicitly
   recommends the named drug in ACC.
2. Classify evidence domain, design, drug specificity, exposure context, and
   direction using the exact vocabulary in `review_protocol.md`.
3. Decide strict clinical eligibility by applying every listed gate.
4. Extract population/model, sample size, treatment/comparator, relevant
   numerical outcomes, and an exact source locator.
5. Record whether the judgment used full text or only an abstract.
6. Use `uncertain` or `NA` rather than guessing.
7. Give one specific exclusion reason whenever strict eligibility is `no`.

Return:

- a completed CSV preserving the exact row order and columns of
  `output_template.csv`;
- a short separate list of low-confidence or unresolved records;
- no aggregate performance analysis and no changes to the source set.

Do not search for replacement papers during the primary review. If an
additional potentially relevant source is noticed, mention it only in
`reviewer_notes` so that a separate source-gap review can be run later.

