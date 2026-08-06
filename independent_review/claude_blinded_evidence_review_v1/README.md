# Claude blinded external-evidence review package

## Purpose

This package records a protocol-locked, blinded **model-based** classification
of 19 prespecified drug–source pairs used in the ACC external-evidence audit.
The classification was generated with Anthropic Claude; it was not an
independent human rereview. The exact Claude model/version was not recorded
and could not be recovered from the saved files.

The original prompt and locked output are retained unchanged for provenance.
The prompt's phrases “independent second reviewer” and “reviewer” are role
instructions to the model and must not be interpreted as evidence that a
second human reviewer participated.

The model task is a **closed-set, blinded classification**. The model should
judge only the supplied drug–source pairs and should not see the first
reviewer's labels, scores, source excerpts, model ranks, manuscript
interpretation, or aggregate results.

## Files

- `CLAUDE_PROMPT.md`: task prompt to paste into a new Claude conversation.
- `blind_review_pairs.csv`: the 19 prespecified drug–source review units.
- `review_protocol.md`: controlled vocabulary and decision rules.
- `output_template.csv`: prefilled row identifiers and empty review fields.
- `sources/source_manifest.csv`: 12 unique sources and local-access status.
- `sources/missing_full_text.csv`: sources that still need legally obtained
  full text before the strongest possible review.
- `sources/pubmed_abstracts.xml`: official PubMed records for sources with a
  PMID.
- `sources/*.pdf`: publicly retrievable primary-source full texts.
- `sources/README.md`: instructions for unavailable full texts.

## Blinding rules

Do not add any of the following project files to this package:

- `data/evidence/evidence_labels_v2.csv`
- `results/evidence_audit/source_locator_audit.md`
- `results/evidence_audit/evidence_audit_metrics.json`
- `results/evidence_audit/C4_change_report.md`
- Figure 4 or its caption
- the manuscript Discussion
- ADRS scores, ranks, or legacy evidence scores

The historical workflow used a new conversation without prior project
context. Claude's first completed output was preserved unchanged before
comparison with the primary human curator.

## Recommended workflow

1. Add any legally obtained missing full texts listed in
   `sources/source_manifest.csv`.
2. Upload this directory to a new Claude conversation.
3. Paste `CLAUDE_PROMPT.md`.
4. Ask Claude to return a completed copy of `output_template.csv`.
5. Save the unmodified result as `claude_review_locked.csv`.
6. Compare the model output with the primary human labels only after the model
   output file is locked; report this as model–human traceability, not human
   inter-rater reliability.
