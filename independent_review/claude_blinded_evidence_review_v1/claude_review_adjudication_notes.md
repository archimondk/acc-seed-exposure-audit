# Blinded Claude model classification — completion note and adjudication list

**Classifier role:** Anthropic Claude acting under a blinded role prompt; not a human rereviewer.
**Input:** `blind_review_pairs.csv` (19 drug–source pairs), sources in `sources/`.
**Output:** `claude_review_locked.csv` (19 rows, original order preserved, one row per input record).
**Blinding statement:** no network scores, drug ranks, manuscript text, or primary-curator labels were supplied; no target number of positive/negative labels was pursued. Classifications rest only on the supplied sources.
**Version limitation:** the exact Claude model/version was not recorded and could not be recovered from the saved files.

## Access reality
- **Full text available (4 sources):** S05 (CDK4/6 in vitro), S06 (lenvatinib+pembrolizumab), S07 (1978 infant case), S09 (cixutumumab+temsirolimus).
- **Official abstract only (7 sources):** S01, S03, S04, S08, S10, S12 (+S05 abstract).
- **Not available at all:** S02 (cabozantinib) — no full text and no abstract record in the package; S11 (gefitinib conference abstract) — abstract body not supplied.

## Outcome summary
- **Strict-eligible records: 1/19** — B01 (mitotane, ESE-ENSAT guideline) → `positive`.
- **All other 18 records: `strict_binary_label = NA`.**

Exclusion reasons (each record carries exactly one primary reason):

| Reason | Records | n |
|---|---|---|
| `combination_regimen_not_decomposable` | B03, B04, B05, B11, B14, B19 | 6 |
| `named_drug_not_studied_in_source` | B08, B10, B15, B17, B18 | 5 |
| `preclinical_only_not_human_clinical` | B07, B09, B13 | 3 |
| `combination_confounded_and_modest_effect` | B06 | 1 |
| `pharmacokinetic_confounding_and_non_binary_direction` | B16 | 1 |
| `isolated_case_report_and_combination` | B12 | 1 |
| `source_not_accessible` | B02 | 1 |

## Verified drug-absence findings (checked directly against the source text)
These drove five exclusions and were confirmed programmatically, not assumed:
- **B08 Carboplatin / FIRM-ACT:** the regimen is etoposide + doxorubicin + **cisplatin** vs streptozocin. Carboplatin does not appear in the record.
- **B10 Abemaciclib / S05:** full-text search returned **zero** occurrences of "abemaciclib" or "LY2835219"; only palbociclib and ribociclib were tested.
- **B15 Sirolimus / S08:** the drugs tested are pasireotide, everolimus and mitotane; "sirolimus"/"rapamycin" do not appear.
- **B17 Osimertinib and B18 Afatinib / S11:** the source evaluates **gefitinib** monotherapy (2007); neither named drug is an intervention in it.

## Records requiring human adjudication

**1. B02 — Cabozantinib (highest priority; access, not science).**
The source could not be opened at all, so the protocol forced `uncertain`/`NA`. The title indicates a single-arm phase 2 monotherapy trial, which *could* satisfy the strict criteria, but the protocol forbids labelling from a title. **Action:** supply the full text or official abstract and re-rate. This row should not be treated as a substantive negative.

**2. B16 — Sunitinib (genuine judgement call).**
A direct monotherapy phase II — the only excluded record that is otherwise structurally eligible. I excluded it because the authors themselves report (i) a major pharmacokinetic confounder (mitotane-induced CYP3A4 induction reducing sunitinib exposure; r = −0.650 between drug and mitotane levels) and (ii) an interpretation of "modest activity… compares favourably with other targeted treatments", which is not an unambiguous negative. A reasonable second rater could label this `negative` on the basis of median PFS 2.8 months and 24/35 progressive disease. **Action:** adjudicate whether protocol criterion 4 (no major pharmacokinetic/interpretive confounder) or the raw endpoint should dominate.

**3. B06 — Gemcitabine (two exclusion criteria collide).**
132/145 patients also received capecitabine (regimen confounding) *and* the effect is modest (PR 4.9%, median PFS 12 weeks). Both `combination_regimen_not_decomposable` and `mixed_or_limited` apply; the template allows only one reason. **Action:** confirm which primary reason the project wants recorded for double-excluded rows.

**4. B17 / B18 — Osimertinib and Afatinib (medium confidence for access reasons only).**
The exclusion itself is secure (the source studies gefitinib; a 2007 trial cannot have used osimertinib or afatinib), but the abstract body was unavailable so no outcome data could be extracted. **Action:** confirm acceptability of a title-plus-manifest basis for a *non-substantive* exclusion, or supply the official abstract.

**5. B01 — Mitotane (only positive; abstract-only access).**
The recommendation wording is quoted verbatim in the official abstract ("we suggest adjuvant mitotane treatment…", "we recommend either mitotane monotherapy or mitotane, etoposide, doxorubicin and cisplatin"), so I rated `positive` with high confidence. **Action:** optional verification against the full guideline (GRADE strength and evidence quality per recommendation) before locking.

## Observations for a later source-gap search (not used in the ratings)
Several records failed strict eligibility only because the *assigned* source is not a drug-specific clinical study. If the project later wants drug-specific clinical evidence for these agents, more suitable sources would need to be identified for: cabozantinib (full text of S02), carboplatin, abemaciclib, sirolimus, osimertinib and afatinib in ACC. Whether such evidence exists was **not** assessed here, as the protocol restricts the primary rating to the closed source set.

## Confidence distribution
high 16 · medium 2 (B17, B18) · low 1 (B02).
