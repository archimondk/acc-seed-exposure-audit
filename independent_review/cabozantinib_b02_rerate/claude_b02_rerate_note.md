# B02 cabozantinib — access-completion rerating note

**Scope:** record B02 only. The locked 19-row second-review file
(`claude_blinded_evidence_review_v1/claude_review_locked.csv`) was **not opened for editing
and not modified** (verified unchanged, 20 lines incl. header). No first-review label, model
score, ranking or agreement result was consulted. Output: `claude_b02_rerated.csv` (one row,
original column order).

## What changed and why
At the first pass S02 was absent from the package (no full text, no abstract), so the protocol
required `uncertain`/`NA` with the exclusion reason `source_not_accessible` — an **access**
limitation, explicitly not a scientific judgement. The official PubMed record (PMID 38608694,
*Lancet Oncol* 2024) now supplies the design, exposure and effect estimates, so the record can be
rated on its merits.

## Rating
`clinical` · `prospective_single_arm_phase_2` · `direct` · `monotherapy` · direction `positive`
· **strict_candidate_eligible = yes** · **strict_binary_label = positive** · confidence `high`.

Check against the five strict conditions:

| Condition | Assessment |
|---|---|
| 1. Human ACC clinical evidence | Met — adults with advanced ACC, histologically confirmed |
| 2. Source directly evaluates the named drug | Met — cabozantinib is the sole prospectively evaluated intervention |
| 3. Exposure attributable to the drug | Met — true monotherapy (60 mg daily); the entry rule requiring serum mitotane <2 mg/L removes the CYP3A4 pharmacokinetic confounder |
| 4. Direction assignable without major confounder | Met — prespecified primary endpoint (PFS at 4 months) met, 13/18 = 72.2% (95% CI 46.5–90.3), median PFS 6 months; authors interpret as promising efficacy |
| 5. Stronger than an isolated case report | Met — registered prospective phase 2 trial (NCT03370718), n = 18, median follow-up 36.8 months |

## Consistency with how the other rows were treated
The differential outcome versus **B16 (sunitinib)** — the other direct single-arm phase 2
monotherapy in the set, which I excluded — is traceable to condition 4 rather than to the drug:

| | B02 cabozantinib | B16 sunitinib |
|---|---|---|
| Pharmacokinetic confounder | Excluded by design (mitotane <2 mg/L required) | Present and author-reported (mitotane-induced CYP3A4; r = −0.650) |
| Primary endpoint result | Met: 72.2% PFS at 4 months | 14.3–15.4% response (PFS ≥12 weeks) |
| Authors' interpretation | "promising efficacy" | "modest activity" |

## Transparency caveats (recorded, but not label-changing under the current protocol)
- Single-arm, no comparator — permitted by the protocol's controlled vocabulary
  (`prospective_single_arm_phase_2`) and not excluded by the strict conditions.
- Small sample (n = 18) and single centre, as expected in a rare cancer.
- Industry funding (Exelixis); the trial is nonetheless investigator-initiated and registered.
- Access is abstract-only, but the structured record contains design, n, intervention, primary
  endpoint, effect estimate with 95% CI, safety and interpretation.
- **Conditional statement:** if the project later requires a randomised comparator for a
  `positive` strict label, this row should become `NA`. That requirement is not in the present
  protocol, so it was not applied.

## Net effect on the second review
Strict-eligible records rise from 1/19 to 2/19 (B01 mitotane, B02 cabozantinib), both `positive`.
The B02 adjudication item raised in `claude_review_adjudication_notes.md` is now resolved; the
remaining adjudication items (B16, B06, B17/B18, B01 full-guideline verification) are unaffected.
