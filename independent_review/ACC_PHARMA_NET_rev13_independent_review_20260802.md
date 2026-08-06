# Independent Peer Review — ACC-PHARMA-NET (rev13)

- **Manuscript**: Seed exposure and normalization sensitivity in network-propagation drug ranking by gene-set averaging: an adrenocortical carcinoma audit
- **Target journal**: Pharmaceutics (MDPI)
- **Manuscript source reviewed**: `ACC_PHARMA_NET_Pharmaceutics_submission_rev13.md` (2026-08-01)
- **Reviewer stance**: independent reviewer; no access to prior review rounds or the authors' internal review records was used
- **Review date**: 2026-08-02
- **Recommendation**: **Major Revision**

---

## 1. Summary of the manuscript

The authors audit a network-propagation drug-ranking construction applied to adrenocortical carcinoma (ACC). A previously published 124-drug/399-gene pharmacogenomic association network (Pharmaceutics 2025, 17, 1421) is contextualized with a 45-gene ACC disease-seed restart vector; gene-set-averaged drug scores (C_ACC, ADRS_comp) are computed in a locked 108-drug complete-case universe. The paper asks six bounded questions (reproducibility, baseline agreement, topology/normalization dependence, positive-control exposure symmetry, seed-specific influence, and the inferential limits of the available literature) and reports: (i) exact recovery of the source network; (ii) strong dependence of the ranking on direct seed overlap (weighted-overlap baseline rho = 0.828, Top-20 Jaccard = 1.000); (iii) material normalization sensitivity; (iv) failure of an outcome-blind breast-cancer positive control, attributed to asymmetric RB1 seed exposure; (v) a frozen RB1 add/remove intervention that passes its full rule set under only 1/4 normalizations (verdict PARTIAL_OR_NOT_SUPPORTED); and (vi) a full 45-seed leave-one-out influence scan (RB1 largest single-drug shift, 57 ranks).

The manuscript is unusually transparent: every analysis is labeled original / exploratory / frozen-before-result / result-known post-hoc; negative results are retained; AI use is disclosed; and an isolated reproduction workflow with SHA-256-frozen inputs and byte-for-byte output comparison is described.

## 2. Verification performed by the reviewer

The following claims were independently re-derived or cross-checked against the frozen data files in the project directory (all consistent with the manuscript):

- CDK4/6 exact one-sided rank-sum P: manuscript 0.0764; independent enumeration of all C(108,3) = 204,156 three-drug subsets gives P = 15589/204156 = 0.07636. ✓
- Section 3.4 hypergeometric P: manuscript 0.17 / 0.07; independent calculation gives 0.1673 / 0.0703. ✓
- Section 3.3 normalization variants: rho 0.622 / 0.641 / 0.640 and Top-20 Jaccard 0.379 / 0.429 / 0.429 match `normalization_sensitivity_summary.csv`. ✓
- Section 3.8 variant-level CDK4/6 P (0.2839 / 0.0104 / 0.0260 / 0.0214) and four-variant family q (0.2839 / 0.0347 / 0.0347 / 0.0347) match `normalization_sensitivity_summary.csv`. ✓
- Section 3.10 breast-arm CDK4/6 group mean 0.5888 vs null 0.5678, empirical P = 0.4538, four-variant q = 0.8486, and ranks (palbociclib 14, abemaciclib 53, ribociclib 68; only palbociclib in Top 25%) match `leakage_audit/arms/B1/*`. ✓
- Section 3.11 intervention numbers (delta-z_ACC = 3.907, delta-z_Breast = 1.961, NC1 ribociclib ≤ 0.015, NC2 maxima 0.319/0.367, F1 not triggered, L1–L4 verdict 1/4 variants) match `leakage_audit/verdict.json` and `arm_summary.csv`. ✓
- Section 3.11.2 RB1 worst-case 57-rank shift and the TERT/BRCA1/MSH6/CHEK2/MEN1 ordering match `results/leave_one_seed_out/leave_one_seed_out_audit.md`. ✓
- Section 3.5 seed-excluded scoring (rho 0.2955 / 0.6459; Top-20 3/20 and 7/20) matches `results/seed_excluded_scoring/seed_excluded_scoring_audit.md`. ✓
- Section 2.11 resolution claims: minimum attainable BH q = 0.0108 over 108 tests, and B ≥ 2160 required for any q < 0.05, are arithmetically correct. ✓
- Abstract numbers agree with the main text. ✓
- Reference list (39 entries) spot-checked; entry 7 corresponds to the source package `pharmaceutics-17-01421-s001.zip` present in the project. ✓

## 3. Major comments

### M1. Table 2 reports a Spearman rho under the wrong label (internal inconsistency)
Table 2's last column is headed "Spearman rho, L_d(count) versus C_ACC percentile". The ACC row reports 0.828, but the frozen observational baseline (`leakage_audit/observational/O1_replacement_baseline.csv`) gives L_count vs C_ACC percentile = **0.822**; the value 0.828 is the **weighted** overlap baseline (L_weighted, Jaccard 1.000). The main text (Section 3.11) correctly attributes rho = 0.828 to L_d(weighted). The breast row (0.577) appears to be the count-based value. The table therefore mixes the two definitions, or its header is wrong.

*Required action*: state which definition is tabulated (count or weighted), use the same definition for both rows, and correct the value/header accordingly. The main text, Table 2, and Figure 5A legend must use one consistent convention.

### M2. Section 3.9: method for the 95% confidence intervals is not described, and the reported intervals are inconsistent with Fisher z-transform values
The manuscript reports: CU-ACC1 rho = 0.84 (95% CI 0.65–0.92, n = 21); CU-ACC2 rho = 0.58 (0.08–0.90, n = 20); NCI-H295R rho = 0.61 (0.14–0.88, n = 16). Independent Fisher z calculation gives (0.64, 0.93), (0.18, 0.81) and (0.16, 0.85) respectively. The discrepancy is large for the n = 20 comparison. Methods (Section 2.11) describe a 2000-draw bootstrap only for the gene-level r_ACC–degree/PageRank correlations, not for the drug-level cross-platform correlations of Section 3.9.

*Required action*: state explicitly how these CIs were computed (bootstrap percentile? Fisher z?), and verify the reported numbers against that method.

### M3. The central positive finding is exploratory; the frozen confirmation is only partially supported — strengthen the framing in the Discussion
The claim that seed exposure dominates the ranking rests on O1–O4, which the authors themselves label exploratory. The only frozen confirmatory design (A2/B2 intervention) passed L1–L3 in just 1/4 normalizations, and the post-hoc per-variant stability analysis (3.11.1) shows NC2 fails in 4/8 cells. The Abstract and Results are carefully worded, but the Discussion (Section 4, first two paragraphs: "seed membership and the fixed drug–gene map determine most decision-relevant ordering in this construction"; "Removing RB1 erased the primary abemaciclib signal") reads more strongly than the frozen evidence supports. Given that the intervention effect is normalization-dependent, the Discussion should state explicitly that the dominance claim is an exploratory, primary-normalization-dependent observation, not a normalization-invariant property.

### M4. Novelty and fit with the journal's scope need a stronger explicit case
This is an audit paper whose principal message is cautionary: for gene-set-averaged propagation rankings built on a fixed pharmacogenomic map, direct seed overlap can dominate the ordering, so the ranking should not be read as an efficacy predictor. References 8–11 already establish degree-aware scoring, normalization sensitivity, benchmarking, and validation-beyond-ranking as known issues. The stated increment ("trace these known risks into a downstream drug ranking in which disease seeds can also be genes in a fixed pharmacogenomic association map", Section 1) is concrete, but the Discussion should (i) explicitly separate which findings are specific to this construction versus re-confirmations of known risks, and (ii) clarify the relationship to the source B-index publication (ref 7): the audit targets the gene-set-averaging scoring layer, not necessarily the B-index similarity network itself (the B_ACC vs B-index comparison is only shown at lambda = 0). Without this explicit separation, readers may infer a broader invalidation of ref 7 than the data show.

### M5. Submission package is incomplete (supplementary assets and metadata)
- Supplementary Figure S6 (referenced as Figure S6A–C in Section 3.9) has no final asset; the legend file itself notes "final asset renumbering deferred until text freeze", and no Fig6 asset exists under `figures/revision/`. Figures S1–S2 are listed as "legacy exploratory analysis" with no source files.
- The Data Availability Statement says the repository URL and DOI "will be inserted before submission". For a manuscript whose selling point is reproducibility, reviewers should be able to access the code repository during review (anonymized link).
- HTML comments in the manuscript (AUTHOR ACTION REQUIRED, ALL-AUTHOR CONFIRMATION REQUIRED) and the missing "All authors have read and agreed..." statement indicate the draft is not submission-ready; this is acceptable for review but must be completed before resubmission.

## 4. Minor comments

1. **Figure asset numbering**: filenames still carry former numbering ("new main Figure 3 uses the former Figure 5 asset; new main Figure 5 uses the former Figure 7 asset"). Rename and re-export all assets to match the final numbering before submission.
2. **Project documentation drift**: README.md still identifies rev5 as the current manuscript; update to rev13 so the repository is self-consistent.
3. **AI tool version**: the Claude version used for the blinded classification is unrecorded (disclosed honestly). Please attempt to recover it from saved metadata; if unrecoverable, state the recovery attempt explicitly in the Supplementary chronology (Table S28).
4. **Table 1 is a genuine contribution**: the "minimum audit requirements" list is transferable and valuable. Consider citing Table 1 explicitly in the Discussion as the paper's actionable output and stating its scope boundaries (which scoring families it applies to).
5. **Abstract (optional)**: the breast-arm counterfactual (13/94 unexposed drugs in the top quartile) is persuasive; adding the breast-arm overlap rho (0.577) would let readers compare the two regimes directly in the Abstract.
6. **Section 3.11.1**: reporting the range of the eight one-sided empirical P values (all < 0.05) would help readers gauge how marginal the post-hoc signals are.

## 5. Strengths worth acknowledging

- Reproducibility practice is exemplary: frozen inputs with SHA-256, isolated runner, byte-for-byte comparison of 31 scientific outputs, 49 passing tests.
- The status-label system (original/exploratory/frozen-before-result/result-known post-hoc) is a model of analytical honesty.
- Negative results (positive-control failure, non-significant CDK4/6 trend in the primary normalization, normalization-dependent intervention) are retained and interpreted rather than discarded.
- The seed-exclusion mitigation (Amendment 5) and the full leave-one-seed-out scan directly address the reviewer-expected questions.
- All numbers I could independently verify match the frozen outputs; the internal consistency of the numerical claims is high.

## 6. Recommendation

**Major revision.** The scientific core (an honest, reproducible audit with a conservative verdict) is sound and the numerical work is verifiable, but the following must be addressed before acceptance: M1 (Table 2 label/value inconsistency), M2 (CI method and values), M3 (Discussion framing matched to the frozen evidence), M4 (novelty and relationship to ref 7 made explicit), and M5 (complete the submission package). With these resolved, the paper would be a useful methodological cautionary contribution; without them, the internal inconsistencies and over-strong Discussion framing would undermine the credibility the paper otherwise earns.
