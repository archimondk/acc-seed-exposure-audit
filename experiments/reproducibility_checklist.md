# Network-propagation seed-exposure audit: reproducibility checklist

Status date: 2026-08-05 (rev13).

## Environment and locked inputs

- [x] Python version is pinned by `.python-version`; dependencies are pinned
  by `requirements-lock.txt`.
- [x] The isolated ACC workflow verifies 14 frozen inputs against
  `reproducibility/input_manifest.json`.
- [x] The blinded Claude-classification inputs remain read-only; adjudication is stored
  only in `data/evidence/evidence_labels_v3_adjudicated.csv`.
- [x] Isolated manifests record operating system, processor identifiers,
  logical CPU count and wall-clock time.

## ACC analysis

- [x] Primary universe fixed at 108 MIPE/NCI-60 complete-case drugs.
- [x] Primary shrinkage fixed at `k=3`; sensitivity covers `k={1,3,5,10}`.
- [x] RWR restart probability fixed at `alpha=0.4`.
- [x] Degree-matched null fixed at `B=10000`, RNG seed `20260727`, with the
  add-one empirical-P and BH-resolution gates.
- [x] Four normalization branches were run: column-min–max,
  column+gene-rank, uniform-ratio+gene-rank and symmetric+gene-rank.
- [x] The primary ACC null produced 2/108 drugs at BH `q<0.05`; the primary
  ACC CDK4/6 group result was `P=0.2839`.
- [x] The final isolated rev13 reproduction verified 38/38 scientific outputs,
  14/14 inputs and 9/9 expected core figure files in 1130.6 seconds:
  `repro_outputs/run_20260803T151144Z/results/reproducibility/run_manifest.json`.

## Evidence and reporting boundaries

- [x] Same-resource NCI-CCR/MIPE results are called cross-platform
  concordance, not independent validation.
- [x] ROC-AUC, PR-AUC and calibration are not reported for the two-positive,
  zero-negative strict clinical subset.
- [x] The local GDSC2 and PRISM resources contain no usable ACC/NCI-H295R
  continuous response endpoint.
- [x] Categorical strict evidence labels agreed in 19/19 records between the
  primary human curator and blinded Claude classification; this is explicitly
  described as model–human concordance, not human inter-rater reliability.
- [x] The continuous `S_external` score was not independently rescored.
- [x] On 1 August 2026, the author team decided not to pursue an independent
  human rereview. The manuscript retains the Claude provenance, describes the
  19/19 comparison as model–human traceability, and makes no human-agreement or
  human inter-rater reliability claim.

## Outcome-blind positive control

- [x] Control disease and accepted mechanism frozen as
  ER-positive/HER2-negative breast cancer and CDK4/6 inhibition.
- [x] Clinical and RPPA files pinned to cBioPortal DataHub commit
  `58341090c8bf0368ebe03f7aa95ec5137a8def25`; API snapshots and all inputs
  recorded by URL, byte count and SHA-256 in
  `data/positive_control/erpos_her2neg_tcga/source_manifest.json`.
- [x] Cohort rules, mutation/RPPA thresholds, direct-target exclusions and
  four recovery criteria frozen in `experiments/positive_control_protocol.md`
  before control-rank inspection.
- [x] Frozen outcome-blind seed: 24 genes; CDK4 and CDK6 excluded; SHA-256
  `a0de1628c847012e0070d52f437cc0e13227dafc22fc8b62e43ded4e397f9a70`.
- [x] Four normalization variants and 10,000 degree-matched draws completed
  with RNG seed `20260727`.
- [x] Full artifacts retained: 108 x 10,000 drug-null matrices, 40,000
  CDK4/6 group-null values and 240,000 degree-matching records.
- [x] Independent recomputation from the raw null matrices reproduced every
  group P value.
- [x] Prespecified result: **fail**. Primary CDK4/6 `P=0.4538`,
  four-variant BH `q=0.8486`; only palbociclib was in the Top 25%; direction
  concordance was 1/4, below the required 3/4; all four recovery criteria
  failed.
- [x] The current formal run completed without stderr in 727.1 seconds and is recorded in
  `results/positive_control/run_manifest.json`.
- [x] No post-result seed, threshold or success-criterion revision was made.

## Seed-exposure intervention audit

- [x] Exploratory O1–O4 analyses are explicitly labeled as observed before
  the interventional protocol freeze.
- [x] Interventional protocol frozen at SHA-256
  `36c9638ade80bb761f6e8481889575b1b80feb460bbd51c9d5d68617a4155e85`;
  `FREEZE.txt` records that A2/B2 results were unobserved.
- [x] Frozen arms: A1 (ACC, 45 seeds), A2 (A1−RB1, 44), B1 (breast, 24),
  B2 (B1+RB1, 25), and B2 RB1-weight sensitivities at 0.5×/1.5×.
- [x] Every arm regenerated 10,000 degree-matched sets with RNG seed
  `20260727`; row counts are 450,000/440,000/240,000/250,000/250,000/250,000.
- [x] Every arm contains 432 drug–variant rows (108 drugs × 4 variants);
  stderr is empty and matched-seed hashes differ across arms.
- [x] A2 includes RB1 in the eligible replacement pool; A1/A2 matched-seed
  hashes differ.
- [x] Frozen verdict: `PARTIAL_OR_NOT_SUPPORTED`. L1 = 3.907, L2 = 1.961
  and L3 = 1.947 passed under `column_minmax`, but L4 failed because only
  1/4 variants passed all three rules (required ≥3/4).
- [x] F1 did not trigger (`z_A2(abemaciclib) = −1.173`). Under the primary
  column-min–max construction, NC1 and NC2 passed; among 106 unexposed drugs,
  maximum absolute shifts were 0.319 (ACC) and 0.367 (breast).
- [x] The seed-exposure figure (rev13 main Figure 5; source asset formerly
  Figure 7) is generated from frozen CSV/JSON outputs in PDF, SVG and
  1000-dpi PNG at 170 mm width; visual and integrity checks passed.
- [x] Amendment 1 was registered on 2026-07-29 after all six arm results and
  the frozen verdict were observed. The frozen protocol file and its SHA-256
  remain unchanged; the amendment is stored in a separate sidecar.
- [x] Post-hoc within-variant analysis reproduced all eight disease × variant
  cells. Abemaciclib was directionally concordant in 8/8 cells and had
  one-sided empirical `P<0.05` in 8/8; ribociclib ranged from −0.150 to
  +0.015. NC2 passed in 4/8 cells and failed in 4/8.
- [x] In the four NC2-pass cells abemaciclib ranked first by absolute
  delta-z; in the four NC2-fail cells it ranked second to fifth and maximum
  unexposed movement reached 2.572 z.
- [x] Pairwise association-set resolution was computed for all 5,778 locked
  drug pairs. Seventeen drugs have `n_d=1`, two pairs are exact duplicates,
  and three pairs have Jaccard similarity `>=0.8`.
- [x] Rev8 synchronized the amendment boundary, scale-free results, network
  resolution limitation, supplementary index and terminology without changing
  `PARTIAL_OR_NOT_SUPPORTED`.

## Dirichlet disease-component-weight sensitivity

- [x] Protocol Amendment 2 was frozen before any Dirichlet result was
  generated at SHA-256
  `206e00aa4dd2fbd13bbaaf75f345b5dd7f5de8d1f2eddaf12e3493e04de8802b`.
- [x] The analysis retained the 45 disease-biology-only seeds and varied only
  the five active components G/R/P/L/S; T and the five
  therapy-response/vulnerability genes remained excluded.
- [x] Exactly 1,000 `Dirichlet(1,1,1,1,1)` vectors were generated with NumPy
  `default_rng` and RNG seed `20260729`; every component and restart vector
  was positive and summed to one.
- [x] Every draw recomputed the primary RWR, min–max `r_ACC`, shrinkage
  `C_ACC` percentile and `ADRS_comp`, while retaining the locked activity
  residual and 0.50/0.50 downstream ADRS weights.
- [x] Baseline reconstruction passed the frozen gate: maximum `r_ACC`
  difference `4.99548e-7`, Spearman `0.999999906`; all 1,000 draws converged
  and yielded complete rank permutations 1–108.
- [x] Median ADRS-rank Spearman correlation with the locked ranking was
  `0.9989` (5th–95th percentile `0.9973–0.9996`); median Top-20 Jaccard was
  `1.000` (5th–95th percentile `0.9048–1.000`).
- [x] Abemaciclib remained Top 20 in 100.0% and Top 10 in 98.8% of draws;
  palbociclib and ribociclib had median ranks 27 and 49.
- [x] All 1,000 component vectors, 108,000 draw-by-drug rows, per-drug
  summaries, run manifest and Figure S3 are retained.
- [x] The analysis is labelled post hoc and descriptive; no empirical-P/FDR
  claim was made, degree-matched nulls were not regenerated, and
  `PARTIAL_OR_NOT_SUPPORTED` was not revised.

## Seed-weight magnitude and gene-assignment sensitivity

- [x] Protocol Amendment 3 was frozen before W1/W2 output generation at
  SHA-256
  `f5fbdc5b2054dba37bb15f0c036f8797af7ac9706e4c851d2e0c515aeffc915d`.
- [x] The analysis used the five-component weights that actually enter the
  disease-only primary restart vector, not the CSV `ACC_weight` field that
  includes the excluded therapeutic component T.
- [x] The 45 baseline raw weights ranged from `0.260` to `0.770`; normalized
  restart mass ranged from `0.0114085` to `0.0337867`.
- [x] W1 assigned exactly `1/45` to every retained seed. W2 generated exactly
  1,000 independent permutations with NumPy `default_rng` and RNG seed
  `20260729`; every draw preserved the baseline weight multiset exactly.
- [x] Every arm recomputed the primary RWR, min–max `r_ACC`, shrinkage
  `C_ACC` percentile and `ADRS_comp`, while retaining the locked activity
  residual and 0.50/0.50 downstream ADRS weights.
- [x] Baseline reconstruction passed the frozen gate: maximum `r_ACC`
  difference `4.99548e-7`, Spearman `0.999999906`; all 1,001 new arm columns
  converged and yielded complete rank permutations 1–108.
- [x] W1 rank Spearman was `0.9915` and Top-20 Jaccard was `0.9048`; 81 exact
  ranks changed, with mean absolute shift `2.37` and maximum shift `19`.
- [x] W2 median rank Spearman was `0.9842` (5th–95th percentile
  `0.9738–0.9928`), entirely below the Amendment 2 interval
  `0.9973–0.9996`. Median Top-20 Jaccard was `0.8182`
  (5th–95th percentile `0.6667–1.000`).
- [x] Abemaciclib remained Top 20 in 99.9% of W2 draws but Top 10 in 47.9%;
  its median W2 rank was 11 (5th–95th percentile 7–16).
- [x] All W1 ranks, 1,000 W2 summaries, 108,000 W2 draw-by-drug rows,
  per-drug summaries, run manifest and Figure S4 are retained.
- [x] The result is labelled post hoc and descriptive. No equivalence margin,
  empirical-P/FDR claim or regenerated degree-matched null was introduced,
  and `PARTIAL_OR_NOT_SUPPORTED` was not revised.
- [x] The matching Amendment 2 minimum rho and W1 rho were checked against
  the original draw-by-drug files. The vectors differ for 85/108 drugs but
  both have a sum of squared rank deviations of 1,774, giving
  `rho=0.9915497252`; this is not a transcription copy.
- [x] A result-known post-hoc W2–MIPE diagnostic used all 108 drugs in every
  draw and the locked `MIPE_potency_pct` orientation. The locked
  `C_ACC`–MIPE rho was `0.0607098455`; the W2 median was `0.0629629453`
  (5th–95th percentile `0.0325075143–0.0913400942`).
- [x] The locked value was at the 45.8th W2 percentile and 542/1,000
  permutations had rho at least as large. The manuscript treats this as
  absence of a detectable advantage on the available criterion, not proof of
  equivalence or absence of all biological information.
- [x] The approximate Fisher z-transformed 95% CI for the locked rho at
  `n=108` was independently recomputed as `−0.1297526959 to +0.2468517626`.
  The manuscript therefore retains a power limitation and does not claim
  demonstrated equivalence or exclusion of a small association.
- [x] `results/primary_analysis/ADRS_context_only_16.csv` confirms that
  trilaciclib is one of the 16 exclusions because MIPE activity was
  unavailable (context-only C_ACC rank 32/124); the RB1 edge is present in
  `data/bindex_network/bindex_edges_1304.csv`.
- [x] W2–MIPE inputs and calculation details are recorded in
  `results/seed_weight_assignment_sensitivity/W2_MIPE_external_criterion_audit.md`.

## Full leave-one-seed-out influence scan

- [x] Protocol Amendment 4 was frozen before any full-scan output at SHA-256
  `6de390b48185acc90703b249d444586ca0680850d9180126102afeb847d753d4`.
- [x] All 45 disease-only seeds were removed one at a time, remaining weights
  were renormalized, and all four declared normalization variants were rerun
  in the locked 108-drug universe (180 seed-by-variant runs).
- [x] The graph, restart probability, shrinkage, activity residual and drug
  universe remained fixed; no degree-matched null, P value or outcome-based
  branch selection was introduced.
- [x] Column, symmetric and uniform-reference propagation runs converged below
  the frozen numerical tolerance.
- [x] Minimum global ADRS-rank Spearman correlation was `0.9167`; minimum
  Top-20 Jaccard was `0.7391`.
- [x] RB1 ranked 1/45 for both maximum single-drug movement and maximum
  directly exposed-drug movement (57 ranks). TERT, BRCA1, MSH6, CHEK2 and
  MEN1 also produced large local shifts.
- [x] The complete 45-seed-by-four-variant results, input hashes, convergence
  values and interpretation boundary are retained in
  `results/leave_one_seed_out/leave_one_seed_out_audit.md`.
- [x] The scan is labelled reviewer-requested, result-known post-hoc and
  descriptive in Methods, Results, Discussion and Supplementary Table S28.
- [x] `analysis/leave_one_seed_out.py` and its regression tests reproduce the
  full audit in the isolated core and complete rev13 orchestration.

## Seed-excluded scoring sensitivity

- [x] Protocol Amendment 5 was frozen before any seed-excluded output at
  SHA-256
  `6c124dc5ee2f3c77634623858a8afa9453bafb3b20d48bea66526b0b587354c2`.
- [x] For each of the 108 locked drugs, directly associated disease seeds
  were removed while `k=3`, the association-edge-weighted locked
  `mu_0=0.052502712423312885`, the activity residual and the drug universe
  remained fixed.
- [x] Drugs with no remaining non-seed association received the prespecified
  neutral prior `mu_0`; this occurred for doxorubicin and pralatrexate.
- [x] Seed-excluded `C_ACC` had rank Spearman correlation `0.2955` with
  locked `C_ACC` and retained 3/20 locked context Top-20 drugs
  (Jaccard `0.0811`).
- [x] Seed-excluded `ADRS` had rank Spearman correlation `0.6459` with
  locked `ADRS_comp` and retained 7/20 locked Top-20 drugs
  (Jaccard `0.2121`).
- [x] Abemaciclib moved from composite rank 8 to 24, palbociclib from 26 to
  69 and ribociclib from 51 to 20; the latter lost no gene and moved because
  percentiles were recomputed within the fixed universe.
- [x] The complete formula, input hashes, coverage, Top-20 turnover and
  focal-drug shifts are retained in
  `results/seed_excluded_scoring/seed_excluded_scoring_audit.md`.
- [x] The primary Amendment 5 calculation reads the archived 12-decimal
  `residual_pct` values named in the frozen input list. A full-precision
  upstream reconstruction changed rho from `0.6459001256` to `0.6449169739`
  and exchanged only the exact-tied ribociclib/decitabine boundary ranks;
  Top-20 overlap and the interpretation were unchanged.
- [x] The analysis is labelled result-known post-hoc and is interpreted as a
  direct mitigation/sensitivity analysis, not proof that the replacement
  ordering is biologically superior.
- [x] `analysis/seed_excluded_scoring.py` and its regression tests reproduce
  Amendment 5 in the isolated core and complete rev13 orchestration.

## Remaining author-controlled submission items

- [x] Verified names, degrees, affiliations, emails and corresponding-author
  contact details supplied on 1 August 2026 were inserted into the rev13
  working manuscript. The user-requested one-time random order for authors
  3–5 was Zhenhan Mo, Yujing Zhang and Haixia Yang.
- [x] AI-use disclosure now identifies OpenAI Codex (GPT-5) and Anthropic
  Claude, states that Claude performed model-based rather than human review,
  and records that the exact Claude model/version could not be recovered.
- [x] Added official-site-supported street addresses, postcodes and province
  fields to both affiliations; institutional English naming still requires
  corresponding-author confirmation.
- [x] Added an order-based CRediT draft to the manuscript. No role is assigned
  for an independent human rereview because that task was not performed.
- [x] Corresponding-author confirmation on 1 August 2026 fixed the final
  author order as Han Zhang, Yuhang Xia, Zhenhan Mo, Yujing Zhang and Haixia
  Yang.
- [x] Funding status was confirmed and the manuscript now states: “This
  research received no external funding.”
- [ ] Obtain all-author approval of every CRediT role and the final manuscript;
  confirm the institutional English forms and any additional acknowledgments
  using `AUTHOR_METADATA_REQUIRED_rev13.md`.
- [ ] Publish the authorized code/data package to GitHub/Zenodo, then insert
  only the real URL and DOI; see
  `release/PHARMACEUTICS_DOI_RELEASE_READINESS_rev13.md`.
- [ ] GitHub upload is author-planned; record the real repository owner, URL,
  license, release/tag and commit after they exist.
- [x] Remade Figure 1 from `analysis/figure1_rev13.py` with a visible four-item
  status key, per-block text badges, explicit 124/108/16 universe boundaries,
  direct-overlap and seed-excluded arms, and the bounded output statement.
- [x] Exported editable SVG, vector PDF and 6692 × 6456 px PNG at 1000 dpi;
  integrity lint and original-size visual QA passed.
- [x] Regenerated final Figure S1 and Figure S2 assets in PDF/SVG/PNG with
  plotted source CSVs and deterministic regression tests.
- [x] Exported the former Figure 6 content as Supplementary Figure S6 in
  PDF/SVG/PNG and verified its manifest and rev13 panel mapping.

## Full rev13 orchestration

- [x] `python -m scripts.reproduce_rev13 --project-root .` completed all 21
  stages with zero non-zero return codes.
- [x] The schema-2 orchestration manifest safely resumed at B2 after a Windows
  host interruption and preserved only prior successful stages.
- [x] Summed stage time was 6046.3 seconds (100.8 minutes); B2, B2_lo and
  B2_hi all completed after resume.
- [x] Final scientific regression suite: 87 passed, 1 obsolete rev5 Word-layout
  skip and 0 failures.
- [x] Human-readable verification is archived at
  `results/reproducibility/REV13_VERIFICATION_REPORT_20260805.md`.
