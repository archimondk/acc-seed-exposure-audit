# Independent evidence-review protocol

## 1. Review design

- Unit of assessment: one prespecified drug–source pair.
- Primary review type: closed-set and blinded.
- Evidence condition: use the primary source and its supplement when
  available. Record `abstract_only` when full text is unavailable.
- Do not infer a label from the article title alone.
- Do not add new literature to the primary rating. Potential missing sources
  may be listed separately in `reviewer_notes` for a later source-gap search.
- Use `uncertain` or `NA` when the available source is insufficient.

## 2. Controlled fields

### evidence_domain

- `clinical`: human ACC guideline, trial, cohort, case series, or case report.
- `preclinical`: direct experiments in ACC cells, organoids, xenografts, or
  other non-human ACC models.
- `class_extrapolation`: no direct evidence for the named drug; inference is
  transferred from a related drug, target, or mechanism.
- `uncertain`

### evidence_design

Choose the closest value:

- `guideline_and_clinical_standard`
- `prospective_single_arm_phase_2`
- `randomized_phase_3_regimen`
- `retrospective_multicenter_regimen`
- `retrospective_case_series`
- `multicenter_case_series`
- `single_pediatric_case_report`
- `phase_1_expansion`
- `in_vitro_ACC_cell_lines`
- `no_direct_ACC_study`
- `other`
- `uncertain`

### drug_specificity

- `direct`: the named drug was administered, experimentally tested, or
  explicitly recommended for ACC.
- `class_only`: only a related drug/class/target was studied.
- `none`: the named drug was not studied and no valid class inference is
  established by the source.
- `uncertain`

Mere mention in the introduction, discussion, or reference list does not
count as direct evidence.

### exposure_context

- `monotherapy`
- `monotherapy_or_standard`
- `combination_regimen`
- `none`
- `uncertain`

If multiple agents were administered and the source does not isolate the
named drug's contribution, use `combination_regimen`.

### direction

- `positive`: attributable human clinical evidence or a guideline supports
  use/benefit of the named drug.
- `negative`: attributable human clinical evidence clearly shows lack of
  activity or failure of the relevant clinical endpoint without a major
  exposure or attribution confounder.
- `mixed_or_limited`: clinical results are modest, heterogeneous,
  pharmacokinetically confounded, or otherwise not cleanly binary.
- `regimen_effect`: the evidence applies to a combination regimen and the
  named drug's individual effect cannot be isolated.
- `positive_preclinical`: direct ACC preclinical activity only.
- `contextual_case`: isolated case-level evidence that cannot support a
  general drug-effect label.
- `class_extrapolated`: evidence is transferred from another drug or class.
- `no_direct_evidence`: the source does not provide evidence for the named
  drug.
- `uncertain`

### strict_candidate_eligible

Use `yes` only when all conditions hold:

1. The evidence is human ACC clinical evidence.
2. The source directly evaluates or explicitly recommends the named drug.
3. Exposure is attributable to the drug, or it is an established ACC
   standard explicitly covered by a guideline.
4. The clinical direction can be assigned as positive or negative without a
   major regimen, exposure, pharmacokinetic, population, or interpretive
   confounder.
5. The evidence is stronger than an isolated case report.

Otherwise use `no`. A direct monotherapy study is not automatically eligible
when its result is mixed or materially confounded.

### strict_binary_label

- `positive`: only when `strict_candidate_eligible=yes` and direction is
  unambiguously positive.
- `negative`: only when `strict_candidate_eligible=yes` and direction is
  unambiguously negative.
- `NA`: all other cases.

## 3. Required extraction

For every row record:

- exact treatment and comparator/regimen;
- ACC population or model;
- sample size;
- primary efficacy or activity outcome;
- numerical results used for the judgment;
- source locator: page, section, table, figure, supplement, or paragraph;
- full-text versus abstract-only access;
- one concise exclusion reason when strict eligibility is `no`;
- confidence: `high`, `medium`, or `low`.

## 4. Independence rules

- Do not use model scores or drug ranks.
- Do not attempt to reproduce a target number of positive or negative labels.
- Do not read any first-review audit or manuscript interpretation.
- Complete all 19 rows before seeing any prior rating.
- The first completed file must be preserved unchanged for agreement
  analysis.

