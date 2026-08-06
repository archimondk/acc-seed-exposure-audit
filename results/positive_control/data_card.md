# Positive-control data card

## Purpose

This dataset defines an outcome-blind ER-positive/HER2-negative breast-cancer
disease context for the reviewer-requested CDK4/6 positive control. It is not
a treatment-response dataset and was not used to fit drug-score weights.

## Provenance and versioning

- Clinical and RPPA inputs: cBioPortal DataHub `brca_tcga`, commit
  `58341090c8bf0368ebe03f7aa95ec5137a8def25`.
- Mutation records and sequenced-sample IDs: public cBioPortal API snapshots
  retrieved 2026-07-28 and frozen by SHA-256.
- Biological source: TCGA Network, *Comprehensive molecular portraits of
  human breast tumours*, Nature 2012, PMID 23000897, PMCID PMC3465532.
- Every URL, request body, byte count and SHA-256 is recorded in
  `data/positive_control/erpos_her2neg_tcga/source_manifest.json`.

## Cohort

- 1,097 clinical patients screened.
- 585 met the locked ER-positive/HER2-negative definition.
- 367 had another definitive ER/HER2 classification and formed the RPPA
  comparator.
- 145 were excluded because ER or HER2 could not be definitively classified.
- 449 control and 292 comparator patients contributed one primary RPPA sample.
- 522 control patients had mutation data.

## Seed construction

- 21 mutation candidates were specified from the TCGA paper's main-text
  significantly mutated genes; 6 passed the locked recurrence thresholds.
- 225 RPPA features were tested; 20 features representing 19 genes passed the
  locked effect-size and BH thresholds.
- 24 genes were retained after union and deduplication.
- CDK4 and CDK6 were force-excluded to prevent copying the accepted targets
  directly into the restart vector.
- Frozen seed SHA-256:
  `a0de1628c847012e0070d52f437cc0e13227dafc22fc8b62e43ded4e397f9a70`.

## Missingness and bias boundary

The cohort uses retrospective TCGA clinical annotations. HER2 calls are
heterogeneous across IHC and FISH; the locked hierarchy uses definitive FISH
when available and excludes unresolved equivocal or indeterminate cases.
Assay availability differs across patients. RPPA measures a limited antibody
panel, and the mutation candidate set is constrained to genes described in
the 2012 study. These limits make the control conservative but do not permit
post-result seed revision.

## Result boundary

The prespecified CDK4/6 control failed all four recovery criteria. This is a
transportability diagnostic of the implemented network procedure, not
evidence against the clinical efficacy of CDK4/6 inhibitors in
ER-positive/HER2-negative breast cancer.
