# Pharmaceutics DOI release readiness — rev13

Status updated: 2026-08-05

## Current decision

The project is **not yet DOI-publishable without further author action**. Five creator names, affiliations and contact emails were supplied on 1 August 2026, official-site-supported addresses/postcodes were added, and the corresponding author confirmed the creator order. MIT was selected for author-owned code and CC BY 4.0 for author-generated data/documentation, with explicit third-party exclusions. The repository owner/URL, institutional-name approval, release tag, Zenodo archival record and DOI remain unavailable. ORCID identifiers are desirable but optional. A DOI must not be invented, and a public archival record must not be created before its remaining metadata and redistribution boundaries are approved.

## Release scope

### Include in the DOI-bearing release

- analysis and reproduction code;
- exact dependency lock and execution instructions;
- derived drug–gene network tables that the upstream publication permits to be redistributed;
- derived ACC seed, score, ranking, null-summary and sensitivity tables;
- protocol files and their SHA-256 values;
- blinded prompt/source manifests, locked Claude classifications and adjudication records, excluding any source full text that cannot legally be redistributed;
- figure-source tables and final manuscript/supplement Markdown;
- a machine-readable file manifest and checksums;
- a clear AI-use and analysis-timing record.

### Link, hash and document rather than redistribute

- STRING raw protein-info and links archives;
- NCATS MIPE raw downloads;
- CellMiner/NCI-60 raw downloads;
- cBioPortal/TCGA source exports;
- publisher PDFs or other copyrighted full text;
- any upstream B-index supplement whose license does not authorize relicensing.

For these inputs, the release should provide the official source URL, version or accession, retrieval date, byte count and SHA-256. This preserves reproducibility without claiming redistribution rights.

## Required author decisions

1. All-author approval of the release metadata and official English affiliation forms; the creator order is fixed; supply ORCID identifiers if available.
2. Repository owner or organization.
3. Confirmation that derived B-index tables may be redistributed under the upstream article's licence.
4. Whether the archival release should be public immediately or reserved until manuscript submission.

## Exact publication sequence

1. Create the version-controlled repository and exclude non-redistributable raw inputs.
2. Run the isolated reproduction and Amendment 4 and Amendment 5 implementations from a clean checkout.
3. Generate the final checksum manifest.
4. Tag the exact submitted state, for example `v1.0.0`.
5. Archive that tag in Zenodo or an equivalent DOI-minting repository.
6. Record the concept DOI and version DOI.
7. Replace the manuscript's Data Availability submission-gate sentence with the real version DOI URL.
8. Verify that the DOI landing page exposes the code, derived data, licenses, creators, version and citation.

## Data Availability text after DOI minting

> Analysis code, derived data, protocol files, blinded-review records, figure-source tables and reproduction manifests are available in the versioned archival release at **[insert the real version DOI URL]**. Third-party raw data that cannot be redistributed are documented by provider URL, version, retrieval date and SHA-256 in the release manifest. Primary sources include the B-index publication and its Supplementary Materials, NCATS MIPE 5.0, STRING v12.0, CellMiner NCI-60, ACC_CellMinerCDB, TCGA and cBioPortal.

Do not place this final wording in the submitted manuscript until the DOI resolves publicly.

## Hard-gate status

| Requirement | Status | Evidence/action |
|---|---|---|
| Version-controlled repository | planned/open | Author plans a GitHub upload; no repository URL or `.git` directory is yet present in the workspace |
| Authenticated publishing route | planned/unverified | GitHub is planned, but no real repository or Zenodo/OSF archival link has been supplied |
| Verified creators | partial | Five names, degrees, affiliations, emails and official-site-supported addresses/postcodes are present; creator order is confirmed; institutional English-name and final release-metadata approval remain open; ORCIDs are optional |
| Code/data licence | ready | MIT is recorded in `LICENSE`; CC BY 4.0 scope and third-party exclusions are recorded in `DATA_LICENSE.md` |
| Redistribution review | partial | Third-party raw inputs identified for exclusion; upstream derived-table rights still require confirmation |
| Reproducible scientific outputs | ready | The isolated core verified 14/14 input hashes and 38/38 byte-identical outputs; the 21-stage rev13 orchestration passed, including executable Amendments 4-5, and the final suite reported 87 passed/1 obsolete layout skip |
| DOI | blocked | Cannot mint truthfully until the preceding items are resolved |

This is an external publication gate, not a scientific-analysis failure.
