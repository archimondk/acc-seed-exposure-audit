# Pharmaceutics DOI release readiness — rev13

Status updated: 2026-08-07

## Current status

The public repository is <https://github.com/archimondk/acc-seed-exposure-audit>.
Zenodo concept DOI `10.5281/zenodo.21824511` identifies all versions; v1.0.2
version DOI `10.5281/zenodo.21825627` remains the previous immutable snapshot.
This v1.0.3 corrective release fixes the corresponding author's given name to
`Yuhang Xia`, restores compact author-generated outputs omitted from v1.0.2,
adds deterministic line-ending rules for frozen hashes and does not change any
scientific result.

## Release boundary

### Included

- analysis and reproduction code, tests and exact dependency declarations;
- author-generated derived tables, compact null outputs and sensitivity draws;
- protocols, input hashes, run manifests and verification reports;
- figure-source tables and author-created figure assets;
- citation metadata, licence files and a root `.zenodo.json` file.

### Linked and hashed, not redistributed

- STRING raw protein-info and links archives;
- NCATS MIPE and CellMiner/NCI-60 raw downloads;
- cBioPortal/TCGA, ACC_CellMinerCDB, DepMap/PRISM and GDSC source exports;
- publisher PDFs and upstream supplementary files without explicit
  redistribution permission.

`DATA_SOURCES.md` records providers and retrieval requirements. The repository
downloader fetches the two STRING inputs from the official provider and checks
their frozen SHA-256 values.

## Verification gates

| Requirement | Status | Evidence |
|---|---|---|
| Public version-controlled repository | pass | GitHub repository and v1.0.2 release resolve publicly |
| Correct creators | pass for v1.0.3 | Five-author order is fixed; `Yuhang Xia` is the corresponding author |
| Code/data licence | pass | MIT for code; CC BY 4.0 for author-generated data/documentation; explicit third-party exclusions |
| Reproducible package | pass | Clean release worktree with official STRING inputs: 65/65 tests passed |
| Sensitive data / credentials | pass | No access token, password, personal local path or identifiable participant data detected |
| New archival DOI | pending external automation | Zenodo will mint the v1.0.3 version DOI after the GitHub release is published |

## Publication sequence

1. Commit the verified v1.0.3 release candidate and publish GitHub release
   `v1.0.3`.
2. Confirm that Zenodo archived the same tag with creators headed by Han Zhang
   and including `Xia, Yuhang`.
3. Record the new version DOI and verify the archive checksum and file list.
4. Insert that version DOI into the local manuscript and regenerate the final
   submission PDF.

The manuscript must cite the version DOI, while the repository citation file
may retain the concept DOI so citations continue to resolve across versions.
