# External data sources and redistribution policy

Last compliance check: 2026-08-05. This file distinguishes analysis inputs
that may be included in the public repository from third-party raw data that
must be downloaded by users. SHA-256 values identify the exact local snapshots
used by the project; they do not grant redistribution rights.

## Release rule

- Include project-authored code, protocols, small derived analysis tables,
  frozen input hashes and result tables.
- Include third-party raw data only where the provider explicitly permits
  redistribution and attribution is preserved.
- For large files or data with account-, citation- or use-specific terms,
  publish the official download route, exact filename and SHA-256 only.
- Do not upload publisher PDFs. Open-access PMC/Europe PMC text supplied in the
  Claude audit package must retain its source identifier and provenance.

## Core network and scoring inputs

| Source | Snapshot used | Official access route | Local input / SHA-256 | License or terms status | Public-release decision |
|---|---|---|---|---|---|
| STRING human protein network | STRING v12.0, taxon 9606; combined score threshold ≥400 | <https://string-db.org/cgi/download> and <https://string-db.org/cgi/access>; exact files: <https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz> and <https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz> | `9606.protein.info.v12.0.txt.gz` — `144de4b0d98c6a7dfde6ddc2591cf88657f27b989eadff4f501450c3ed1f0f1c`; `9606.protein.links.v12.0.txt.gz` — `3e22f32572211aa341d5b4bd08d30c32e693e294603202120936872f87719d4f` | STRING states that download data are CC BY 4.0; attribution and disclosure of modifications are required | Omit the 83 MB raw links file from GitHub and let `python -m scripts.fetch_external_inputs` download and verify both official files; include hashes and the small project-derived network tables |
| NCI/NCATS MIPE ACC screen | Project-frozen ACC activity table used by the original workflow | NCI/NCATS source described in the manuscript and repository metadata | `data/ACC_MIPE5_ranked_by_mean_ZAUC.csv` — `8f2c8ae98be1b854b5d85dbef94329cddf5949ee871d9a6e1034d4f2d5dd8559`; `data/bindex_network/Sactivity_124_v1.csv` — `5d371ff5d8b8261f1b3a091131506be67fb0a39f39563c333a0bcbd9f5ade475` | Provider-specific reuse conditions were not captured in the original download record | Include only the minimal frozen project-derived tables required for audit/reproduction, with source citation; do not claim that the project relicenses the underlying screen |
| NCI-60 / CellMiner DTP activity | DTP mean Z-score and raw activity downloads | <https://discover.nci.nih.gov/cellminer/loadDownload.do> and direct processed-data route documented in `data/bindex_network/GET_NCI60_zscore.html` | `DTP_NCI60_ZSCORE.zip` — `17e6a62597a32caa5d43d2e8f81422c108b47bcefd51a3b2269834001da1d7aa`; `DTP_NCI60_RAW.zip` — `f9df751ddacb2565462182176622db145f506256ac57038598f4b45c7956942b`; derived `data/bindex_network/NCI60_potency_124.csv` — `60e5e60e5211990b659b807181d968cb3c06e9c0889635d12d5fea15fb2bb367` | U.S. government source, but the project did not archive a provider license statement with the download | Omit raw ZIP/XLS files from GitHub; include the 124-drug derived table, citation, filenames and hashes |
| Curated ACC seed and drug–gene association tables | Project-frozen curated inputs | Provenance fields and references are carried in the manuscript/supplement and project tables | `data/ACC_P0.5C_gene_weights_v1.csv` — `a513b96cf50ee10041068ab621ba978acfdf0c89116df32f86a9fcbc043a3897`; `data/bindex_network/bindex_edges_1304.csv` — `0ac37e507e763b5d20e78571df372ac14350e95b4d19c0360d1324ed35018d63`; `data/bindex_network/rACC_399_fullSTRING.csv` — `9c9a18d8937b83f1da79cbe8f8e6d6870c30b267ed410fbcc283edc0507c4ab8` | Project-authored derived/curated data are CC BY 4.0 only to the extent that the authors hold the necessary rights; underlying source databases retain their own terms | Include with the attribution and exclusions in `DATA_LICENSE.md`; final confirmation of upstream B-index redistribution rights remains required |
| Curated external ACC evidence | Nineteen prespecified drug–source records, adjudicated v3 | DOI, PubMed and PMC locators are listed in `independent_review/claude_blinded_evidence_review_v1/sources/source_manifest.csv` | `data/bindex_network/S_external_curated.csv` — `31b01d479ef24d51d2b92de3345ddb6050fdd7bba2756fa5fde47ded6a1f35d5`; primary labels `8e01712a00c5f25a5b20d1fd61e81583e8244859b79cdda36d613ac1a2c918f5`; locked Claude output `02d9e1c0d987b7e092f901c8bb3e69be9a6cb38e9cf88b12d5f6ea8be3da2c09` | Project labels are project-authored; article text remains subject to publisher/PMC terms | Include labels, prompt, locked model output, source manifest and open-access text already obtained through PMC/Europe PMC; exclude publisher PDFs and identify unavailable full texts by locator only |

The authoritative frozen list for the core runner is
`reproducibility/input_manifest.json`. If a value here and that manifest ever
differ, stop the release and resolve the mismatch rather than updating a hash
silently.

## Supplementary validation sources

These sources support supplementary or exploratory analyses. They are not
required by the 14-input core runner.

| Source | Version / route | Exact local snapshot | Redistribution decision |
|---|---|---|---|
| DepMap PRISM Repurposing secondary screen | 20Q2 secondary-screen dose-response parameters; DepMap data portal | `prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv`, 290,170,269 bytes, SHA-256 `2ac69a21f1d681fe7447689262b82ca6e3dc90bfef0bd96eb5479b96f424e43d` | Do not upload the raw 290 MB file. Direct users to DepMap, cite the release, and publish only derived ACC-relevant result tables |
| DepMap CRISPR gene effect and model annotation | DepMap 24Q4 Public, Figshare item 27993248 | `CRISPRGeneEffect.csv`, SHA-256 `3d8f3ec6dbf2db7ff834b79b508622ec0b226f3518003fe96ecf5a4fcf167e3b`; `Model.csv`, SHA-256 `ea4e0b2a3bc806f81df62689a5ae75f1a100135727a3d7b8a4c7ccc8815183f8` | Do not redistribute raw files; provide the official Figshare/DepMap route and derived summary tables |
| GDSC2 fitted dose response | Release dated 24 July 2022; <https://www.cancerrxgene.org/downloads/bulk_download> | `GDSC2_fitted_dose_response_24Jul22.xlsx`, SHA-256 `6757b6d2bda06d810dcb67ebdd40e13a8634e86ea5d568710608773d8a05a186` | Do not redistribute the raw workbook; provide official route, citation, filename, hash and derived summaries |
| ACC CellMinerCDB | <https://discover.nci.nih.gov/acc_cellminercdb/> | `ACCcmdb_NCICCR_act.txt` — `9f31e47d7531f99343fbb0f25276d72010068ea1bc543679d45b7608b7f230b8`; `ACCcmdb_NCATS_act.txt` — `de5ce001bff1e2b0b60267d2666ebfb631f6f79b659ac527251eae8e5d66b8a5`; RNA-seq — `1c04288fceafa25cb644262f2fd31f56b675d440eeea6b4b41f0569c8d6d1d7d`; copy number — `97199cebfe800896c43c03b60fd7c409bba5f6b77da6db6b6e4641d987174bd6` | Include project-derived focused tables and the download instructions; omit raw portal downloads unless the authors confirm provider redistribution terms |
| cBioPortal breast-cancer positive control | TCGA BRCA datahub commit `58341090c8bf0368ebe03f7aa95ec5137a8def25`, snapshot 2026-07-28 | Exact URLs, request body, byte counts and hashes are frozen in `data/positive_control/erpos_her2neg_tcga/source_manifest.json` | Do not bundle the 46 MB mutation response in the GitHub source release; retain the manifest and downloader/reproduction route. Dataset-level terms are distinct from the AGPL license of cBioPortal software |

## Literature-evidence package boundary

The Claude package contains four openly accessible full-text extracts obtained
from PMC/Europe PMC and official PubMed metadata for other records. The source
manifest is the controlling provenance record. The historical prompt is kept
unchanged, but its “independent reviewer” wording is a model role instruction;
no independent human rereview was performed. Non-open full texts must be
obtained legally by each user and are not part of the release.

## Citation and modification notices

The manuscript reference list supplies the scientific citations for STRING,
NCI/CellMiner, DepMap/PRISM, GDSC, ACC CellMinerCDB, cBioPortal and each ACC
evidence source. Any redistributed STRING-derived tables should state that
they were filtered and transformed by the authors and are not an official
STRING product. Repository-level code and project-data licences are MIT and
CC BY 4.0, respectively, within the boundaries recorded in `LICENSE`,
`DATA_LICENSE.md` and `release/LICENSE_DECISION_RECORDED.md`.
