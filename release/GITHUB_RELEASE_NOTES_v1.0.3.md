# v1.0.3 — author metadata and reproducibility-package correction

This corrective release fixes the corresponding author's given name from
**Yuxing** to **Yuhang Xia** in citation, licence, reproducibility and Zenodo
metadata. It also restores compact author-generated result files that were
unintentionally absent from the v1.0.2 source archive.

Additional release hardening:

- adds root-level `.zenodo.json` metadata with the verified five-author order;
- adds deterministic LF rules so frozen protocol hashes verify across Windows
  and Linux clones;
- restores the final rev13 orchestration manifest and verification report;
- includes the degree-matched seed table and the two 108,000-row sensitivity
  draw tables required by the regression suite;
- removes machine-specific paths from downloader instructions;
- keeps all provider-controlled third-party raw inputs outside the archive.

Verification: **65 tests passed** in a clean release worktree after fetching
the two STRING v12.0 inputs from the official provider and checking their
frozen SHA-256 values.

No scientific result, ranking or conclusion changed. Code is MIT licensed;
author-generated derived data and documentation are CC BY 4.0 within the scope
defined in `DATA_LICENSE.md`.
