# Positive-control amendment log

## 2026-07-28: provenance documentation only

After the formal 10,000-draw analysis completed, the already downloaded input
URLs, API request body, byte counts and SHA-256 values were copied into
`data/positive_control/erpos_her2neg_tcga/source_manifest.json`, and a human-
readable data card was added. No cohort rule, seed row, seed weight, analysis
parameter, random seed, null draw or result was changed.

The pre-analysis seed freeze remains the controlling record:

- Seed file:
  `results/positive_control/positive_control_seed_frozen.csv`
- Pre-analysis SHA-256:
  `a0de1628c847012e0070d52f437cc0e13227dafc22fc8b62e43ded4e397f9a70`
- Formal-run manifest:
  `results/positive_control/run_manifest.json`
