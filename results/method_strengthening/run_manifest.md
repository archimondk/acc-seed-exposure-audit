# Method-strengthening run manifest

- Analysis version: `method-strengthening-v2-reviewer-null-resolution`
- Command: `python -m analysis.method_strengthening --project-root .`
- Python: `3.14.3`
- NumPy: `2.4.4`
- SciPy: `1.17.1`
- Platform: `Windows-11-10.0.26200-SP0`
- RNG seed: `20260727`
- Null draws: `10000`

## Input SHA-256

- `data/bindex_network/bindex_edges_1304.csv`: `0ac37e507e763b5d20e78571df372ac14350e95b4d19c0360d1324ed35018d63`
- `data/bindex_network/rACC_399_fullSTRING.csv`: `9c9a18d8937b83f1da79cbe8f8e6d6870c30b267ed410fbcc283edc0507c4ab8`
- `data/bindex_network/Sactivity_124_v1.csv`: `5d371ff5d8b8261f1b3a091131506be67fb0a39f39563c333a0bcbd9f5ade475`
- `data/bindex_network/NCI60_potency_124.csv`: `60e5e60e5211990b659b807181d968cb3c06e9c0889635d12d5fea15fb2bb367`
- `data/bindex_network/S_external_curated.csv`: `75e8f4bb3d6ba5d6415647164ac26bd81aa2e479c352650a66ac885eea67aea8`
- `data/ACC_P0.5C_gene_weights_v1.csv`: `a513b96cf50ee10041068ab621ba978acfdf0c89116df32f86a9fcbc043a3897`
- `9606.protein.info.v12.0.txt.gz`: `144de4b0d98c6a7dfde6ddc2591cf88657f27b989eadff4f501450c3ed1f0f1c`
- `9606.protein.links.v12.0.txt.gz`: `3e22f32572211aa341d5b4bd08d30c32e693e294603202120936872f87719d4f`

## Output SHA-256

- `baseline_comparison_primary108.csv`: `1199d2e864a7da01c87682bff636036499dd396a73dae289dd7fa87a5b5da307`
- `centrality_gene399.csv`: `bf7ce717b27fd8d6844657d510995b0a5b21bfadf7137b9794cd7738e2e99cee`
- `centrality_drug108.csv`: `dd57959d829893c810d24c4be8e88ac5a04c0de9208eb629a61268236ed89d12`
- `degree_matched_seed_sets.csv`: `f4b8815e730f84f5f2dbc6f9f8cbc11f8258ab9ac40f8e08e8d13d9c491eb733`
- `random_seed_null_primary108.csv`: `7bc202b8bee56d181e396805049916d663305df1a1dbefb2c284321ee1e47722`
- `method_strengthening_metrics.json`: `227e0ca6d2937f0a286946527674725c079d4a2559bf2106d5a9f4ece53afaa1`
- `claim_evidence_table.md`: `f5dd5bed407dfd6beb5c1f7c298c711a285f95dd9198920c710ff60219dd8ded`
- `method_strengthening_report.md`: `85a1cacfcaae72e29e4c4a55bfef216f3929a5face13b9e878a875c03b10e83f`

## Guardrails

- The C1 primary score is read and compared, not refitted.
- No AUC, PR-AUC or calibration is computed without negative clinical comparators.
- Random-seed empirical P values use add-one correction and BH adjustment across 108 drugs.
