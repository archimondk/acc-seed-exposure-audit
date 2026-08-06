# Seed-weight assignment sensitivity run manifest

## Frozen design

- Analysis version: `seed-weight-assignment-sensitivity-v1`
- Protocol: `seed_weight_assignment_sensitivity_v1`
- Post-hoc status: `true`
- Fixed disease-only seeds: `45`
- W1: uniform restart weight `1/45`
- W2 permutations: `1000`
- RNG seed: `20260729`
- Primary drug universe: `108`
- Frozen leakage verdict revised: `false`

## Key descriptive results

- W1 ADRS-rank Spearman versus locked:
  `0.991550`.
- W1 Top-20 Jaccard versus locked:
  `0.904762`.
- W2 ADRS-rank Spearman versus locked:
  median `0.984166`, 5th–95th percentile
  `0.973818`–`0.992771`.
- W2 Top-20 Jaccard versus locked:
  median `0.818182`, 5th–95th percentile
  `0.666667`–`1.000000`.

## Quality control

- Baseline r_ACC maximum absolute difference:
  `4.9954807259e-07`.
- Baseline r_ACC Spearman:
  `0.999999905543`.
- Every W2 draw preserves the baseline weight multiset:
  `true`.
- Complete rank permutations:
  `true`.
- Maximum RWR iterations:
  `30`.
- Maximum final RWR L1 delta:
  `7.79881877284e-11`.

## Runtime environment

- Python: `3.14.3`
- NumPy: `2.4.4`
- SciPy: `1.17.1`
- Platform: `Windows-11-10.0.26200-SP0`
- Logical CPUs: `12`
- Wall-clock seconds: `41.360`

## Input SHA-256

- `experiments/seed_weight_assignment_sensitivity_protocol_v1.md`: `f5fbdc5b2054dba37bb15f0c036f8797af7ac9706e4c851d2e0c515aeffc915d`
- `experiments/SEED_WEIGHT_ASSIGNMENT_SENSITIVITY_FREEZE.txt`: `0c4bb118071c72c829270ae706f8f2f6f9c0d4552aa2f87b4cd33cec36f8e0b1`
- `data/ACC_P0.5C_gene_weights_v1.csv`: `a513b96cf50ee10041068ab621ba978acfdf0c89116df32f86a9fcbc043a3897`
- `data/bindex_network/bindex_edges_1304.csv`: `0ac37e507e763b5d20e78571df372ac14350e95b4d19c0360d1324ed35018d63`
- `data/bindex_network/rACC_399_fullSTRING.csv`: `9c9a18d8937b83f1da79cbe8f8e6d6870c30b267ed410fbcc283edc0507c4ab8`
- `data/bindex_network/Sactivity_124_v1.csv`: `5d371ff5d8b8261f1b3a091131506be67fb0a39f39563c333a0bcbd9f5ade475`
- `data/bindex_network/NCI60_potency_124.csv`: `60e5e60e5211990b659b807181d968cb3c06e9c0889635d12d5fea15fb2bb367`
- `results/dirichlet_weight_sensitivity/draw_summary.csv`: `31af55f1172eaefb3f004b4d1db0b6e119197f5cb421de99a01fe1349c2cbe54`
- `9606.protein.info.v12.0.txt.gz`: `144de4b0d98c6a7dfde6ddc2591cf88657f27b989eadff4f501450c3ed1f0f1c`
- `9606.protein.links.v12.0.txt.gz`: `3e22f32572211aa341d5b4bd08d30c32e693e294603202120936872f87719d4f`

## Output SHA-256

- `results/seed_weight_assignment_sensitivity/uniform_drug_ranks.csv`: `b9b2ea7d367e0f0b531fc8807ec3a684ec235cfee5917acddd6b3cb00627cbd3`
- `results/seed_weight_assignment_sensitivity/permutation_draw_summary.csv`: `cec6ae36a6cca2761833b3641f691146cda23a1a02f4275ebcb46dc8b44709d3`
- `results/seed_weight_assignment_sensitivity/permutation_drug_rank_draws.csv`: `386fb9c24a911ab33edd2ae3514ad750d948ede80e857725d7af401a611097b3`
- `results/seed_weight_assignment_sensitivity/permutation_drug_rank_summary.csv`: `3cb5ea09d459f848b7e7d025f4cfc76bfa80def782cdda91eb4840f635a360e0`
- `results/seed_weight_assignment_sensitivity/seed_weight_assignment_summary.json`: `358ac467ba3e22952f0cc49d8123c0c7beb6f179b1491651d84539529b5be6da`
- `figures/revision/FigS4_seed_weight_assignment_sensitivity.png`: `1690582ba48c590b21cd8f4bbe15b5ba9eade7df3aa72bd4768955df5ada94b3`
- `figures/revision/FigS4_seed_weight_assignment_sensitivity.pdf`: `933c2af7e6785911c00914f3922f830f46a0086583ad988579ad9f8e3a6377ba`
- `figures/revision/FigS4_seed_weight_assignment_sensitivity.svg`: `a8d91daa575e1929fdca8b9c1af6e12385d4d1df2160fed5c5c0211a15a87dcc`
