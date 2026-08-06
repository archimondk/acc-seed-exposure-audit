# Dirichlet component-weight sensitivity run manifest

## Frozen design

- Analysis version: `dirichlet-component-weight-sensitivity-v1`
- Protocol: `dirichlet_component_weight_sensitivity_v1`
- Post-hoc status: `true`
- Active components: `G, R, P, L, S`
- Fixed disease-only seeds: `45`
- Draws: `1000`
- Dirichlet alpha: `[1.0, 1.0, 1.0, 1.0, 1.0]`
- RNG seed: `20260729`
- Primary drug universe: `108`
- Frozen leakage verdict revised: `false`

## Key descriptive results

- ADRS-rank Spearman versus locked ranking:
  median `0.998904`, 5th–95th percentile
  `0.997294`–`0.999628`.
- Top-20 Jaccard versus locked ranking:
  median `1.000000`, 5th–95th percentile
  `0.904762`–`1.000000`.
- Drugs with Top-20 probability >=0.80:
  `20`.

## Quality control

- Baseline r_ACC maximum absolute difference:
  `4.99548072595e-07`.
- Baseline r_ACC Spearman:
  `0.999999905543`.
- Complete rank permutations:
  `true`.
- Maximum RWR iterations:
  `30`.
- Maximum final RWR L1 delta:
  `6.18104087149e-11`.

## Runtime environment

- Python: `3.14.3`
- NumPy: `2.4.4`
- SciPy: `1.17.1`
- Platform: `Windows-11-10.0.26200-SP0`
- Logical CPUs: `12`
- Wall-clock seconds: `42.028`

## Input SHA-256

- `experiments/dirichlet_component_weight_sensitivity_protocol_v1.md`: `206e00aa4dd2fbd13bbaaf75f345b5dd7f5de8d1f2eddaf12e3493e04de8802b`
- `experiments/DIRICHLET_WEIGHT_SENSITIVITY_FREEZE.txt`: `87c016a006eeb9252290612bcd1b699d063bb0b4928983c41efdea42feebd4e6`
- `data/ACC_P0.5C_gene_weights_v1.csv`: `a513b96cf50ee10041068ab621ba978acfdf0c89116df32f86a9fcbc043a3897`
- `data/bindex_network/bindex_edges_1304.csv`: `0ac37e507e763b5d20e78571df372ac14350e95b4d19c0360d1324ed35018d63`
- `data/bindex_network/rACC_399_fullSTRING.csv`: `9c9a18d8937b83f1da79cbe8f8e6d6870c30b267ed410fbcc283edc0507c4ab8`
- `data/bindex_network/Sactivity_124_v1.csv`: `5d371ff5d8b8261f1b3a091131506be67fb0a39f39563c333a0bcbd9f5ade475`
- `data/bindex_network/NCI60_potency_124.csv`: `60e5e60e5211990b659b807181d968cb3c06e9c0889635d12d5fea15fb2bb367`
- `9606.protein.info.v12.0.txt.gz`: `144de4b0d98c6a7dfde6ddc2591cf88657f27b989eadff4f501450c3ed1f0f1c`
- `9606.protein.links.v12.0.txt.gz`: `3e22f32572211aa341d5b4bd08d30c32e693e294603202120936872f87719d4f`

## Output SHA-256

- `results/dirichlet_weight_sensitivity/component_weight_draws.csv`: `a5cbfe6a8cea0211a057ed849145af4d81768616a07a2d95d0878d339dc43a9c`
- `results/dirichlet_weight_sensitivity/drug_rank_draws.csv`: `1ee218868815a3364e0e82182f1a57e05be09441abed312f65185b68a48fd6d8`
- `results/dirichlet_weight_sensitivity/drug_rank_summary.csv`: `e4b73d158ef7852bf5ea6f94f6462df8fbb852695e6e1b91609687b14947d4a6`
- `results/dirichlet_weight_sensitivity/draw_summary.csv`: `31af55f1172eaefb3f004b4d1db0b6e119197f5cb421de99a01fe1349c2cbe54`
- `results/dirichlet_weight_sensitivity/dirichlet_weight_sensitivity_summary.json`: `90af1920baa08781fb0e0b6c1d372dfd6059e11bf6f53b76abf74f5e8ef14f5d`
- `figures/revision/FigS3_dirichlet_weight_sensitivity.png`: `ebac0fd94707c0926cd192f654d57720746b33889e57c9f7c766b92f0e852d0d`
- `figures/revision/FigS3_dirichlet_weight_sensitivity.pdf`: `a56b45c5c4e1a9de169e64f8bbbbdc302951b9619eef4e6f0e61ea063f834ef6`
- `figures/revision/FigS3_dirichlet_weight_sensitivity.svg`: `6bbf02be019befb0591229ebe242366d64d8d390115800152bb8380d4b888761`
