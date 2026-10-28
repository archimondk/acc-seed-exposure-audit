# ACC parameter sensitivity: local rerun

- Date: 2026-10-19 (Asia/Singapore)
- Working directory: `D:/03-SCI-ACC药物重定位`
- Python: 3.14.3
- Command: `python -m analysis.acc_parameter_sensitivity --string-dir . --out results/acc_parameter_sensitivity_R1_20261019_rerun`
- Exit code: 0
- Script SHA-256: `adc5a8121dd1d5c8ff47fe8f76471545c0268c2729f0c678ddf7a900534e12ee`
- Baseline maximum absolute difference versus archived 108-drug column-minmax C_ACC percentiles: `4.953260024365136e-13`.
- Acceptance tolerance: `5e-11`; baseline gate passed. This is agreement within the declared tolerance, not bitwise identity to the historical scores.
- Output: 88 rows = 11 settings (baseline plus 10 single-factor perturbations) × 2 drug universes × 4 normalizations.
- Files: `acc_parameter_sensitivity.csv`, `acc_parameter_sensitivity_manifest.json`.
- No new null-model draws or inferential P values were calculated. This is result-known post hoc sensitivity analysis.
- `graph_nodes_edge_rows` records source edge-row counts, not necessarily the number of unique undirected edges.
