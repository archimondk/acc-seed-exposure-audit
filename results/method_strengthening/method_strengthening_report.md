# Method strengthening report

## Design lock

- Version: `method-strengthening-v2-reviewer-null-resolution`.
- Primary universe: 108 drugs.
- Full STRING: 19,701 symbols and 929,472 aggregated undirected edges at combined score ≥400.
- Disease-only ACC seeds: 45.
- Degree-matched null: 10,000 draws; RNG seed 20260727; empirical P floor 0.000100.
- Best-case BH q floor across 108 drugs: 0.010799; resolution adequate for q<0.05: True.
- Clinical AUC/PR-AUC/calibration: not estimable (2 strict positives, 0 strict negatives).

## 1. Simple baselines

| Ranking | Spearman vs ADRS | Top-20 Jaccard | CDK4/6 mean rank | Exact P | BH q across rankings |
|---|---:|---:|---:|---:|---:|
| ADRS_comp | 1.000 | 1.000 | 28.33 | 0.0764 | 0.3054 |
| raw_MIPE_potency | 0.689 | 0.333 | 36.33 | 0.1643 | 0.3093 |
| residual_alone | 0.707 | 0.429 | 32.67 | 0.1185 | 0.3093 |
| C_ACC_alone | 0.700 | 0.290 | 39.67 | 0.2145 | 0.3093 |
| association_count | -0.123 | 0.026 | 44.33 | 0.2940 | 0.3093 |
| direct_seed_overlap_fraction | 0.537 | 0.290 | 42.33 | 0.2491 | 0.3093 |
| S_external | 0.026 | 0.143 | 8.67 | 0.0010 | 0.0080 |
| degree_matched_random_seed_mean | 0.252 | 0.143 | 45.00 | 0.3093 | 0.3093 |

Former H1 decision: **retired**. Correlation and Top-20 overlap are retained only as a descriptive composition audit because ADRS_comp is an equal-weight sum of two ranked components.

## 2. STRING centrality

- `r_ACC` vs degree: rho=0.797, 95% bootstrap CI [0.753, 0.833].
- `r_ACC` vs strength: rho=0.798.
- `r_ACC` vs PageRank: rho=0.789, 95% bootstrap CI [0.744, 0.826].
- Partial Spearman with PageRank controlling log-degree: rho=-0.017.

Bootstrap intervals are descriptive because network nodes are not independent observational units.

## 3. Degree-matched random-seed null

- Drugs with BH q<0.05: 2.
- Significant drugs: Pralatrexate (Z=2.73, q=0.0108), Pemetrexed (Z=3.27, q=0.0270).
- CDK4/6 observed mean C_ACC percentile: 0.639.
- CDK4/6 null mean ± SD: 0.560 ± 0.125; empirical P=0.2839.
- H2 decision: **partially_supported**.
- H3 decision: **not_supported**.

## 4. Quality control

- Recomputed vs frozen r_ACC: maximum absolute difference 5e-07; Spearman 0.999999906.
- Effective degree bins: 10 (requested 10).
- Null RWR iterations: 30–38; maximum final L1 delta 9.9e-11.

## 5. Interpretation boundary

These analyses can show non-redundant ranking structure and quantify network/seed dependence. They cannot establish drug efficacy, causal targets or prospective generalization. `S_external` remains a non-independent exploratory baseline because the same literature informed evidence-aware reprioritization.
