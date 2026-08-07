# Protocol Amendment 5: seed-excluded scoring audit

- Analysis implementation: `seed-excluded-scoring-v1`.
- Protocol SHA-256: `6c124dc5ee2f3c77634623858a8afa9453bafb3b20d48bea66526b0b587354c2`.
- Status: result-known post-hoc; protocol frozen before outputs.

## Primary results

| Comparison with locked ranking | Spearman rho | Top-20 intersection | Top-20 Jaccard |
|---|---:|---:|---:|
| Seed-excluded C_ACC versus locked C_ACC | 0.2955 | 3/20 | 0.0811 |
| ADRS_nonseed versus locked ADRS_comp | 0.6459 | 7/20 | 0.2121 |

Entered: Acalabrutinib, Afatinib, Bosutinib, Chlorambucil, Fulvestrant, Lapatinib, Melphalan, Osimertinib, Panobinostat, Ribociclib, Romidepsin, Tivozanib, Tucatinib.

Left: Abemaciclib, Actinomycin D, Brigatinib, Crizotinib, Cytarabine, Daunorubicin, Docetaxel, Fluorouracil, Ibrutinib, Irinotecan, Mercaptopurine, Mitomycin, Paclitaxel.

## Coverage and focal drugs

Among 108 drugs, 46 had at least one direct seed association and 62 had none. After exclusion, 2 drugs had no remaining gene, 17 had one and 89 had at least two. The locked reference mean was mu_0 = 0.0525027124.

| Drug | Seed genes removed | Non-seed genes retained | Locked ADRS rank | ADRS_nonseed rank | Locked C_ACC rank | Seed-excluded C_ACC rank |
|---|---:|---:|---:|---:|---:|---:|
| Abemaciclib | 1 | 4 | 8 | 24 | 8 | 40 |
| Palbociclib | 1 | 25 | 26 | 69 | 32 | 79 |
| Ribociclib | 0 | 8 | 51 | 20 | 79 | 43 |

## Seed-excluded Top 20

| New rank | Drug | Locked rank | Seed genes removed | Non-seed genes retained | C_ACC,nonseed | ADRS_nonseed | Rank change |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | Olaparib | 12 | 0 | 1 | 0.049785 | 0.953271 | -11 |
| 2 | Celecoxib | 13 | 0 | 2 | 0.040677 | 0.939252 | -11 |
| 3 | Homoharringtonine | 15 | 0 | 1 | 0.041069 | 0.911215 | -12 |
| 4 | Mitotane | 17 | 0 | 1 | 0.041092 | 0.892523 | -13 |
| 5 | Axitinib | 19 | 0 | 1 | 0.039897 | 0.880841 | -14 |
| 6 | Romidepsin | 21 | 0 | 1 | 0.040875 | 0.878505 | -15 |
| 7 | Bosutinib | 28 | 0 | 1 | 0.040281 | 0.810748 | -21 |
| 8 | Panobinostat | 29 | 0 | 1 | 0.043244 | 0.808411 | -21 |
| 9 | Lapatinib | 33 | 0 | 1 | 0.040281 | 0.787383 | -24 |
| 10 | Chlorambucil | 32 | 0 | 6 | 0.023885 | 0.771028 | -22 |
| 11 | Tivozanib | 36 | 0 | 1 | 0.039419 | 0.766355 | -25 |
| 12 | Osimertinib | 35 | 0 | 2 | 0.033761 | 0.761682 | -23 |
| 13 | Fulvestrant | 39 | 0 | 2 | 0.039470 | 0.757009 | -26 |
| 14 | Ixazomib | 3 | 3 | 2 | 0.036252 | 0.757009 | +11 |
| 15 | Acalabrutinib | 37 | 0 | 1 | 0.044529 | 0.757009 | -22 |
| 16 | Tucatinib | 45 | 0 | 2 | 0.035955 | 0.742991 | -29 |
| 17 | Melphalan | 38 | 0 | 6 | 0.023037 | 0.733645 | -21 |
| 18 | Belinostat | 6 | 1 | 4 | 0.029922 | 0.728972 | +12 |
| 19 | Afatinib | 42 | 0 | 9 | 0.023031 | 0.719626 | -23 |
| 20 | Ribociclib | 51 | 0 | 8 | 0.026686 | 0.710280 | -31 |

## Largest composite-rank movements

| Drug | Locked rank | ADRS_nonseed rank | Absolute shift | Seed genes removed | Non-seed genes retained |
|---|---:|---:|---:|---:|---:|
| Mercaptopurine | 20 | 80 | 60 | 3 | 24 |
| Idarubicin | 24 | 73 | 49 | 1 | 10 |
| Cyclophosphamide | 50 | 96 | 46 | 1 | 15 |
| Dabrafenib | 64 | 108 | 44 | 4 | 36 |
| Allopurinol | 43 | 86 | 43 | 1 | 12 |
| Hydroxyurea | 31 | 74 | 43 | 1 | 6 |
| Ibrutinib | 16 | 59 | 43 | 1 | 12 |
| Palbociclib | 26 | 69 | 43 | 1 | 25 |
| Daunorubicin | 4 | 46 | 42 | 2 | 10 |
| Irinotecan | 1 | 42 | 41 | 1 | 10 |

## Interpretation boundary

This deterministic sensitivity arm shows that direct seed contributions can be removed computationally. It does not establish that the remaining score is more biologically valid or clinically predictive, and it must be reported alongside rather than silently substituted for the locked analysis.

## Input SHA-256

- `data/bindex_network/bindex_edges_1304.csv`: `0ac37e507e763b5d20e78571df372ac14350e95b4d19c0360d1324ed35018d63`
- `data/bindex_network/rACC_399_fullSTRING.csv`: `9c9a18d8937b83f1da79cbe8f8e6d6870c30b267ed410fbcc283edc0507c4ab8`
- `data/ACC_P0.5C_gene_weights_v1.csv`: `a513b96cf50ee10041068ab621ba978acfdf0c89116df32f86a9fcbc043a3897`
- `results/primary_analysis/ADRS_comp_primary_108.csv`: `6e0c45dbb193d1660862b91d540d236e85e35a0359405cb6718a71c8b8318c8f`
