# Figure source data — ACC-PHARMA-NET

Self-contained CSVs holding exactly the data behind each manuscript figure, so the figures can be re-plotted independently. The plotting script is `data/bindex_network/make_figures.py`.

| File | Figure / panel | Columns |
|---|---|---|
| *(none — schematic)* | **Figure 1** | Framework diagram; no underlying data |
| `Fig2_lambda0_reproduction.csv` | **Figure 2a,b** | drugA, drugB, published_Bindex, recovered_BACC_lambda0, abs_diff (1337 drug pairs) |
| `Fig3_component_values.csv` | **Figure 3** | per-drug: C_ACC_pct, ACCrelative_resid_pct, S_Bneighbor_pct, NCI60_potency, MIPE_potency, S_external (Spearman matrix computed from these) |
| `Fig4a_benchmark_labels_scores.csv` | **Figure 4** | drug, comp_score (0.5·C_ACC+0.5·residual), label (positive/negative) |
| `Fig4b_permutation_null_AUC.csv` | **Figure 4** | permutation_AUC — all 1001 label-permutation AUCs (observed AUC = 0.40) |
| `Fig5a_weightgrid_comp_ranks.csv` | **Figure 5a** | drug × wC∈{0..1 step .05}: rank of each drug in the 2-component comp score (108 complete-case drugs) |
| `Fig5b_CDK46_enrichment_null.csv` | **Figure 5b** | mean_rank_random3 — 50,000 random 3-drug mean ranks (null) |
| `Fig5b_CDK46_observed.txt` | **Figure 5b** | observed CDK4/6 mean rank + P/FDR |
| `Fig6abc_crossplatform_concordance.csv` | **Figure 6a–c** | cell_line, drug, NCICCR_neglogIC50, NCATS_neglogIC50 |
| `Fig6d_biomarker_expression.csv` | **Figure 6d** | category, gene, then 4 cell lines + 6 surgical tumors (log2 FPKM+1) |

## Re-plotting notes
- **Fig 3** heatmap = Spearman correlation among the columns of `Fig3_component_values.csv` (pairwise-complete).
- **Fig 4** = histogram of `Fig4b...` with a vertical line at the observed AUC (0.40) from `Fig4a...`.
- **Fig 5a** = boxplot per drug across the wC columns of `Fig5a...`.
- **Fig 5b** = histogram of `Fig5b...null.csv` with a vertical line at the observed value in `...observed.txt`.
- **Fig 6a–c** = scatter of the two IC50 columns, split by `cell_line`.
- All values are the exact numbers used in the published figures; empty cells = not measured for that drug.
