# Leakage-audit claim–evidence table

| Claim ID | Candidate manuscript claim | Evidence | Design/status | Permitted wording |
|---|---|---|---|---|
| LA-C01 | The primary ACC propagation score is strongly structured by direct seed adjacency. | L_weighted vs C_ACC: Spearman ρ = 0.828, P = 2.4×10^-28; Top-20 Jaccard = 1.000. | Exploratory O1; observed before interventional freeze. | “Strongly structured by” or “closely reproduced by”; do not call causal. |
| LA-C02 | Drugs without a direct ACC-seed association could not enter the primary top quartile. | 0/62 entered; highest percentile 0.673. | Exploratory O2; hard partition in this dataset. | Restrict to the locked 108-drug primary construction. |
| LA-C03 | The association between C_ACC and matched-null z is not independent of seed overlap. | Partial Spearman ρ = 0.042, P = 0.667 after controlling ranked L_count. | Exploratory O3. | “No residual monotonic association was detected”; do not claim proof of zero effect. |
| LA-C04 | RB1 exposure drives the primary min–max abemaciclib signal. | z_A1 = 2.734; z_A2 = −1.173; Δz_ACC = 3.907. B1 = 0.021; B2 = 1.982; Δz_Breast = 1.961. | Frozen 2×2 intervention; L1–L3 pass under primary variant. | “Produced a large, symmetric shift under the primary construction.” |
| LA-C05 | The RB1 effect is normalization dependent. | L1–L3 jointly passed in 1/4 variants; Δz_ACC = 0.680–0.731 and Δz_Breast = 0.405–1.144 in the three rank-based variants. | Frozen L4 failed. | “Failed the prespecified cross-normalization robustness criterion.” |
| LA-C06 | The intervention did not globally move drug z scores. | NC2 passed for 106 unexposed drugs; max \|Δz\| = 0.319 ACC and 0.367 breast. | Frozen negative control. | “Unexposed-drug shifts remained below 0.5.” |
| LA-C07 | Ribociclib behaved as an unexposed mechanism control. | \|Δz\| = 0.014 ACC and 0.015 breast. | Frozen NC1 passed. | “Ribociclib was essentially unchanged.” |
| LA-C08 | The breast-arm primary effect was insensitive to moderate RB1-weight perturbation. | Δz_Breast = 1.772, 1.961 and 1.965 at 0.5×, 1.0× and 1.5× median weight. | Prespecified B2 weight sensitivity. | “Direction and approximate magnitude were stable across ±50% weight perturbation.” |
| LA-C09 | Seed–target leakage is the normalization-invariant explanation of the framework. | L4 = 1/4, threshold ≥3/4. | **Not supported.** | Must not be claimed. |
| LA-C10 | The original positive control failed under all four directions. | Primary observed group mean 0.5888 vs null mean 0.5678; other three variants below null. | Rev6 factual statement was incorrect. | “Only the primary variant exceeded its null mean (1/4), below the ≥3/4 rule.” |

