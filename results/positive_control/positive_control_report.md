# ER+/HER2-negative breast-cancer positive-control report

- Prespecified decision: **fail**.
- Frozen disease seeds: 24.
- Degree-matched null draws: 10,000.
- Primary CDK4/6 group P: 0.4538.
- Primary CDK4/6 group q across variants: 0.8486.

## Primary CDK4/6 ranks

| Drug | Rank / 108 | Top quartile | C_ACC percentile |
|---|---:|:---:|---:|
| Palbociclib | 14 | yes | 0.879 |
| Abemaciclib | 53 | no | 0.514 |
| Ribociclib | 68 | no | 0.374 |

## Prespecified criteria

- FAIL: `primary_group_p_lt_0_05`
- FAIL: `primary_group_q_lt_0_05`
- FAIL: `at_least_two_cdk46_drugs_top_quartile`
- FAIL: `direction_concordant_at_least_three_variants`

Recovery tests implementation transportability for one established disease-mechanism pair. It does not validate ACC candidates or establish efficacy.
