# Full 45-seed leave-one-out influence audit

- Implementation: `leave-one-seed-out-v1`.
- Protocol SHA-256: `6de390b48185acc90703b249d444586ca0680850d9180126102afeb847d753d4`.
- Classification: reviewer-requested, result-known post-hoc descriptive analysis.
- Seeds: 45; locked drugs: 108; variants: 4; deterministic runs: 180.

## Headline result

Across variants, the minimum ADRS-rank Spearman correlation was 0.9167 and the minimum Top-20 Jaccard overlap was 0.7391. RB1 produced a worst-case shift of 57 ranks and ranked 1/45 by maximum single-drug movement and 1/45 by maximum directly exposed-drug movement.

28 of 45 seeds had no direct association to a drug in the locked universe. Largest worst-case shifts: RB1 (57), TERT (54), BRCA1 (51), MSH6 (51), CHEK2 (50), MEN1 (49), EGFR (43), BRCA2 (40).

## Variant-level summary

| Variant | Minimum rho | Median rho | Minimum Top-20 Jaccard | Largest absolute rank shift | Seeds causing any >=10-rank shift |
|---|---:|---:|---:|---:|---:|
| column_minmax | 0.9167 | 0.9998 | 0.7391 | 57 | 13 |
| column_gene_rank | 0.9945 | 0.9988 | 0.8182 | 19 | 12 |
| uniform_ratio_gene_rank | 0.9877 | 0.9982 | 0.8182 | 25 | 18 |
| symmetric_gene_rank | 0.9923 | 0.9986 | 0.8182 | 20 | 11 |

## Seed-level worst case across four variants

| Seed | Normalized weight | Directly exposed drugs | Minimum rho | Minimum Top-20 Jaccard | Maximum absolute shift | Maximum number shifting >=10 | Maximum exposed shift | Maximum unexposed shift |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| AKT1 | 0.01316 | 0 | 0.9986 | 0.9048 | 6 | 0 | 0 | 6 |
| APC | 0.02896 | 0 | 0.9981 | 0.9048 | 9 | 0 | 0 | 9 |
| ATM | 0.01711 | 1 | 0.9990 | 0.9048 | 6 | 0 | 1 | 6 |
| ATR | 0.01711 | 3 | 0.9929 | 0.9048 | 37 | 1 | 37 | 13 |
| ATRX | 0.03050 | 0 | 0.9944 | 0.9048 | 12 | 2 | 0 | 12 |
| BRCA1 | 0.01755 | 7 | 0.9601 | 0.8182 | 51 | 5 | 51 | 7 |
| BRCA2 | 0.01755 | 4 | 0.9912 | 0.9048 | 40 | 1 | 40 | 9 |
| CCNE1 | 0.02896 | 0 | 0.9979 | 1.0000 | 13 | 1 | 0 | 13 |
| CDK4 | 0.03028 | 2 | 0.9982 | 0.9048 | 10 | 1 | 7 | 10 |
| CDKN2A | 0.03203 | 0 | 0.9970 | 0.9048 | 10 | 1 | 0 | 10 |
| CHEK2 | 0.01711 | 4 | 0.9736 | 0.9048 | 50 | 3 | 50 | 6 |
| CTNNB1 | 0.03379 | 0 | 0.9944 | 0.9048 | 17 | 1 | 0 | 17 |
| CYP11A1 | 0.01404 | 0 | 0.9996 | 1.0000 | 6 | 0 | 0 | 6 |
| CYP17A1 | 0.01404 | 0 | 0.9984 | 1.0000 | 15 | 1 | 0 | 15 |
| CYP21A2 | 0.01404 | 0 | 0.9995 | 0.9048 | 5 | 0 | 0 | 5 |
| DAXX | 0.02742 | 1 | 0.9949 | 0.9048 | 14 | 2 | 0 | 14 |
| EGFR | 0.01272 | 5 | 0.9683 | 0.8182 | 43 | 5 | 43 | 7 |
| FRS2 | 0.01887 | 0 | 0.9961 | 0.9048 | 12 | 2 | 0 | 12 |
| HSD3B2 | 0.01404 | 0 | 0.9995 | 1.0000 | 4 | 0 | 0 | 4 |
| IGF1R | 0.01711 | 0 | 0.9976 | 0.9048 | 8 | 0 | 0 | 8 |
| IGF2 | 0.03247 | 0 | 0.9963 | 0.8182 | 13 | 1 | 0 | 13 |
| IL7R | 0.01931 | 1 | 0.9958 | 0.9048 | 18 | 2 | 14 | 18 |
| KDR | 0.01141 | 1 | 0.9982 | 1.0000 | 10 | 1 | 10 | 8 |
| KRAS | 0.01975 | 0 | 0.9966 | 0.8182 | 14 | 1 | 0 | 14 |
| LRP1B | 0.01887 | 0 | 0.9946 | 0.8182 | 13 | 2 | 0 | 13 |
| MED12 | 0.02304 | 1 | 0.9932 | 0.9048 | 36 | 3 | 36 | 15 |
| MEN1 | 0.03050 | 7 | 0.9778 | 0.8182 | 49 | 5 | 49 | 25 |
| MLH1 | 0.02391 | 0 | 0.9992 | 0.9048 | 5 | 0 | 0 | 5 |
| MSH2 | 0.02391 | 4 | 0.9990 | 0.9048 | 6 | 0 | 6 | 5 |
| MSH6 | 0.02391 | 3 | 0.9862 | 0.9048 | 51 | 1 | 51 | 6 |
| MTOR | 0.01316 | 0 | 0.9985 | 0.9048 | 10 | 1 | 0 | 10 |
| NF1 | 0.02501 | 0 | 0.9977 | 0.9048 | 9 | 0 | 0 | 9 |
| NR5A1 | 0.02260 | 0 | 0.9978 | 1.0000 | 10 | 1 | 0 | 10 |
| PIK3CA | 0.01536 | 0 | 0.9972 | 0.9048 | 16 | 1 | 0 | 16 |
| PMS2 | 0.02391 | 0 | 0.9992 | 0.9048 | 5 | 0 | 0 | 5 |
| PRKAR1A | 0.03093 | 0 | 0.9902 | 0.8182 | 20 | 8 | 0 | 20 |
| PTCH1 | 0.01887 | 3 | 0.9873 | 0.9048 | 30 | 3 | 30 | 9 |
| RB1 | 0.03203 | 2 | 0.9731 | 0.9048 | 57 | 5 | 57 | 19 |
| RPL22 | 0.02523 | 0 | 0.9984 | 0.9048 | 15 | 1 | 0 | 15 |
| STAR | 0.01404 | 0 | 0.9996 | 0.9048 | 6 | 0 | 0 | 6 |
| TERF2 | 0.02764 | 0 | 0.9974 | 0.8182 | 13 | 1 | 0 | 13 |
| TERT | 0.03072 | 16 | 0.9167 | 0.7391 | 54 | 16 | 54 | 19 |
| TP53 | 0.03379 | 0 | 0.9973 | 0.9048 | 8 | 0 | 0 | 8 |
| VEGFA | 0.01141 | 0 | 1.0000 | 1.0000 | 0 | 0 | 0 | 0 |
| ZNRF3 | 0.03181 | 0 | 0.9990 | 1.0000 | 9 | 0 | 0 | 9 |

## Complete seed-by-variant results

| Seed | Variant | rho | Top-20 Jaccard | Median abs. shift | Max abs. shift | n >=5 | n >=10 | exp n | exp median | exp max | unexp median | unexp max |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| AKT1 | column_gene_rank | 0.9994 | 1.0000 | 0.0 | 4 | 0 | 0 | 0 | NA | NA | 0.0 | 4 |
| AKT1 | column_minmax | 1.0000 | 1.0000 | 0.0 | 0 | 0 | 0 | 0 | NA | NA | 0.0 | 0 |
| AKT1 | symmetric_gene_rank | 0.9986 | 0.9048 | 1.0 | 6 | 4 | 0 | 0 | NA | NA | 1.0 | 6 |
| AKT1 | uniform_ratio_gene_rank | 0.9991 | 1.0000 | 0.5 | 4 | 0 | 0 | 0 | NA | NA | 0.5 | 4 |
| APC | column_gene_rank | 0.9982 | 0.9048 | 1.0 | 9 | 5 | 0 | 0 | NA | NA | 1.0 | 9 |
| APC | column_minmax | 0.9998 | 1.0000 | 0.0 | 3 | 0 | 0 | 0 | NA | NA | 0.0 | 3 |
| APC | symmetric_gene_rank | 0.9986 | 0.9048 | 1.0 | 6 | 4 | 0 | 0 | NA | NA | 1.0 | 6 |
| APC | uniform_ratio_gene_rank | 0.9981 | 1.0000 | 1.0 | 7 | 4 | 0 | 0 | NA | NA | 1.0 | 7 |
| ATM | column_gene_rank | 0.9990 | 0.9048 | 0.0 | 6 | 2 | 0 | 1 | 0.0 | 0 | 0.0 | 6 |
| ATM | column_minmax | 0.9999 | 1.0000 | 0.0 | 2 | 0 | 0 | 1 | 1.0 | 1 | 0.0 | 2 |
| ATM | symmetric_gene_rank | 0.9992 | 0.9048 | 0.0 | 6 | 2 | 0 | 1 | 0.0 | 0 | 0.0 | 6 |
| ATM | uniform_ratio_gene_rank | 0.9992 | 1.0000 | 1.0 | 3 | 0 | 0 | 1 | 0.0 | 0 | 1.0 | 3 |
| ATR | column_gene_rank | 0.9983 | 0.9048 | 1.0 | 13 | 2 | 1 | 3 | 1.0 | 1 | 1.0 | 13 |
| ATR | column_minmax | 0.9929 | 0.9048 | 1.0 | 37 | 1 | 1 | 3 | 2.0 | 37 | 1.0 | 4 |
| ATR | symmetric_gene_rank | 0.9984 | 0.9048 | 1.0 | 9 | 3 | 0 | 3 | 1.0 | 2 | 1.0 | 9 |
| ATR | uniform_ratio_gene_rank | 0.9991 | 1.0000 | 1.0 | 5 | 1 | 0 | 3 | 1.0 | 1 | 1.0 | 5 |
| ATRX | column_gene_rank | 0.9966 | 0.9048 | 1.0 | 12 | 9 | 1 | 0 | NA | NA | 1.0 | 12 |
| ATRX | column_minmax | 0.9997 | 1.0000 | 0.0 | 4 | 0 | 0 | 0 | NA | NA | 0.0 | 4 |
| ATRX | symmetric_gene_rank | 0.9944 | 0.9048 | 1.0 | 12 | 17 | 2 | 0 | NA | NA | 1.0 | 12 |
| ATRX | uniform_ratio_gene_rank | 0.9948 | 0.9048 | 1.0 | 12 | 13 | 2 | 0 | NA | NA | 1.0 | 12 |
| BRCA1 | column_gene_rank | 0.9989 | 0.9048 | 0.0 | 7 | 3 | 0 | 7 | 0.0 | 7 | 0.0 | 6 |
| BRCA1 | column_minmax | 0.9601 | 0.8182 | 2.0 | 51 | 16 | 5 | 7 | 38.0 | 51 | 2.0 | 7 |
| BRCA1 | symmetric_gene_rank | 0.9975 | 1.0000 | 1.0 | 7 | 8 | 0 | 7 | 0.0 | 2 | 1.0 | 7 |
| BRCA1 | uniform_ratio_gene_rank | 0.9987 | 1.0000 | 1.0 | 6 | 2 | 0 | 7 | 3.0 | 6 | 1.0 | 6 |
| BRCA2 | column_gene_rank | 0.9988 | 0.9048 | 0.5 | 7 | 2 | 0 | 4 | 0.5 | 1 | 0.5 | 7 |
| BRCA2 | column_minmax | 0.9912 | 1.0000 | 1.0 | 40 | 3 | 1 | 4 | 6.5 | 40 | 1.0 | 4 |
| BRCA2 | symmetric_gene_rank | 0.9988 | 1.0000 | 1.0 | 9 | 2 | 0 | 4 | 0.5 | 1 | 1.0 | 9 |
| BRCA2 | uniform_ratio_gene_rank | 0.9988 | 0.9048 | 1.0 | 4 | 0 | 0 | 4 | 1.0 | 3 | 1.0 | 4 |
| CCNE1 | column_gene_rank | 0.9986 | 1.0000 | 0.0 | 9 | 3 | 0 | 0 | NA | NA | 0.0 | 9 |
| CCNE1 | column_minmax | 0.9999 | 1.0000 | 0.0 | 3 | 0 | 0 | 0 | NA | NA | 0.0 | 3 |
| CCNE1 | symmetric_gene_rank | 0.9987 | 1.0000 | 1.0 | 6 | 3 | 0 | 0 | NA | NA | 1.0 | 6 |
| CCNE1 | uniform_ratio_gene_rank | 0.9979 | 1.0000 | 1.0 | 13 | 5 | 1 | 0 | NA | NA | 1.0 | 13 |
| CDK4 | column_gene_rank | 0.9990 | 0.9048 | 0.0 | 5 | 1 | 0 | 2 | 0.0 | 0 | 0.5 | 5 |
| CDK4 | column_minmax | 0.9995 | 1.0000 | 0.0 | 7 | 1 | 0 | 2 | 3.5 | 7 | 0.0 | 2 |
| CDK4 | symmetric_gene_rank | 0.9982 | 1.0000 | 1.0 | 7 | 4 | 0 | 2 | 0.5 | 1 | 1.0 | 7 |
| CDK4 | uniform_ratio_gene_rank | 0.9987 | 1.0000 | 1.0 | 10 | 2 | 1 | 2 | 0.5 | 1 | 1.0 | 10 |
| CDKN2A | column_gene_rank | 0.9979 | 1.0000 | 1.0 | 8 | 7 | 0 | 0 | NA | NA | 1.0 | 8 |
| CDKN2A | column_minmax | 0.9998 | 1.0000 | 0.0 | 3 | 0 | 0 | 0 | NA | NA | 0.0 | 3 |
| CDKN2A | symmetric_gene_rank | 0.9970 | 0.9048 | 1.0 | 10 | 8 | 1 | 0 | NA | NA | 1.0 | 10 |
| CDKN2A | uniform_ratio_gene_rank | 0.9974 | 1.0000 | 1.0 | 9 | 7 | 0 | 0 | NA | NA | 1.0 | 9 |
| CHEK2 | column_gene_rank | 0.9994 | 1.0000 | 0.0 | 6 | 2 | 0 | 4 | 0.5 | 1 | 0.0 | 6 |
| CHEK2 | column_minmax | 0.9736 | 0.9048 | 1.0 | 50 | 6 | 3 | 4 | 35.5 | 50 | 1.0 | 6 |
| CHEK2 | symmetric_gene_rank | 0.9993 | 1.0000 | 1.0 | 5 | 1 | 0 | 4 | 1.5 | 3 | 0.5 | 5 |
| CHEK2 | uniform_ratio_gene_rank | 0.9997 | 1.0000 | 0.0 | 2 | 0 | 0 | 4 | 1.0 | 1 | 0.0 | 2 |
| CTNNB1 | column_gene_rank | 0.9982 | 0.9048 | 1.0 | 6 | 5 | 0 | 0 | NA | NA | 1.0 | 6 |
| CTNNB1 | column_minmax | 0.9998 | 1.0000 | 0.0 | 3 | 0 | 0 | 0 | NA | NA | 0.0 | 3 |
| CTNNB1 | symmetric_gene_rank | 0.9944 | 0.9048 | 1.5 | 17 | 14 | 1 | 0 | NA | NA | 1.5 | 17 |
| CTNNB1 | uniform_ratio_gene_rank | 0.9972 | 1.0000 | 1.0 | 9 | 7 | 0 | 0 | NA | NA | 1.0 | 9 |
| CYP11A1 | column_gene_rank | 0.9997 | 1.0000 | 0.0 | 4 | 0 | 0 | 0 | NA | NA | 0.0 | 4 |
| CYP11A1 | column_minmax | 1.0000 | 1.0000 | 0.0 | 0 | 0 | 0 | 0 | NA | NA | 0.0 | 0 |
| CYP11A1 | symmetric_gene_rank | 0.9996 | 1.0000 | 0.0 | 6 | 1 | 0 | 0 | NA | NA | 0.0 | 6 |
| CYP11A1 | uniform_ratio_gene_rank | 0.9996 | 1.0000 | 0.0 | 5 | 1 | 0 | 0 | NA | NA | 0.0 | 5 |
| CYP17A1 | column_gene_rank | 0.9996 | 1.0000 | 0.0 | 5 | 2 | 0 | 0 | NA | NA | 0.0 | 5 |
| CYP17A1 | column_minmax | 1.0000 | 1.0000 | 0.0 | 0 | 0 | 0 | 0 | NA | NA | 0.0 | 0 |
| CYP17A1 | symmetric_gene_rank | 0.9996 | 1.0000 | 0.0 | 4 | 0 | 0 | 0 | NA | NA | 0.0 | 4 |
| CYP17A1 | uniform_ratio_gene_rank | 0.9984 | 1.0000 | 0.0 | 15 | 1 | 1 | 0 | NA | NA | 0.0 | 15 |
| CYP21A2 | column_gene_rank | 0.9995 | 1.0000 | 0.0 | 4 | 0 | 0 | 0 | NA | NA | 0.0 | 4 |
| CYP21A2 | column_minmax | 1.0000 | 1.0000 | 0.0 | 0 | 0 | 0 | 0 | NA | NA | 0.0 | 0 |
| CYP21A2 | symmetric_gene_rank | 0.9996 | 0.9048 | 0.0 | 5 | 1 | 0 | 0 | NA | NA | 0.0 | 5 |
| CYP21A2 | uniform_ratio_gene_rank | 0.9995 | 1.0000 | 0.0 | 5 | 1 | 0 | 0 | NA | NA | 0.0 | 5 |
| DAXX | column_gene_rank | 0.9966 | 0.9048 | 1.0 | 14 | 5 | 2 | 1 | 0.0 | 0 | 1.0 | 14 |
| DAXX | column_minmax | 0.9996 | 1.0000 | 0.0 | 4 | 0 | 0 | 1 | 0.0 | 0 | 0.0 | 4 |
| DAXX | symmetric_gene_rank | 0.9967 | 1.0000 | 1.0 | 9 | 10 | 0 | 1 | 0.0 | 0 | 1.0 | 9 |
| DAXX | uniform_ratio_gene_rank | 0.9949 | 1.0000 | 1.0 | 14 | 13 | 2 | 1 | 0.0 | 0 | 1.0 | 14 |
| EGFR | column_gene_rank | 0.9989 | 1.0000 | 1.0 | 5 | 2 | 0 | 5 | 1.0 | 4 | 1.0 | 5 |
| EGFR | column_minmax | 0.9683 | 0.8182 | 1.0 | 43 | 16 | 5 | 5 | 33.0 | 43 | 1.0 | 7 |
| EGFR | symmetric_gene_rank | 0.9983 | 0.9048 | 1.0 | 6 | 3 | 0 | 5 | 3.0 | 4 | 1.0 | 6 |
| EGFR | uniform_ratio_gene_rank | 0.9960 | 1.0000 | 1.0 | 18 | 6 | 1 | 5 | 9.0 | 18 | 1.0 | 5 |
| FRS2 | column_gene_rank | 0.9985 | 1.0000 | 1.0 | 5 | 3 | 0 | 0 | NA | NA | 1.0 | 5 |
| FRS2 | column_minmax | 0.9998 | 1.0000 | 0.0 | 4 | 0 | 0 | 0 | NA | NA | 0.0 | 4 |
| FRS2 | symmetric_gene_rank | 0.9978 | 0.9048 | 1.0 | 12 | 3 | 1 | 0 | NA | NA | 1.0 | 12 |
| FRS2 | uniform_ratio_gene_rank | 0.9961 | 1.0000 | 1.0 | 12 | 8 | 2 | 0 | NA | NA | 1.0 | 12 |
| HSD3B2 | column_gene_rank | 0.9998 | 1.0000 | 0.0 | 4 | 0 | 0 | 0 | NA | NA | 0.0 | 4 |
| HSD3B2 | column_minmax | 1.0000 | 1.0000 | 0.0 | 2 | 0 | 0 | 0 | NA | NA | 0.0 | 2 |
| HSD3B2 | symmetric_gene_rank | 0.9995 | 1.0000 | 0.0 | 4 | 0 | 0 | 0 | NA | NA | 0.0 | 4 |
| HSD3B2 | uniform_ratio_gene_rank | 0.9996 | 1.0000 | 0.0 | 4 | 0 | 0 | 0 | NA | NA | 0.0 | 4 |
| IGF1R | column_gene_rank | 0.9976 | 0.9048 | 1.0 | 7 | 8 | 0 | 0 | NA | NA | 1.0 | 7 |
| IGF1R | column_minmax | 0.9998 | 1.0000 | 0.0 | 3 | 0 | 0 | 0 | NA | NA | 0.0 | 3 |
| IGF1R | symmetric_gene_rank | 0.9984 | 1.0000 | 1.0 | 6 | 5 | 0 | 0 | NA | NA | 1.0 | 6 |
| IGF1R | uniform_ratio_gene_rank | 0.9977 | 1.0000 | 1.0 | 8 | 7 | 0 | 0 | NA | NA | 1.0 | 8 |
| IGF2 | column_gene_rank | 0.9967 | 0.9048 | 1.0 | 13 | 8 | 1 | 0 | NA | NA | 1.0 | 13 |
| IGF2 | column_minmax | 0.9997 | 1.0000 | 0.0 | 3 | 0 | 0 | 0 | NA | NA | 0.0 | 3 |
| IGF2 | symmetric_gene_rank | 0.9974 | 0.8182 | 1.0 | 13 | 6 | 1 | 0 | NA | NA | 1.0 | 13 |
| IGF2 | uniform_ratio_gene_rank | 0.9963 | 0.9048 | 1.0 | 13 | 6 | 1 | 0 | NA | NA | 1.0 | 13 |
| IL7R | column_gene_rank | 0.9960 | 0.9048 | 1.0 | 18 | 8 | 1 | 1 | 1.0 | 1 | 1.0 | 18 |
| IL7R | column_minmax | 0.9987 | 1.0000 | 0.0 | 14 | 1 | 1 | 1 | 14.0 | 14 | 0.0 | 4 |
| IL7R | symmetric_gene_rank | 0.9976 | 0.9048 | 1.0 | 7 | 5 | 0 | 1 | 1.0 | 1 | 1.0 | 7 |
| IL7R | uniform_ratio_gene_rank | 0.9958 | 0.9048 | 1.0 | 13 | 12 | 2 | 1 | 0.0 | 0 | 1.0 | 13 |
| KDR | column_gene_rank | 0.9983 | 1.0000 | 0.0 | 8 | 4 | 0 | 1 | 1.0 | 1 | 0.0 | 8 |
| KDR | column_minmax | 0.9992 | 1.0000 | 0.0 | 10 | 1 | 1 | 1 | 10.0 | 10 | 0.0 | 3 |
| KDR | symmetric_gene_rank | 0.9990 | 1.0000 | 0.0 | 5 | 1 | 0 | 1 | 1.0 | 1 | 0.0 | 5 |
| KDR | uniform_ratio_gene_rank | 0.9982 | 1.0000 | 1.0 | 7 | 6 | 0 | 1 | 4.0 | 4 | 1.0 | 7 |
| KRAS | column_gene_rank | 0.9986 | 1.0000 | 1.0 | 6 | 4 | 0 | 0 | NA | NA | 1.0 | 6 |
| KRAS | column_minmax | 0.9999 | 1.0000 | 0.0 | 3 | 0 | 0 | 0 | NA | NA | 0.0 | 3 |
| KRAS | symmetric_gene_rank | 0.9966 | 0.8182 | 1.0 | 14 | 9 | 1 | 0 | NA | NA | 1.0 | 14 |
| KRAS | uniform_ratio_gene_rank | 0.9973 | 0.9048 | 1.0 | 9 | 7 | 0 | 0 | NA | NA | 1.0 | 9 |
| LRP1B | column_gene_rank | 0.9968 | 0.9048 | 1.0 | 11 | 9 | 1 | 0 | NA | NA | 1.0 | 11 |
| LRP1B | column_minmax | 0.9998 | 1.0000 | 0.0 | 3 | 0 | 0 | 0 | NA | NA | 0.0 | 3 |
| LRP1B | symmetric_gene_rank | 0.9982 | 0.8182 | 1.0 | 8 | 3 | 0 | 0 | NA | NA | 1.0 | 8 |
| LRP1B | uniform_ratio_gene_rank | 0.9946 | 0.8182 | 1.0 | 13 | 13 | 2 | 0 | NA | NA | 1.0 | 13 |
| MED12 | column_gene_rank | 0.9976 | 0.9048 | 1.0 | 12 | 3 | 1 | 1 | 12.0 | 12 | 1.0 | 8 |
| MED12 | column_minmax | 0.9932 | 1.0000 | 0.0 | 36 | 1 | 1 | 1 | 36.0 | 36 | 0.0 | 4 |
| MED12 | symmetric_gene_rank | 0.9945 | 0.9048 | 1.0 | 20 | 9 | 2 | 1 | 20.0 | 20 | 1.0 | 15 |
| MED12 | uniform_ratio_gene_rank | 0.9944 | 1.0000 | 1.0 | 13 | 13 | 3 | 1 | 13.0 | 13 | 1.0 | 13 |
| MEN1 | column_gene_rank | 0.9945 | 1.0000 | 1.0 | 13 | 13 | 4 | 7 | 3.0 | 7 | 1.0 | 13 |
| MEN1 | column_minmax | 0.9778 | 0.8182 | 1.0 | 49 | 9 | 2 | 7 | 4.0 | 49 | 1.0 | 8 |
| MEN1 | symmetric_gene_rank | 0.9923 | 0.9048 | 1.0 | 20 | 15 | 3 | 7 | 1.0 | 5 | 1.0 | 20 |
| MEN1 | uniform_ratio_gene_rank | 0.9877 | 0.8182 | 2.0 | 25 | 21 | 5 | 7 | 2.0 | 10 | 2.0 | 25 |
| MLH1 | column_gene_rank | 0.9998 | 0.9048 | 0.0 | 2 | 0 | 0 | 0 | NA | NA | 0.0 | 2 |
| MLH1 | column_minmax | 0.9998 | 1.0000 | 0.0 | 3 | 0 | 0 | 0 | NA | NA | 0.0 | 3 |
| MLH1 | symmetric_gene_rank | 0.9992 | 0.9048 | 0.0 | 5 | 2 | 0 | 0 | NA | NA | 0.0 | 5 |
| MLH1 | uniform_ratio_gene_rank | 0.9992 | 0.9048 | 1.0 | 4 | 0 | 0 | 0 | NA | NA | 1.0 | 4 |
| MSH2 | column_gene_rank | 0.9990 | 0.9048 | 0.0 | 6 | 2 | 0 | 4 | 2.0 | 6 | 0.0 | 5 |
| MSH2 | column_minmax | 0.9991 | 0.9048 | 0.5 | 5 | 1 | 0 | 4 | 2.0 | 4 | 0.0 | 5 |
| MSH2 | symmetric_gene_rank | 0.9990 | 0.9048 | 1.0 | 5 | 2 | 0 | 4 | 1.5 | 2 | 1.0 | 5 |
| MSH2 | uniform_ratio_gene_rank | 0.9993 | 1.0000 | 0.0 | 5 | 1 | 0 | 4 | 0.5 | 3 | 0.0 | 5 |
| MSH6 | column_gene_rank | 0.9994 | 0.9048 | 0.0 | 4 | 0 | 0 | 3 | 0.0 | 2 | 0.0 | 4 |
| MSH6 | column_minmax | 0.9862 | 1.0000 | 1.0 | 51 | 3 | 1 | 3 | 6.0 | 51 | 1.0 | 5 |
| MSH6 | symmetric_gene_rank | 0.9989 | 0.9048 | 1.0 | 6 | 2 | 0 | 3 | 1.0 | 4 | 1.0 | 6 |
| MSH6 | uniform_ratio_gene_rank | 0.9991 | 0.9048 | 1.0 | 5 | 1 | 0 | 3 | 0.0 | 1 | 1.0 | 5 |
| MTOR | column_gene_rank | 0.9991 | 1.0000 | 0.0 | 4 | 0 | 0 | 0 | NA | NA | 0.0 | 4 |
| MTOR | column_minmax | 1.0000 | 1.0000 | 0.0 | 2 | 0 | 0 | 0 | NA | NA | 0.0 | 2 |
| MTOR | symmetric_gene_rank | 0.9987 | 0.9048 | 1.0 | 7 | 2 | 0 | 0 | NA | NA | 1.0 | 7 |
| MTOR | uniform_ratio_gene_rank | 0.9985 | 1.0000 | 1.0 | 10 | 2 | 1 | 0 | NA | NA | 1.0 | 10 |
| NF1 | column_gene_rank | 0.9992 | 1.0000 | 0.0 | 5 | 1 | 0 | 0 | NA | NA | 0.0 | 5 |
| NF1 | column_minmax | 0.9998 | 1.0000 | 0.0 | 3 | 0 | 0 | 0 | NA | NA | 0.0 | 3 |
| NF1 | symmetric_gene_rank | 0.9979 | 0.9048 | 1.0 | 9 | 5 | 0 | 0 | NA | NA | 1.0 | 9 |
| NF1 | uniform_ratio_gene_rank | 0.9977 | 1.0000 | 1.0 | 8 | 5 | 0 | 0 | NA | NA | 1.0 | 8 |
| NR5A1 | column_gene_rank | 0.9992 | 1.0000 | 0.0 | 5 | 2 | 0 | 0 | NA | NA | 0.0 | 5 |
| NR5A1 | column_minmax | 0.9999 | 1.0000 | 0.0 | 2 | 0 | 0 | 0 | NA | NA | 0.0 | 2 |
| NR5A1 | symmetric_gene_rank | 0.9990 | 1.0000 | 0.0 | 6 | 4 | 0 | 0 | NA | NA | 0.0 | 6 |
| NR5A1 | uniform_ratio_gene_rank | 0.9978 | 1.0000 | 1.0 | 10 | 5 | 1 | 0 | NA | NA | 1.0 | 10 |
| PIK3CA | column_gene_rank | 0.9982 | 0.9048 | 1.0 | 8 | 4 | 0 | 0 | NA | NA | 1.0 | 8 |
| PIK3CA | column_minmax | 0.9999 | 1.0000 | 0.0 | 2 | 0 | 0 | 0 | NA | NA | 0.0 | 2 |
| PIK3CA | symmetric_gene_rank | 0.9972 | 0.9048 | 1.0 | 16 | 4 | 1 | 0 | NA | NA | 1.0 | 16 |
| PIK3CA | uniform_ratio_gene_rank | 0.9976 | 1.0000 | 1.0 | 11 | 5 | 1 | 0 | NA | NA | 1.0 | 11 |
| PMS2 | column_gene_rank | 0.9995 | 0.9048 | 0.0 | 5 | 1 | 0 | 0 | NA | NA | 0.0 | 5 |
| PMS2 | column_minmax | 0.9998 | 1.0000 | 0.0 | 3 | 0 | 0 | 0 | NA | NA | 0.0 | 3 |
| PMS2 | symmetric_gene_rank | 0.9992 | 0.9048 | 0.0 | 5 | 2 | 0 | 0 | NA | NA | 0.0 | 5 |
| PMS2 | uniform_ratio_gene_rank | 0.9994 | 1.0000 | 1.0 | 4 | 0 | 0 | 0 | NA | NA | 1.0 | 4 |
| PRKAR1A | column_gene_rank | 0.9948 | 0.9048 | 1.0 | 19 | 12 | 1 | 0 | NA | NA | 1.0 | 19 |
| PRKAR1A | column_minmax | 0.9998 | 1.0000 | 0.0 | 4 | 0 | 0 | 0 | NA | NA | 0.0 | 4 |
| PRKAR1A | symmetric_gene_rank | 0.9936 | 1.0000 | 1.0 | 19 | 13 | 2 | 0 | NA | NA | 1.0 | 19 |
| PRKAR1A | uniform_ratio_gene_rank | 0.9902 | 0.8182 | 1.0 | 20 | 15 | 8 | 0 | NA | NA | 1.0 | 20 |
| PTCH1 | column_gene_rank | 0.9981 | 0.9048 | 1.0 | 8 | 4 | 0 | 3 | 1.0 | 8 | 1.0 | 6 |
| PTCH1 | column_minmax | 0.9873 | 1.0000 | 1.0 | 30 | 6 | 3 | 3 | 29.0 | 30 | 1.0 | 7 |
| PTCH1 | symmetric_gene_rank | 0.9990 | 1.0000 | 1.0 | 6 | 2 | 0 | 3 | 2.0 | 3 | 1.0 | 6 |
| PTCH1 | uniform_ratio_gene_rank | 0.9985 | 1.0000 | 0.0 | 9 | 3 | 0 | 3 | 1.0 | 4 | 0.0 | 9 |
| RB1 | column_gene_rank | 0.9960 | 0.9048 | 1.0 | 11 | 10 | 2 | 2 | 5.5 | 8 | 1.0 | 11 |
| RB1 | column_minmax | 0.9731 | 0.9048 | 1.0 | 57 | 3 | 2 | 2 | 51.5 | 57 | 1.0 | 7 |
| RB1 | symmetric_gene_rank | 0.9923 | 1.0000 | 1.0 | 19 | 13 | 4 | 2 | 6.0 | 9 | 1.0 | 19 |
| RB1 | uniform_ratio_gene_rank | 0.9933 | 0.9048 | 1.5 | 14 | 16 | 5 | 2 | 5.0 | 7 | 1.0 | 14 |
| RPL22 | column_gene_rank | 0.9984 | 0.9048 | 0.0 | 15 | 2 | 1 | 0 | NA | NA | 0.0 | 15 |
| RPL22 | column_minmax | 1.0000 | 1.0000 | 0.0 | 2 | 0 | 0 | 0 | NA | NA | 0.0 | 2 |
| RPL22 | symmetric_gene_rank | 0.9992 | 1.0000 | 0.0 | 6 | 2 | 0 | 0 | NA | NA | 0.0 | 6 |
| RPL22 | uniform_ratio_gene_rank | 0.9990 | 1.0000 | 1.0 | 8 | 1 | 0 | 0 | NA | NA | 1.0 | 8 |
| STAR | column_gene_rank | 0.9996 | 1.0000 | 0.0 | 4 | 0 | 0 | 0 | NA | NA | 0.0 | 4 |
| STAR | column_minmax | 1.0000 | 1.0000 | 0.0 | 0 | 0 | 0 | 0 | NA | NA | 0.0 | 0 |
| STAR | symmetric_gene_rank | 0.9996 | 0.9048 | 0.0 | 6 | 1 | 0 | 0 | NA | NA | 0.0 | 6 |
| STAR | uniform_ratio_gene_rank | 0.9996 | 1.0000 | 0.0 | 3 | 0 | 0 | 0 | NA | NA | 0.0 | 3 |
| TERF2 | column_gene_rank | 0.9974 | 0.8182 | 0.0 | 13 | 6 | 1 | 0 | NA | NA | 0.0 | 13 |
| TERF2 | column_minmax | 0.9997 | 1.0000 | 0.0 | 3 | 0 | 0 | 0 | NA | NA | 0.0 | 3 |
| TERF2 | symmetric_gene_rank | 0.9987 | 1.0000 | 1.0 | 6 | 2 | 0 | 0 | NA | NA | 1.0 | 6 |
| TERF2 | uniform_ratio_gene_rank | 0.9985 | 1.0000 | 1.0 | 8 | 2 | 0 | 0 | NA | NA | 1.0 | 8 |
| TERT | column_gene_rank | 0.9992 | 0.9048 | 1.0 | 4 | 0 | 0 | 16 | 1.0 | 3 | 1.0 | 4 |
| TERT | column_minmax | 0.9167 | 0.7391 | 3.0 | 54 | 37 | 16 | 16 | 16.5 | 54 | 3.0 | 19 |
| TERT | symmetric_gene_rank | 0.9983 | 1.0000 | 1.0 | 6 | 3 | 0 | 16 | 0.0 | 2 | 1.0 | 6 |
| TERT | uniform_ratio_gene_rank | 0.9973 | 1.0000 | 1.0 | 10 | 8 | 1 | 16 | 2.0 | 10 | 1.0 | 8 |
| TP53 | column_gene_rank | 0.9994 | 1.0000 | 0.0 | 5 | 1 | 0 | 0 | NA | NA | 0.0 | 5 |
| TP53 | column_minmax | 0.9999 | 1.0000 | 0.0 | 3 | 0 | 0 | 0 | NA | NA | 0.0 | 3 |
| TP53 | symmetric_gene_rank | 0.9973 | 0.9048 | 1.0 | 8 | 8 | 0 | 0 | NA | NA | 1.0 | 8 |
| TP53 | uniform_ratio_gene_rank | 0.9980 | 1.0000 | 1.0 | 8 | 5 | 0 | 0 | NA | NA | 1.0 | 8 |
| VEGFA | column_gene_rank | 1.0000 | 1.0000 | 0.0 | 0 | 0 | 0 | 0 | NA | NA | 0.0 | 0 |
| VEGFA | column_minmax | 1.0000 | 1.0000 | 0.0 | 0 | 0 | 0 | 0 | NA | NA | 0.0 | 0 |
| VEGFA | symmetric_gene_rank | 1.0000 | 1.0000 | 0.0 | 0 | 0 | 0 | 0 | NA | NA | 0.0 | 0 |
| VEGFA | uniform_ratio_gene_rank | 1.0000 | 1.0000 | 0.0 | 0 | 0 | 0 | 0 | NA | NA | 0.0 | 0 |
| ZNRF3 | column_gene_rank | 0.9990 | 1.0000 | 1.0 | 6 | 1 | 0 | 0 | NA | NA | 1.0 | 6 |
| ZNRF3 | column_minmax | 0.9998 | 1.0000 | 0.0 | 3 | 0 | 0 | 0 | NA | NA | 0.0 | 3 |
| ZNRF3 | symmetric_gene_rank | 0.9991 | 1.0000 | 0.0 | 5 | 1 | 0 | 0 | NA | NA | 0.0 | 5 |
| ZNRF3 | uniform_ratio_gene_rank | 0.9990 | 1.0000 | 0.0 | 9 | 2 | 0 | 0 | NA | NA | 0.0 | 9 |

## Convergence and reproducibility

- Python: `3.14.3`; NumPy: `2.4.4`; SciPy: `1.17.1`.
- Column-stochastic batch: 30 iterations; final maximum L1 delta `6.996e-11`.
- Symmetric batch: 31 iterations; final maximum L1 delta `5.690e-11`.
- Uniform-reference run: 30 iterations; final maximum L1 delta `7.119e-11`.

### Input SHA-256

- `data/ACC_P0.5C_gene_weights_v1.csv`: `a513b96cf50ee10041068ab621ba978acfdf0c89116df32f86a9fcbc043a3897`
- `data/bindex_network/bindex_edges_1304.csv`: `0ac37e507e763b5d20e78571df372ac14350e95b4d19c0360d1324ed35018d63`
- `data/bindex_network/rACC_399_fullSTRING.csv`: `9c9a18d8937b83f1da79cbe8f8e6d6870c30b267ed410fbcc283edc0507c4ab8`
- `9606.protein.info.v12.0.txt.gz`: `144de4b0d98c6a7dfde6ddc2591cf88657f27b989eadff4f501450c3ed1f0f1c`
- `9606.protein.links.v12.0.txt.gz`: `3e22f32572211aa341d5b4bd08d30c32e693e294603202120936872f87719d4f`

## Interpretation boundary

This scan measures deterministic ranking sensitivity when one curated seed is removed. It does not test whether an omitted seed is biologically correct and does not provide external efficacy validation.
