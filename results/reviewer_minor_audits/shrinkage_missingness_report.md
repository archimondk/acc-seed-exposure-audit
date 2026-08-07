# Shrinkage and MIPE-missingness audit

## Pseudo-count sensitivity

| k | Spearman vs k=3 | Top-20 Jaccard | Median absolute rank shift | Maximum absolute rank shift |
|---:|---:|---:|---:|---:|
| 1 | 0.979 | 1.000 | 2.0 | 19 |
| 3 | 1.000 | 1.000 | 0.0 | 0 |
| 5 | 0.998 | 1.000 | 1.0 | 7 |
| 10 | 0.985 | 1.000 | 3.0 | 20 |

## MIPE missingness

- Missing drugs: 16/124.
- Targeted-class missingness Fisher exact odds ratio: 3.627; two-sided P=0.0237.
- Median all-124 C_ACC rank: missing 72.0, observed 60.5; Mann-Whitney P=0.4230.

These descriptive tests do not establish a missingness mechanism. The 16-drug table is retained so readers can assess the excluded high-context candidates directly.
