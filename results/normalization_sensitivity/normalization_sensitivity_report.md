# Propagation-normalization sensitivity report

- Version: `normalization-sensitivity-v1-reviewer-m4`.
- Degree-matched draws per variant: 10,000.
- Best-case BH q floor: 0.010799.

| Variant | rho(gene, degree) | rho(drug, locked C_ACC) | Top-20 Jaccard | Drugs q<0.05 | Minimum q | CDK4/6 P | CDK4/6 q across variants |
|---|---:|---:|---:|---:|---:|---:|---:|
| column_minmax | 0.797 | 1.000 | 1.000 | 2 | 0.0108 | 0.2839 | 0.2839 |
| column_gene_rank | 0.797 | 0.622 | 0.379 | 0 | 0.0756 | 0.0104 | 0.0347 |
| uniform_ratio_gene_rank | 0.517 | 0.641 | 0.429 | 0 | 0.1080 | 0.0260 | 0.0347 |
| symmetric_gene_rank | 0.624 | 0.640 | 0.429 | 1 | 0.0432 | 0.0214 | 0.0347 |

The alternatives are sensitivity analyses, not outcome-optimized replacement models. A drug-level FDR signal means that its network context exceeded the matched-seed null under that transformation; it does not establish efficacy or clinical validity.
