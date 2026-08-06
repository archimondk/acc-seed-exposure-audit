# Post-hoc W2–MIPE external-criterion audit

Date: 2026-07-29

## Status and estimand

This diagnostic was proposed and executed after the Amendment 3 W1/W2 results
had been inspected. It is therefore post hoc and result-known, is not a new
frozen amendment, and does not revise `PARTIAL_OR_NOT_SUPPORTED`.

The question is narrower than ranking sensitivity: does the locked curated
weight-to-gene assignment show a stronger monotonic association with measured
MIPE potency than assignments obtained by randomly permuting the same weights
over the same 45 seed genes?

For each W2 draw \(p\), the diagnostic calculated

\[
\rho_p=\operatorname{Spearman}(C_{\mathrm{ACC},p},
P_{\mathrm{MIPE\ potency}})
\]

across the locked 108-drug complete-case universe. The locked statistic used
the corresponding frozen `C_ACC_pct` and `MIPE_potency_pct` columns. Higher
`MIPE_potency_pct` denotes greater measured potency, so a larger positive
correlation is the direction that could favour the curated assignment.

The locked percentile was defined descriptively as

\[
100\times \#\{\rho_p\leq\rho_{\mathrm{locked}}\}/1000.
\]

No threshold, equivalence margin, P value or FDR claim was prespecified.
For descriptive power context only, the locked-rho interval used
\(z=\operatorname{atanh}(\rho)\), standard error \(1/\sqrt{108-3}\), and
back-transformation of \(z\pm1.959964\,SE\).

## Inputs

- `results/seed_weight_assignment_sensitivity/permutation_drug_rank_draws.csv`
  - 108,000 rows: 1,000 draws × 108 drugs
  - SHA-256:
    `386fb9c24a911ab33edd2ae3514ad750d948ede80e857725d7af401a611097b3`
- `results/primary_analysis/ADRS_comp_primary_108.csv`
  - 108 locked complete-case drugs
  - SHA-256:
    `6e0c45dbb193d1660862b91d540d236e85e35a0359405cb6718a71c8b8318c8f`

Drug names were stripped, lower-cased and matched one-to-one. All 108,000 W2
rows mapped to a non-missing locked MIPE potency percentile.

## Result

| Quantity | Value |
|---|---:|
| Locked Spearman rho | 0.0607098455 |
| Approximate Fisher z 95% CI for locked rho (n = 108) | −0.1297526959 to +0.2468517626 |
| W2 minimum | 0.0080073166 |
| W2 5th percentile | 0.0325075143 |
| W2 median | 0.0629629453 |
| W2 95th percentile | 0.0913400942 |
| W2 maximum | 0.1191094302 |
| W2 mean (SD) | 0.0620354018 (0.0178771808) |
| W2 values below the locked rho | 458/1,000 |
| W2 values at or below the locked rho | 458/1,000 |
| Locked percentile in W2 distribution | 45.8th |
| W2 values at or above the locked rho | 542/1,000 |

The locked assignment is near the centre of the permutation distribution and
does not show a detectable advantage over random weight-to-gene assignments
with respect to measured MIPE potency. This does not prove that the curated
weights contain no biological information; it shows that the available MIPE
criterion does not distinguish the locked assignment from the W2 assignments.

Accordingly, W2 establishes sensitivity of ranks to weight-to-gene assignment,
not biological informativeness of the curated assignment. In the absence of a
discriminating external criterion, this sensitivity is best treated as an
investigator degree of freedom.

The permutation distribution is itself centred near zero: its median rho is
0.0630 and its 95th percentile is 0.0913. The interval around the locked rho is
nevertheless wide enough that a small phenotype association cannot be
excluded. This diagnostic should therefore not be described as an equivalence
test: it makes a large phenotype correlation for the curated assignment
unlikely but is underpowered to rule out a small one.

## Exact-rho coincidence check

The Amendment 2 minimum rank correlation and the Amendment 3 W1 correlation
are both

`0.9915497251517144`

up to the precision recoverable from their stored outputs. This is not a
copied value:

- the Amendment 2 minimum occurred at draw 253;
- the W1 and draw-253 rank vectors differed for 85 of 108 drugs;
- each nevertheless had the same sum of squared discrepancies from the locked
  ranks: 1,774;
- because both comparisons are permutations of ranks 1–108 without ties, the
  Spearman formula gives the same rho.

Spearman correlation is discrete in this setting because it is computed from
integer rank permutations. The matching value is therefore a verified
combinatorial coincidence rather than a transcription error.
