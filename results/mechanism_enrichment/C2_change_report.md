# C2 change report: mechanism-class enrichment

## Locked analysis

- Version: `mechanism-primary108-v1`.
- Universe: the 108 drugs in the C1 complete-case ADRS_comp ranking.
- Class labels: MIPE `Primary MOA` plus explicit, score-blind drug-name aliases.
- Eligible tests: frozen named classes represented by at least 3 drugs; `Other` is not tested.
- Statistic: within-class rank sum (equivalently mean rank; lower is better).
- Null: every same-size subset of ranks 1–108, counted exactly by dynamic programming.
- Multiplicity: Benjamini–Hochberg across the ten eligible tests.

## Corrected results

| Class | k | Mean rank | Exact P | BH q |
|---|---:|---:|---:|---:|
| ALK/MET | 4 | 30.00 | 0.0576 | 0.3711 |
| CDK4/6 | 3 | 28.33 | 0.0764 | 0.3711 |
| HDAC | 4 | 35.25 | 0.1113 | 0.3711 |
| Antimetabolite | 13 | 49.69 | 0.2812 | 0.6565 |
| Tubulin | 5 | 50.20 | 0.3829 | 0.6565 |
| Topo/anthracycline | 11 | 53.55 | 0.4601 | 0.6565 |
| MEK | 4 | 53.50 | 0.4781 | 0.6565 |
| EGFR/HER | 8 | 55.12 | 0.5252 | 0.6565 |
| Alkylator/platinum | 12 | 61.75 | 0.8015 | 0.8559 |
| Multikinase/VEGFR | 6 | 67.83 | 0.8559 | 0.8559 |

## CDK4/6 correction

The primary universe contains Abemaciclib, Palbociclib and Ribociclib at ranks 8; 26; 51 (mean 28.33). The exact one-sided random-set test gives P=0.0764 and BH q=0.3711.

The old P=0.023/FDR=0.10 statement is not reproducible under the locked C1 model and must be removed. The supported wording is: **“The pre-specified CDK4/6 class showed a non-significant ranking trend.”** This is a candidate hypothesis, not evidence of efficacy.

Trilaciclib is correctly absent because it lacks MIPE data and is one of the 16 context-only drugs excluded by the C1 complete-case rule.

## Classification audit

Two legacy substring rules were removed: `"alk" in MOA` could misclassify “alkaloid” as ALK, and `"egfr" in MOA` also matched the substring in “VEGFR”. Separating VEGFR from EGFR yields ten, rather than nine, eligible families. All name aliases are now listed in the analysis source, and the complete 108-drug membership table records the source and assignment rule for every drug.

## Downstream synchronization

C3 must regenerate Figure 5b from the CSV output and replace the hard-coded P=0.023 annotation. Manuscript-wide wording and figure legends will be synchronized after the C3 figure outputs are frozen.
