# C3 change report: model-matched figures

- Figure version: `c3-primary108-v1`.
- Target: Pharmaceutics/MDPI full width (170 mm).
- Figure 3 was rebuilt from the C1 108-drug table; legacy S_B-neighbor/external columns and pairwise 124-drug mixtures were removed.
- Figure 5a uses 21 two-component weights from w_C=0 to 1 in 0.05 increments.
- Figure 5b uses all 204,156 possible three-drug subsets, not a Monte-Carlo sample.
- CDK4/6 annotation is computed as ranks 8/26/51, mean 28.33, exact P=0.0764, BH q=0.3711.
- The title and caption state “non-significant ranking trend”; the legacy positive-enrichment wording was removed.
- Figure 4 was not regenerated because its evidence-label semantics are the subject of C4; retaining a temporary figure gap is preferable to presenting the mixed-evidence set as clinical validation.

## Files

- `figure_data/revision/Fig3_component_correlation_primary108.csv`
- `figure_data/revision/Fig5a_weight_scan_primary108.csv`
- `figure_data/revision/Fig5b_CDK46_exact_null_primary108.csv`
- `figure_data/revision/C3_figure_stats.json`
- `figures/revision/Fig3_component_correlation_primary108.pdf`
- `figures/revision/Fig3_component_correlation_primary108.svg`
- `figures/revision/Fig3_component_correlation_primary108.png`
- `figures/revision/Fig5_weight_stability_CDK46_primary108.pdf`
- `figures/revision/Fig5_weight_stability_CDK46_primary108.svg`
- `figures/revision/Fig5_weight_stability_CDK46_primary108.png`
- `results/figures/C3_figure_QA.json`
- `projects/ACC-PHARMA-NET/figures/manifest.md`

## QA

- F3: width 169.977 mm; PNG dpi 1000.0; minimum font 8.0 pt.
- F5: width 169.977 mm; PNG dpi 1000.0; minimum font 8.0 pt.
- Vector PDF and SVG retain editable text.
