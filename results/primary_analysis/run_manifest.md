# Primary analysis run manifest

## Locked model

- Model version: `primary-108-v2`
- Primary universe: `108` complete-case drugs
- C_ACC shrinkage pseudo-count: `3.0`
- Percentile method: average rank scaled to [0,1] within analysis universe
- Primary formula: `0.50*C_ACC_pct + 0.50*residual_pct`
- Evidence-informed formula: `0.40*C_ACC_pct + 0.40*residual_pct + 0.20*S_external`
- Clinical efficacy benchmark: retired after the C4 evidence-label audit because
  the strict subset contains two positive and zero negative clinical comparators.

## Environment

- Python: `3.14.3`
- NumPy: `2.4.4`
- SciPy: `1.17.1`
- Platform: `Windows-11-10.0.26200-SP0`
- Command: `python -m analysis.acc_primary_pipeline`

## Input SHA-256

- `data/bindex_network/bindex_edges_1304.csv`: `0ac37e507e763b5d20e78571df372ac14350e95b4d19c0360d1324ed35018d63`
- `data/bindex_network/rACC_399_fullSTRING.csv`: `9c9a18d8937b83f1da79cbe8f8e6d6870c30b267ed410fbcc283edc0507c4ab8`
- `data/bindex_network/Sactivity_124_v1.csv`: `5d371ff5d8b8261f1b3a091131506be67fb0a39f39563c333a0bcbd9f5ade475`
- `data/bindex_network/NCI60_potency_124.csv`: `60e5e60e5211990b659b807181d968cb3c06e9c0889635d12d5fea15fb2bb367`
- `data/bindex_network/S_external_curated.csv`: `31b01d479ef24d51d2b92de3345ddb6050fdd7bba2756fa5fde47ded6a1f35d5`

## Output SHA-256

- `ADRS_comp_primary_108.csv`: `6e0c45dbb193d1660862b91d540d236e85e35a0359405cb6718a71c8b8318c8f`
- `ADRS_evidence_informed_108.csv`: `bbfa47f51b159102bca9057b980eb54ef673368e3c74e3155222056171e9291f`
- `ADRS_context_only_16.csv`: `f34d12712964c7d94369f2026f102d6c2958ec762a7b307004b497e03c7a2c1a`
- `primary_metrics.json`: `947d276310a68244aff31207e82dfabf489abe3bcf736e4f4ae0fb8064645f57`

## Interpretation guardrail

No ROC-AUC or PR-AUC is emitted. The heterogeneous legacy labels are retained
only in the C4 evidence audit and are not a clinical-efficacy validation set.
