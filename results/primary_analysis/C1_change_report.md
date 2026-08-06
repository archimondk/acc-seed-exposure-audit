# C1 change report: locked primary model

## Outcome

The authoritative primary analysis is frozen as `primary-108-v2`:

- 108 drugs with both MIPE and NCI-60 measurements;
- `C_ACC` background-mean shrinkage with `k = 3`;
- average-rank percentiles within the 108-drug universe;
- `ADRS_comp = 0.50 × P(C_ACC) + 0.50 × P(residual)`;
- no secondary `n_assoc` confidence transform;
- 16 MIPE-missing drugs excluded from primary/EI rankings and reported
  separately.

## Corrected results

| Result | Legacy manuscript/output | Locked primary-108-v2 |
|---|---:|---:|
| Mixed-evidence benchmark | AUC 0.40 | Retired by C4; strict clinical AUC is not estimable |
| Mitotane primary rank | 23 | 17 |
| Doxorubicin primary rank | 54 | 23 |
| Abemaciclib primary rank | 14 | 8 |
| Ixazomib primary rank | 6 | 3 |
| Palbociclib primary rank | 31 | 26 |
| Ribociclib primary rank | 47 | 51 |

C4 subsequently separated clinical, regimen/case, preclinical and
class-extrapolated evidence. The strict drug-specific clinical subset has two
positive candidates and no unambiguous negative comparator; therefore ROC-AUC
and PR-AUC are not estimable. The old mixed-evidence AUC is retired and is no
longer emitted by the C1 pipeline.

## Authoritative artifacts

- `analysis/acc_primary_pipeline.py`
- `tests/test_acc_primary_pipeline.py`
- `results/primary_analysis/ADRS_comp_primary_108.csv`
- `results/primary_analysis/ADRS_evidence_informed_108.csv`
- `results/primary_analysis/ADRS_context_only_16.csv`
- `results/primary_analysis/primary_metrics.json`
- `results/primary_analysis/run_manifest.md`

## Verification

- `python -m pytest -q tests/test_acc_primary_pipeline.py` -> 3 passed.
- `python -m py_compile analysis/acc_primary_pipeline.py` -> exit 0.
- `python -m analysis.acc_primary_pipeline` -> status `ok`, `n_primary=108`,
  clinical benchmark `retired_in_C4`.

## Downstream status

C2 mechanism enrichment, C3 figure reconstruction and C4 evidence-label audit
are complete. Manuscript text synchronization remains a later revision stage.
