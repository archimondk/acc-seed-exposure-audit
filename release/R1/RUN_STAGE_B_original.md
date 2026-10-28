# 阶段B本地运行说明

协议：melanoma_R1_v0.3。阶段B已于2026-10-18T22:09:43-04:00冻结（`04_melanoma_protocol/FREEZE_STAGE_B.txt`）。传播需要STRING v12图，只能在本地运行。

## 0. 先打tag

把`04_melanoma_protocol/FREEZE_STAGE_B.txt`、`STAGE_B_PREDICTIONS.md`、`scripts/stageB_exposure_and_anchors.py`和`results/melanoma_stageB_prep/`下的四个文件提交，打tag `melanoma-protocol-stageB`。**tag发布之后再运行第3步以后的内容。**

## 1. 放置代码

| 文件 | 放到ACC仓库的位置 |
|---|---|
| `analysis/melanoma_rules.py` | `analysis/` |
| `analysis/melanoma_audit.py` | `analysis/` |
| `tests/test_melanoma_rules.py`、`tests/test_melanoma_audit.py` | `tests/` |

`--project`指向返修目录，即同时包含`04_melanoma_protocol/`和`results/`的那个目录。阶段B冻结记录里的路径都相对于这个目录。程序启动时会依次核对阶段A、A1附录和阶段B三份冻结记录中所有文件的哈希，任何一个对不上都会拒绝运行。

## 2. 测试

```
pytest tests/test_melanoma_rules.py tests/test_melanoma_audit.py
```

预期42项全部通过（规则32项、合成网络端到端10项）。

## 3. 乳腺癌回归核验（真实STRING）

```
python -m analysis.melanoma_audit --validate-breast --string-dir <STRING目录> --out results/melanoma_stageB/validate
```

用本程序的引擎重算乳腺癌在108药、四种归一化下的观测百分位，与存档`positive_control_drug_ranks_primary108.csv`逐行比对。预期为`PASS`，最大差异不超过5e-11（存档值只保留12位有效数字）。**不是PASS就停下，不要继续。**

本环境做过同样的比对：在一个包含全部399个真实关联基因的合成STRING图上，同时运行原仓库的`run_positive_control_analysis`和本程序，432行结果完全相同，最大差异为0。

## 4. 计时

```
python -m analysis.melanoma_audit --project <返修目录> --string-dir <STRING目录> --out results/melanoma_stageB/benchmark --benchmark
```

只跑100次随机种子（arm_index 99，不属于正式随机流），不产生科学输出。`projected_with_30pct_margin`超过12小时的，按方案第6节先停下报告，不自行删减分支或减少B。

## 5. 正式运行

```
python -m analysis.melanoma_audit --project <返修目录> --string-dir <STRING目录> --out results/melanoma_stageB/run
```

主要输出：

| 文件 | 内容 |
|---|---|
| `decisions.json` | H4（124；108为S1）、H5（108；带算术标注）、H5b、H6（124；108另报） |
| `drug_statistics.csv` | 两臂 × 两个药物集合 × 四种归一化 × 每个药：排名、百分位、零模型均值和SD、z、超越次数、P、Clopper–Pearson区间、组内BH |
| `audit_metrics.csv` | 第4.1节审计指标；M0的行即跨疾病对照表中黑色素瘤的部分 |
| `seed_nonseed_decomposition.csv` | 种子/非种子分解，逐药 |
| `H6_delta_z_124.csv`、`H6_delta_z_108.csv` | M1−M0的Δz |
| `pooled_bh_992.csv` | 两臂 × 四变体 × 124药合并BH（敏感性） |
| `loso_M0.csv` | M0逐个剔除种子 |
| `S10_parameter_sensitivity.csv` | α、k、STRING阈值700 |
| `seed_branches.csv` | 阶段A各种子分支的确定性排序 |
| `null_seed_sets_M*.csv`、`null_C_raw_M*.npz`、`observed_C_raw_M*.npz` | 全部零模型抽样与全精度分数 |
| `run_manifest.json` | 冻结核验结果、代码和STRING哈希、图规模、迭代次数、软件版本、耗时 |

跑完后把`decisions.json`、`audit_metrics.csv`、`run_manifest.json`和`drug_statistics.csv`发给我。

## 实现说明（不改变冻结规则）

1. **零方差判定**：方案要求零模型标准差严格大于0时才计算z，等于0时记为NA，不加ε。实现中改为直接检查零模型的值是否完全相同，因为对一个完全恒定的向量，浮点运算算出的标准差可能是约1e-17而不是0。这样仍然不使用任何容差，而且避免了产生极大的虚假z值。
2. **五个随机流**：通过向原有`generate_degree_matched_seed_sets`传入`SeedSequence([stream_seed, arm_index])`实现，抽样函数本身没有改动。
3. **float64**：零模型百分位全程使用float64；乳腺癌当时存的是float32。
4. **H5的124药判定**：需要ACC和乳腺癌在124药集合中的参照点（属于第5项补算），本次输出`PENDING_124_ANCHORS`，只报告黑色素瘤的数值。
5. **STRING ≥ 700**：S10中的阈值分支会重新载入一次图，运行时间相应增加。
