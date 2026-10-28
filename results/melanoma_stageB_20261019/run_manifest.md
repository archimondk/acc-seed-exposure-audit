# 黑色素瘤阶段 B 运行核验

日期：2026-10-19（Asia/Singapore）。协议：`melanoma_R1_v0.3`，阶段 A 与 A1 加上阶段 B 冻结记录。

## 顺序门槛

1. 阶段 B GitHub Release 已发布：<https://github.com/archimondk/acc-seed-exposure-audit/releases/tag/melanoma-protocol-stageB>。tag `melanoma-protocol-stageB` 指向提交 `bdd0c7030d0cb024427dfa4d7e5faa2a5e92ba23`；tag 下 7 个文件的 Git blob 均与本地冻结材料一致。发布时间 `2026-10-19T02:44:36Z`。发布早于本机正式分析开始。
2. `python -m pytest tests/test_melanoma_rules.py tests/test_melanoma_audit.py -q -p no:cacheprovider`：**42 passed**。首次尝试用 `python -I -m pytest`，因隔离模式不把本地项目目录加入模块路径而在收集阶段报 `No module named analysis`；改按包内说明运行后全部通过，没有修改测试或分析脚本。
3. 用本地真实 STRING v12 图进行乳腺癌回归核验：**PASS**，432 行，最大绝对差 `4.953260024365136e-13`，小于阈值 `5e-11`。文件：`validate/breast_regression.csv`。
4. `--benchmark`：100 次抽样耗时 `5.412989377975464` 秒；两臂各 10,000 次零模型且加 30% 余量的估算为 `1407.3772382736206` 秒，约 23.5 分钟，低于 12 小时门槛。该估算主要覆盖零模型；正式运行全部计算实际耗时 `1075.7169790267944` 秒，约 17.9 分钟。
5. 正式运行退出码 0，要求的四个输出均存在且能解析。两个分析臂各有 10,000 次零模型抽样；`drug_statistics.csv` 含 1,856 行，覆盖 M0/M1、108/124 药和四种归一化；`audit_metrics.csv` 含 16 行。程序在运行前核对阶段 A、A1、阶段 B 冻结哈希，分别为 3、1、6 个条目。

## 主要判定及边界

| 预设判定 | 本次输出 |
|---|---|
| H4 | 124 药与 108 药均 `fail`。124 药主归一化六药组检验 `p=q=0.9653034696530347`，六药没有进入前 25 名。 |
| H5（108 药） | `SUPPORTED`，但必须同时标注 `ARITHMETICALLY_DETERMINED`。本例覆盖度极低，判定由暴露数的算术界限决定，**不能作为种子暴露假说的经验支持证据**。 |
| H5（124 药） | `PENDING_124_ANCHORS`，参照疾病的 124 药锚点尚待补算。 |
| H6 | 124 药与 108 药均 `PARTIAL`，不是全面支持。 |
| H5b | `NOT_EVALUABLE`，六个 BRAF/MEK 药都未直接暴露。 |

上述是冻结规则和程序输出的判定；尚未修改投稿正文、补充材料或回复信。S4、S5、S6 不在收到的阶段 B 实现内，未标记为已完成。

## 环境与输入

- Python `3.14.3`，NumPy `2.4.4`，SciPy `1.17.1`，pytest `9.0.3`。
- STRING info SHA-256：`144de4b0d98c6a7dfde6ddc2591cf88657f27b989eadff4f501450c3ed1f0f1c`。
- STRING links SHA-256：`3e22f32572211aa341d5b4bd08d30c32e693e294603202120936872f87719d4f`。两者与原项目 `reproducibility/input_manifest.json` 一致。
- 本次运行图规模：19,700 节点，主阈值图读入 1,858,944 条边记录。
- 代码直接复制自用户提供的 `C:/Users/008bu/Downloads/3/melanoma_stageB_package.zip`，复制后逐文件哈希一致。
- 为满足脚本的 `--project` 目录结构，冻结文件与阶段 A 已映射种子表按原字节复制到 `manuscript_revision/plos_one_R1_20261018/stageB_runtime_20261019/`；原有 v0.1 草案和 v0.3 冻结源文件保留。

## 19,701 与 19,700 节点的差异核验

本地 STRING v12 `protein.info` 的 `preferred_name` 有 **19,699** 个唯一符号。`load_string_graph` 将这些符号与调用方的 `required_nodes` 做并集，而非只保留原生 STRING 节点。

- B-index 的 399 个关联基因中，`MALAT1` 不在上述 STRING 符号集中；ACC 与黑色素瘤两次调用均将它作为必需关联基因加入，成为度数为 0 的孤立节点。原 ACC 的 `results/method_strengthening/centrality_gene399.csv` 也记录其度数为 0。
- ACC 的 45 个正权重种子中，`VEGFA` 也不在 STRING 符号集中；原 ACC 调用将其额外加入，形成第二个度数为 0 的孤立节点。因此 ACC 为 **19,699 + 2 = 19,701** 节点。
- 本次黑色素瘤各臂的种子全部在 STRING 符号集中，没有额外的 `VEGFA`；因此为 **19,699 + 1 = 19,700** 节点。两次图的主阈值源边记录数均为 **1,858,944**；多出的 ACC 节点没有引入边。

差异基因为 **VEGFA**；共同额外节点是 **MALAT1**。本次黑色素瘤运行没有依赖该 ACC 专有孤立节点，当前结果无需因节点数差异而重跑。不过 VEGFA 在原 ACC 种子权重中的原始值为 `0.26`，占归一化重启权重约 `1.14%`，因此不能把它对原 ACC 数值的影响断言为零；若修改 ACC 图节点或种子定义，须作为单独敏感性分析并保留原结果。

## 执行命令

以下命令从项目根目录 `D:/03-SCI-ACC药物重定位` 运行：

```powershell
python -m pytest tests/test_melanoma_rules.py tests/test_melanoma_audit.py -q -p no:cacheprovider
python -m analysis.melanoma_audit --validate-breast --string-dir . --out results/melanoma_stageB_20261019/validate
python -m analysis.melanoma_audit --project manuscript_revision/plos_one_R1_20261018/stageB_runtime_20261019 --string-dir . --out results/melanoma_stageB_20261019/benchmark --benchmark
python -m analysis.melanoma_audit --project manuscript_revision/plos_one_R1_20261018/stageB_runtime_20261019 --string-dir . --out results/melanoma_stageB_20261019/run
```

要求交付的四个文件在 `run/`：`decisions.json`、`audit_metrics.csv`、`run_manifest.json`、`drug_statistics.csv`。机器可读计时结果在 `benchmark/benchmark.json`，真实乳腺癌比对在 `validate/breast_regression.csv`。其它正式运行输出也全部保存在 `run/`，没有删减。

## 2026-10-19 后续补算：H5 的 124 药次要评估

上表中的 `PENDING_124_ANCHORS` 是**原始阶段 B 运行时**的状态，未回写原始 `run/decisions.json`。此后重新传播并补算 ACC、乳腺癌在 124 药集合中的四种归一化分数；先通过 108 药存档逐药回归核验，再按冻结规则补判 H5。124 药次要评估为 `SUPPORTED` 且必须标注 `ARITHMETICALLY_DETERMINED`，不能解读为独立经验支持。逐药分数、24 行跨疾病对照、16 药明细、补判及哈希见 [`../melanoma_h5_124_supplement_20261019/cross_disease_comparison_124.md`](../melanoma_h5_124_supplement_20261019/cross_disease_comparison_124.md) 和同目录机器可读文件。
