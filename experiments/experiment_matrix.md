# 实验矩阵表：ACC-PHARMA-NET 方法强化

| 实验ID | 对应假设 | 数据集 | baseline | 指标 | 随机性 | 状态 | 完成判定 |
|---|---|---|---|---|---|---|---|
| EXP-01 | 构成审计（原 H1 已退休） | C1 complete-case 108 | raw MIPE、residual、C_ACC、n_assoc、direct seed overlap、S_external | Spearman、Top-20 Jaccard | 无 | 已完成 | 7 个排名均覆盖 108 药物；仅作描述，不再据此推断增量价值 |
| ABL-01 | 构成审计/H3 | 同一 108 药物 | `ADRS_comp` 各单一分量 | CDK4/6 平均秩、精确单侧 P、跨排名 BH q | 无 | 已完成 | 所有分量使用同一 3 个 CDK4/6 药物和 108 药物零集合 |
| EXP-02 | H2 | STRING v12 ≥400 + 399 关联基因 | degree、strength、PageRank | Spearman rho、95% bootstrap CI | bootstrap seed 20260727 | 已完成 | 三种中心性的 Spearman rho 和 95% CI 均输出且 n=399 |
| ABL-02 | H2 | STRING v12 ≥400 + 399 关联基因 | PageRank 去除 degree 解释部分 | 部分 Spearman rho | 无 | 已完成 | n=399 且输出有限的部分 Spearman rho（`-1≤rho≤1`） |
| ROB-01 | H2/H3 | STRING v12 ≥400 + ACC 种子 | 10,000 个 degree-matched random seed sets | 逐药 Z、经验 P、BH q；CDK4/6 组经验 P | RNG seed 20260727；B=10000 | 已完成 | 10,000 次全部完成；同 degree bin；经验 P≥1/10001；BH 分辨率门通过 |
| ROB-02 | H2/H3 | 同上 | 四种传播/归一化分支 | 排名相关、Top-20 Jaccard、中心性相关、逐药 BH q、CDK4/6 P/q | 同 ROB-01 | 已完成 | 四个预定义分支均输出；结果差异按模型依赖性解释 |
| QC-01 | H2 | 冻结 `rACC_399_fullSTRING.csv` | 从原始 STRING/种子重算 | 最大绝对差、Spearman | 无 | 已完成 | 最大绝对差≤1.1e-6 且 Spearman≥0.999999 |

## 派生数据规格

本阶段没有外部派生数据集。ROB-01 的随机种子集由代码在运行时生成，
保存到结果目录并由固定 RNG seed 和输入哈希完全确定；不会回流修改主输入。

## 算力与成本

| 实验ID | 设备 | 重复次数 | 预计时数 | 单价 | 小计 |
|---|---|---:|---:|---:|---:|
| EXP-01/ABL-01 | 本地 CPU | 确定性 | <0.1 h | 已有设备 | 0 |
| EXP-02/QC-01 | 本地 CPU | 确定性 | <0.2 h | 已有设备 | 0 |
| ROB-01/ROB-02 | 本地 CPU | 4 × 10,000 null draws | 以最终 manifest 墙钟时间为准 | 已有设备 | 0 |
| 合计 | | | 以最终 manifest 为准 | | 0 |

统计说明：这里的 10,000 次重复是随机化零模型而非模型训练种子。经验
P 的最小分辨率为 `1/10001`，并报告 Monte-Carlo 标准误；不以 3–5 个
训练种子的口径替代随机化检验所需的重复数。
