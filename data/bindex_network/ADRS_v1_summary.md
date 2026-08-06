# ADRS v1 —— 第一版 ACC Drug Repurposing Score 结果与诊断

## 配置
STRING score≥400 · RWR α=0.4 · C_ACC 收缩 k=3 · B_ACC λ=1 · 权重 {C_ACC 0.35, S_activity 0.40, S_B-neighbor 0.25}（可调）。
所有分量先转 124 药内百分位再加权；分量缺失时按可用权重重归一化。

## 层次
- **Layer 3 S_activity**：108/124 有 MIPE Z-AUC；activity = -mean(Z-AUC) 在全 MIPE 分布的百分位。16 个新型靶向药/老烷化剂无 MIPE 数据。
- **Layer 4 S_B-neighbor（LODO）**：候选药对 108 个活性已知参考药的 B_ACC 加权活性均值，留一法防泄漏。为 16 个无活性药提供网络插补活性。
- **Layer 6 ADRS**：三分量加权。

## Top 结果（前 12）
Actinomycin D、Doxorubicin、Homoharringtonine、Mitomycin、Romidepsin、Afatinib、Abemaciclib、Daunorubicin、Docetaxel、Trilaciclib、Ixazomib、Ixabepilone。

## 机制阳性对照（QC）
- **CDK4/6 抑制剂**：Abemaciclib rank 7、Ribociclib rank 20 —— **真阳性**，契合 ACC 的 CDKN2A 缺失 / CDK4 扩增 / CCNE1-RB1 轴。
- **Proteasome**：Ixazomib rank 11、Carfilzomib rank 24。
- **EGFR**：Afatinib rank 6。
- **MEK**：Cobimetinib 23、Binimetinib 33、Trametinib 39、Selumetinib 54（中高，未进最顶）。

## ⚠ 关键诊断：细胞毒性混杂（cytotoxicity confound）
Top 段被广谱强细胞毒药主导——蒽环类（Doxorubicin #2、Daunorubicin、Idarubicin、Epirubicin）、Actinomycin D、Mitomycin、紫杉类。**Doxorubicin 排到第 2**，而它正是临床失败的 EDP-M 方案成员（P0.5-A 已证 48h 单药体外活性 ≠ 临床疗效）。

原因：
1. S_activity（绝对体外效力）天然奖励"什么都杀"的细胞毒药，不是 ACC 特异性。
2. Doxorubicin 单靶点（BRCA1，r_ACC 高），C_ACC 收缩后仍偏高（0.98），叠加 S_B-neighbor 0.98 被网络分量抬高。

即：**ADRS v1 目前区分不了"体外普遍强效"与"ACC 特异获益"**。这与 EDP-M 反面教训一致，是预期中会暴露、也必须在 v2 修正的问题。

## v2 修正方向（建议）
1. **ACC 选择性而非绝对活性**：用 pan-cancer 基线校正 activity（ACC Z-AUC 相对其他癌系的 selectivity），压制广谱细胞毒。需再拉一层 pan-cancer 药敏（NCI60/MIPE 其他组织）。
2. **单靶点 C_ACC 降权**：对 n_targets=1 的药加更强收缩或置信度惩罚（Doxorubicin/Olaparib 类）。
3. **机制去冗余**：对蒽环/DNA 损伤类做机制聚类去重，避免 Top 被同一类刷屏。
4. **权重按阳性调优**：用已知 ACC-active（MEK、proteasome、CDK4/6、UBA1 方向）与阴性（EDP-M）做有监督权重标定，而非人工固定。
5. δ·S_external：接 PDC/PDO/xenograft/文献（TAK-243 等）。

## 结论
框架已端到端打通并通过机械性验证（λ=0 精确复现原始 B-index；CDK4/6 等机制阳性被召回）。ADRS v1 作为**诊断性基线**是成功的：它既召回了有生物学意义的 ACC 靶点（CDK4/6、EGFR、proteasome、MEK），也清晰暴露了细胞毒性混杂这一必须解决的核心问题。下一步优先做"ACC 选择性校正"，把 v1 的绝对活性升级为疾病特异活性。
