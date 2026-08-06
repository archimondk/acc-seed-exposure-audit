# ADRS v2 —— 选择性校正版 ACC Drug Repurposing Score

## 相对 v1 的两处修正
1. **ACC 选择性（核心）**：引入 NCI-60 pan-cancer 效力基线（CellMiner DTP，绝对 -logGI50，124/124 NSC 全匹配）。
   - `S_ACCselective(d)` = ACC 效力（MIPE -mean ZAUC）对 NCI-60 一般效力做回归后的**残差百分位**——即"超出广谱细胞毒预期的 ACC 效力"。
   - ACC 效力 vs NCI-60 效力 Pearson r=0.42（有真实共享细胞毒方差可去除）；回归 acc=0.29+0.42·nci。
2. **单靶点置信度**：网络分量（C_ACC、S_B-neighbor）按 `conf=n/(n+2)` 向中性 0.5 收缩，压制低靶点药的网络虚高。S_B-neighbor 改为传播 S_ACCselective（非绝对活性）。

配置：STRING 400 · RWR α=0.4 · C_ACC 收缩 k=3 · B_ACC λ=1 · 权重 {C 0.35, Sel 0.40, Bnb 0.25}。

## QC —— 混杂被有效压制
| 药 | v1 rank | v2 rank | 说明 |
|---|---|---|---|
| **Doxorubicin** | **2** | **51** | 单靶点 BRCA1 虚高被置信度修正；选择性仅 0.41 |
| Vinblastine | 高 | 99 | 广谱最强效但 ACC 平庸（sel 0.20）→ 正确压制 |
| Paclitaxel | — | 34 | 紫杉类下沉 |
| Epirubicin | 18 | 52 | 蒽环下沉 |

## QC —— ACC 机制阳性上浮
- **IGF-1R**：Brigatinib #19（IGF2 是 ACC 头号驱动，强机制命中）。
- **CDK4/6**：Ribociclib #7、Palbociclib #28、Abemaciclib #32（CDKN2A/CDK4/RB1 轴）。
- **MEK**：Cobimetinib #15、Binimetinib #30、Trametinib #40。
- **Proteasome**：Ixazomib #23、Carfilzomib #38。
- **EGFR**：Afatinib #2。**PARP**：Olaparib #13（DDR/BRCA）。
- **Mitotane #21**：唯一获批 ACC 药，选择性 0.87，落在合理中高位（未被埋没）。

## Top 20（selectivity-corrected）
Actinomycin D、Afatinib、Trilaciclib、Mitomycin、Crizotinib、Daunorubicin、Ribociclib、Homoharringtonine、Zanubrutinib、Arsenic trioxide、Romidepsin、Ibrutinib、Olaparib、Irinotecan、Cobimetinib、Cytarabine、Docetaxel、Ceritinib、Brigatinib、Celecoxib。

## 残留说明与下一步
- **Actinomycin D #1 / Daunorubicin #6 / Mitomycin #4**：仍为细胞毒化学型，但具**真实 ACC 选择性残差**（多靶点、命中 CHD4/MYB/BRCA1 等 ACC hub），非纯广谱假象；建议标注"cytotoxic-but-ACC-enriched"，正式候选时人工复核。
- **16 个无 MIPE 活性药**（Trilaciclib #3、Zanubrutinib #9、Arsenic trioxide #10 等）仅靠 C_ACC+S_B-neighbor 排名，**低置信度**，已在 activity_known 列标记。
- 下一步（Layer 5 + 权重调优）：
  1. δ·S_external：PDC/PDO/xenograft/文献（TAK-243/UBA1 等）。
  2. 有监督权重标定：用 MEK/proteasome/CDK4-6/IGF-1R 阳性与 EDP-M/纯细胞毒阴性回归 α/β/γ。
  3. 机制聚类去冗余出最终 Top 10–20。
