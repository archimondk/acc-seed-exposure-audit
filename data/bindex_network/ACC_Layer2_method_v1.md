# ACC-PHARMA-NET Layer 2 — PPI Propagation & ACC Contextualization (v1)

## 1. 输入数据
- **Original B-index network**（Table S1，Pharmaceutics 2025;17:1421，DOI 10.3390/pharmaceutics17111421）：124 药 × 399 基因 × 1304 边，QC 完全匹配。
- **ACC seed set**：50 基因 + ACC_weight v1.0（`ACC_P0.5C_gene_weights_v1.csv`）。
- **PPI 网络**：STRING v12 human，取 431 基因（399 宇宙 ∪ 50 seed）诱导子网，`required_score≥700`（高置信度），2437 条边。
  - 47/50 seed、366/399 宇宙基因在图内；ABCB1、SLFN11、VEGFA 无≥700 边（保留其 restart 质量但无扩散）。

## 2. 传播算法：Random Walk with Restart（个性化 PageRank）
- 邻接 W：STRING combined score 作为无向加权边。
- 转移矩阵 T：按列归一化；孤立节点列为 0（质量回到 restart）。
- Restart 向量 r0：seed 的 ACC_weight 归一化到和为 1。
- 迭代：`p = (1-α)·T·p + α·r0`，α=0.40，47 步收敛（L1<1e-12）。
- **r_ACC(g)**：取 399 宇宙基因的稳态访问概率，min-max 归一化到 [0,1]；同时提供 rank-percentile 版本。

产物：`rACC_399_v1.csv`（raw / minmax / percentile / is_seed / in_graph）。

## 3. Drug ACC Context Score
`C_ACC(d) = mean[ r_ACC(g) : g ∈ targets(d) ]`
产物：`CACC_drug_context_124_v1.csv`。

## 4. ACC-weighted B-index
`w_ACC(g) = 1 + λ·r_ACC(g)`
`B_ACC(A,B) = 0.5 · Σ_shared w_ACC · (1/Σ_A w_ACC + 1/Σ_B w_ACC)`
λ 网格 {0, 0.5, 1, 2, 5}，产物 `BACC_pairs_lambda*_v1.csv`。

## 5. 关键 QC —— λ=0 退化验证
λ=0 时 w≡1，B_ACC 应等于原始 B-index。与论文 Table S3 的 1337 对 B_index 逐一比较：
- 匹配对数：**1337/1337**
- **最大绝对误差 9.85e-16**（机器精度），diff>1e-6 的对数=0。

→ B-index 计算引擎实现**精确正确**；λ>0 的 ACC 加权是对已验证基线的可控扰动。

## 6. 结果概览
- r_ACC 顶部为宇宙内 seed（BRCA1=1.00、ATM 0.91、CDK4 0.91、TERT 0.88、MSH2/MSH6、RB1、MEN1、BRCA2…），符合 ACC 生物学（DDR/MMR、cell-cycle、TERT、MEN1）。
- 传播把宇宙外核心驱动（TP53、CTNNB1、IGF2、PRKAR1A、ZNRF3）的 ACC 相关性经 PPI 扩散到宇宙内邻居，解决了 18/50 硬交集过稀的问题。

## 7. 已知局限（v1 caveat，进入 ADRS 前需处理）
1. **单靶点药物 C_ACC 不稳定**：124 药中 17 个只有 1 个基因靶点、30 个 ≤2 个。C_ACC 作为均值对这些药噪声大（如 Doxorubicin 单靶点→C_ACC=1.0）。建议：加入 target-count 收缩项，或对 n_targets 设下限 / 贝叶斯 shrinkage。
2. **score≥700 稀疏**：33/399 宇宙基因在高置信图中孤立，r_ACC≈0。可用 required_score=400 densify 做敏感性分析（URL 已备 `STRING_network_URL_score400.txt`）。
3. **α、归一化方式**未调优：需对 α∈{0.2,0.3,0.4,0.5} 与 min-max vs percentile 做敏感性。
4. seed weight v1.0 仍为专家先验，未接 TCGA/GTEx 实证（P0.5-C 待升级项）。

## 7b. 已处理：C_ACC 收缩 + α/阈值敏感性（v2）
脚本 `run_layer2_sensitivity.py`，产物 `CACC_drug_context_124_v2_shrunk.csv`、`layer2_sensitivity_meta.json`。

**(1) C_ACC 单靶点收缩（empirical-Bayes）**
`C_ACC_shrunk(d) = (n_d·mean_d + k·μ0)/(n_d + k)`，μ0=全体 drug-target 实例的 r_ACC 均值=0.0953。
- k 网格 {0,1,2,3,5,10}：k 越大越把低靶点药拉回背景。
- **默认 k=3**：与 raw 的 Spearman=0.91（保序），top20 里单靶点药从 3 个降到 1 个。
- 效果：Doxorubicin（单靶点，raw C_ACC=1.0 的假象）被收缩到 0.321；收缩后 Top10 变为合理的多靶点抗代谢/抗肿瘤药（Pralatrexate、Ixazomib、Pemetrexed、Ixabepilone、Gemcitabine、Cytarabine、Tioguanine、Fluorouracil）。
- **ADRS 一律采用 shrunk k=3 版本**；raw 仅留档。

**(2) α（restart）敏感性，score≥700**
- r_ACC 基因级 Spearman：所有 α∈{0.2,0.3,0.4,0.5} 两两 ≥0.98（相邻 α ≥0.998）。
- C_ACC 药物级 Spearman：≥0.93（相邻 α ≥0.99）；top20 Jaccard(0.3/0.5 vs 0.4)=0.82。
- 结论：**α 几乎不影响结论，固定 α=0.4**。

**(3) STRING 阈值敏感性（α=0.4，700→800→900 上滤）**
- 边数 2437→1595→1010；有信号基因 364→345→277。
- r_ACC Spearman：700vs800=0.92、700vs900=0.84；C_ACC：0.97/0.94；top20 Jaccard 0.82/0.74。
- 结论：**阈值是主要杠杆**（比 α 敏感），越严越丢基因覆盖。→ 建议补测 densify 方向（score≥400，需一次下载）作为覆盖度/稳健性上界；`STRING_network_URL_score400.txt` 已备。

**(4) Densify score≥400（已补测）**
- 边数 2437→8142；有信号基因 364→**391/399**（补回 27 个孤立基因，如 LCP1、MN1、NFIB、CIC、EXT1/2、PRDM16 等）。
- r_ACC 基因级 Spearman(700vs400)=0.915；C_ACC 药物级=0.969；top20 Jaccard=0.90；Top10 药物几乎不变。
- 结论：densify 不改变药物级结论，但覆盖度接近完整。**采用 score≥400 为 ADRS 主网络，score≥700 作高置信度稳健性对照。**
- 产物：`rACC_399_v2_700vs400.csv`、`CACC_124_v2_700vs400.csv`。

**最终默认配置：STRING score≥400 + RWR α=0.4 + C_ACC 收缩 k=3。**

## 8. 下一步（Layer 4 → Layer 6）
- S_B-neighbor(d)：候选药与 ACC 高置信参考药集的 B_ACC 邻域相似性（需 leave-one-drug-out 防泄漏）。
- 融合 S_ACC-activity（MIPE Z-AUC）。
- 第一版 ADRS 实算 → Top 候选。
