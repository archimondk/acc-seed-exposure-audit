# ADRS v3 —— 六层框架收尾：外部证据整合 + 分层候选

## 完整框架（Layer 1–6 全部实算）
1. Pan-cancer B-index 网络（124药×399基因×1304边，精确复现）
2. ACC 疾病情境化（STRING PPI + RWR → r_ACC → C_ACC / B_ACC）
3. ACC 实验药敏（MIPE Z-AUC，108/124）
4. B-index 邻域（S_B-neighbor，LODO）
5. **外部 ACC 证据（本轮新增，Layer 5）**：文献 curation 的 S_external（19 药有证据，含临床/临床前/失败）
6. ADRS 整合排名

配置：STRING 400 · RWR α=0.4 · C_ACC 收缩 k=3 · B_ACC λ=1 · 单靶点置信度 n/(n+2) · 选择性残差（NCI-60 基线）· 权重 {C 0.30, Sel 0.35, Bnb 0.20, Ext 0.15}。

## ⭐ 核心科学发现：计算相关性 ≠ 临床疗效
用**独立外部标签**（10 个有 ACC 证据的阳性药 vs 4 个 ACC 临床失败的阴性药）验证**纯计算分**（C_ACC+Sel+Bnb，不含外部）：
- **AUC = 0.35（< 0.5，弱反相关）**。阴性药中位排名（30）竟高于阳性药（48）。
- 具体：模型最爱的靶向药 **Afatinib（EGFR）计算 #2**，恰是 ACC 临床**失败**的机制（gefitinib II 期 0% 缓解）；而临床有活性的 **Cabozantinib 计算仅 #90**、Gemcitabine #71。
- 原因：ACC 临床史上，靶向单药大多失败（IGF-1R linsitinib III 期失败、EGFR 失败、sunitinib 令人失望），而广谱化疗（EDP-M/gemcitabine）仅温和有效——计算模型偏好"机制上 ACC 相关的靶向药"，方向恰与这段临床史相左。

> 这不是 n=14 小样本要去"调权重抹平"的噪声，而是**罕见癌纯药理基因组重定位评分的根本局限**：网络/药理基因组相似性刻画的是"机制相关性"，与"已证实临床疗效"是两个轴。这正是必须叠加外部证据（Layer 5）与湿实验的理由，也是本框架方法学讨论的核心卖点。

外部整合后 AUC 0.35→0.48（wD=0.15 温和修正：Afatinib #2→#7，Mitotane #21→#3）。**刻意不过度加权外部**以免 ADRS 退化成 S_external 本身（n=14 会过拟合）。

## 分层最终候选（Top 30 → 四层）
**Tier A｜计算命中 + 独立 ACC 证据支持（最高优先）**
- **CDK4/6：Ribociclib #6、Palbociclib #24、Abemaciclib #28** —— 计算命中 + ACC 细胞系（SW-13/NCI-H295R）临床前证据 + ACC 的 CDKN2A缺失/CDK4扩增/RB 轴生物学三重吻合。**旗舰候选。**
- **Mitotane #3** —— 唯一获批 ACC 药，作为阳性锚点自然上浮，佐证框架合理。

**Tier B｜计算命中，暂无外部证据（新重定位假设，待湿实验）**
- MEK：Cobimetinib #17、Binimetinib #29
- Proteasome/UPS：Ixazomib #23（呼应 P0.5-A 的 MLN-7243/UBA1、proteasome 富集）
- HDAC：Romidepsin #14、Belinostat #22
- PARP/DDR：Olaparib #11（ACC 有 DDR/BRCA 情境）
- ALK/ROS1/MET：Crizotinib #4、Ceritinib #21、Brigatinib #19、Tepotinib #30
- BTK：Ibrutinib #16、Zanubrutinib #15

**Tier C｜计算命中但临床已失败（降级）**
- EGFR/HER：**Afatinib #7** —— EGFR 单抑在 ACC 已证失败，除非联合策略（EGFR+IGF-1R 有临床前协同）否则不建议。

**Cytotoxic chemotype｜ACC 富集但广谱（非新靶向机会）**
- Actinomycin D #1、Mitomycin #2、Daunorubicin #8、Homoharringtonine #9、Irinotecan #10、Cytarabine #13、Docetaxel #18、Fluorouracil #27 —— 有真实 ACC 选择性残差，但属细胞毒化学型；正式候选时人工复核。

## 局限
- 外部标签 n=14 且本身混杂（阳性多为温和化疗、阴性多为失败靶向药）；AUC 为噪声大的指示性结果，不宜过度解读。
- S_external 为人工 curation（含来源），非系统性；proteasome/HDAC/PARP/BTK 暂无 ACC 直接文献，列 Tier B 而非否定。
- 16 药无 MIPE 活性，靠网络插补，低置信度已标记。

## 下一步（P4）
1. **湿实验优先级**：CDK4/6（Tier A）作首选验证；Tier B 的 proteasome/HDAC/MEK 作机制新假设。
2. 扩充 S_external 为系统性证据（ClinicalTrials.gov + PDC/PDO/xenograft 数据库）。
3. seed weight v2 接 TCGA/GTEx 实证（P0.5-C 待升级项）。
4. 组合策略建模（EGFR+IGF-1R 协同等），突破单药靶向在 ACC 的历史失败。
