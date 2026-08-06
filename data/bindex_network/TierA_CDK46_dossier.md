# Tier A 旗舰候选深度证据档案：CDK4/6 抑制剂在 ACC 的重定位

**候选**：Ribociclib、Palbociclib、Abemaciclib（+ Trilaciclib）
**结论**：CDK4/6 抑制是本框架**证据最一致、生物学最扎实**的 ACC 重定位方向——计算命中、独立临床前证据、ACC 基因组学、以及我们自己的 MIPE 药敏四条线索交叉吻合，且带有可操作的 **RB1 状态生物标志**。建议列为首选湿实验验证方向。

---

## 1. 基因组学理据：ACC 的 p53/Rb 细胞周期通路高频改变
TCGA 全基因组 ACC 研究（Zheng et al., *Cancer Cell* 2016，91 例）：
- **p53/Rb 细胞周期通路在 44.9% 的 ACC 中改变**（TP53/CDKN2A/RB1/CDK4/CCNE1 合计）。
- 复发性局灶**扩增**：**CDK4**（12q14.1）、**CCNE1**（19q12）、TERT、TERF2。
- 复发性**缺失**：**RB1**（13q14.2）、**CDKN2A**（9p21.2）、ZNRF3。
- 单基因频率：TP53 21%、CTNNB1 16%、CDKN2A 15%、TERT 14%、PRKAR1A 11%。

→ CDK4 扩增 + CDKN2A（p16）缺失 = CDK4/6 活性去抑制，是 CDK4/6i 的**经典敏感性情境**；该通路近半数 ACC 受累，靶点人群充足。

## 2. 框架内部信号（四条线索之一：计算）
| 药 | ADRS v3 rank | C_ACC | 选择性残差 | S_Bneighbor | 外部证据 | 网络靶点（含）|
|---|---|---|---|---|---|---|
| Ribociclib | **#6** | 0.53 | **0.81** | 0.81 | 0.60 | MYC, TRIM24, SF3B1 |
| Palbociclib | #24 | 0.66 | 0.61 | 0.61 | **0.65** | **RB1**, MYC(ρ=.60), NOTCH1, FBXW7 |
| Abemaciclib | #28 | **0.77** | 0.69 | 0.28 | 0.55 | **RB1**, MET, CUX1, FAS |
| Trilaciclib | #5 | 0.72 | (无MIPE) | 0.82 | 0.50 | **RB1**, MYC, ETV6, RUNX1 |

关键点：B-index 药理基因组网络**自发捕捉到 RB1 为 Palbociclib/Abemaciclib/Trilaciclib 的共享靶基因**，与 CDK4/6–RB 轴机制一致——非人为注入，是数据驱动的机制内证。

## 3. ⚠ 细胞抑制型（cytostatic）药敏解读——为何 MIPE 排名不靠前不是坏消息
MIPE 48h 单药 Z-AUC：Abemaciclib avgZ=-0.49（rank 584）、Palbociclib -0.16（rank 862）、Ribociclib 缺失。NCI-60 -logGI50 仅 5.1–6.0（GI50 ~1–8 µM）。
- CDK4/6i 机制是**诱导细胞周期停滞/衰老**（cytostatic），不是快速杀伤（cytotoxic）。**48h 活力/GI50 assay 系统性低估 cytostatic 药**。
- 因此它们在"绝对活性层"排名中等是**机制预期**，不构成反证；这也正是为什么 ADRS 的**网络/情境层 + 选择性残差**能把它们正确上浮，而纯活性层不能。

## 4. 独立临床前证据（四条线索之二：外部实验）
Fiorentini et al.（*Aging* 2018，PMC5764399）在 ACC 细胞系：
- **SW-13（RB 阳性）**：Palbociclib IC50 15.5 µM、Ribociclib 19.1 µM，均降活力并**诱导衰老**。
- **NCI-H295R（RB1 纯合缺失，RB 阴性）**：**仅 Palbociclib 有效**（IC50 14.1 µM，经**凋亡**；且走 RB 非依赖的 GSK3β↑→活性 β-catenin↓→AXIN2 通路）；**Ribociclib 无效**（H295R 抵抗）。
- 结论：**RB1 功能缺失 = CDK4/6i 敏感性的负预测标志**（canonical）；pRB 磷酸化状态可作靶点调控 biomarker。

## 5. ⭐ 内部一致性亮点（四条线索之三 + 之四交叉验证）
我们的 **MIPE ACC 药敏与 RB1 生物标志高度吻合**：
- **NCI-H295R（RB1 缺失）**：Palbociclib H295R Z=+0.05、Abemaciclib H295R Z=-0.03 —— 近乎无效，**与 RB 缺失致 CDK4/6i 抵抗一致**。
- **CU-ACC2**：Abemaciclib Z=**-1.55**（明显活性）、Palbociclib Z=-0.61 —— 提示 CU-ACC2 更可能 RB 完整/CDK 通路依赖。
→ 独立的 MIPE 数据**自洽地重现了 RB1-依赖的响应模式**，强化生物标志假设。

## 6. 生物标志与患者选择（可操作）
- **敏感亚群**：RB1 完整 + CDK4 扩增 / CDKN2A(p16) 缺失 / CCNE1 扩增、pRB 阳性。
- **抵抗亚群**：RB1 缺失/失活（ACC 中 13q14 缺失常见）——除 Palbociclib 可能经 RB 非依赖机制部分克服外，一般预期抵抗。
- 落地：以 IHC pRB + CDKN2A/CDK4 拷贝数 / RB1 缺失作为入组分层。

## 7. 三药差异化与选择
- **Abemaciclib**：连续给药、CNS 穿透、兼抑 CDK9；本框架 C_ACC 最高、CU-ACC2 活性最强 → **首选深挖**。
- **Palbociclib**：唯一在 RB 阴性 H295R 显示活性（RB 非依赖机制）+ 外部证据分最高 → **RB 缺失亚群的独特候选**。
- **Ribociclib**：选择性残差最高、ADRS 最高（#6），但 H295R 抵抗 → 适合 RB 完整亚群。
- 组合覆盖：Abemaciclib（RB+）+ Palbociclib（RB−）可覆盖两个亚群。

## 8. 组合策略与关键药物相互作用
- ⚠ **Mitotane–CYP3A4 相互作用（关键）**：mitotane 是强 CYP3A4 诱导剂，而 CDK4/6i 多为 CYP3A4 底物 → 与标准 ACC 治疗（mitotane）合用会**显著降低 CDK4/6i 暴露**。sunitinib 在 ACC 失败部分即因此。设计时须剂量/PK 监测或考虑非 mitotane 背景。
- 理性联合：CDK4/6i + IGF-1R/mTOR（ACC 的 IGF2 轴；虽 IGF-1R 单药 III 期失败，联合可能协同）；CDK4/6i + MEK；CDK4/6i 维持 EDP 化疗后。

## 9. 湿实验验证方案（建议）
1. **细胞panel**：CU-ACC1、CU-ACC2、NCI-H295R（RB−）、SW-13（RB+）；先做 RB1/CDKN2A/CDK4 基因型与 pRB IHC 分型。
2. **assay 设计**：**避免只用 48h 活力**；采用长程克隆形成/生长曲线 + 衰老标志（SA-β-gal、p21/p16）+ pRB/cell-cycle（流式）——匹配 cytostatic 机制。
3. **药物**：Abemaciclib、Palbociclib、Ribociclib 剂量-反应，按 RB 状态分层比较。
4. **组合**：± IGF-1R/mTOR 抑制剂；测 mitotane 共处理对 CDK4/6i 暴露/疗效的影响。
5. **进阶**：ACC PDO/PDX（含 RB 完整与缺失模型）验证；类固醇生成读出（steroidogenesis）作 ACC 特异药效标志。

## 10. 风险与局限
- ACC 中 RB1 缺失亚群对经典 CDK4/6i 可能抵抗——须生物标志分层，勿全人群给药。
- CDK4/6i 多为 cytostatic，单药或难显著缩瘤，或更适合维持/联合。
- 与 mitotane 的 CYP3A4 相互作用是真实临床障碍。
- 现有 ACC 临床前证据为 2 株经典细胞系；需 CU-ACC1/2 + PDO/PDX 扩展。

## 来源
- TCGA-ACC：Zheng et al., *Cancer Cell* 2016 — https://www.cell.com/cancer-cell/fulltext/S1535-6108(16)30160-X
- CDK4/6i in ACC 临床前：Fiorentini et al., *Aging* 2018 — https://pmc.ncbi.nlm.nih.gov/articles/PMC5764399/
- CDK4/6i 差异综述：JNCI Cancer Spectr 2023 — https://academic.oup.com/jncics/article/7/4/pkad045/7208873
- ACC 精准治疗综述 — https://www.sciencedirect.com/science/article/pii/S1040842825002276
