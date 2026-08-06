# 最终 ADRS（全 STRING 主分析）—— 定稿数字

配置：**全人类 STRING v12 互作组（≥400，1,858,944 边）** RWR（disease-only 种子，α=0.4）→ r_ACC → C_ACC（收缩 k=3）→ B_ACC（λ=1）→ + ACC-relative potency residual（NCI-60 基线）+ S_B-neighbor（置信度加权，n/(n+2)）+ S_external。权重 {C 0.30, residual 0.35, Bnb 0.20, ext 0.15}。脚本 `run_ADRS_final_fullSTRING.py`，产物 `ADRS_final_fullSTRING_ranked.csv`。

## 与诱导网络版一致性（关键）
- 最终 ADRS 药物排名 **Spearman(full vs induced) = 0.972**；top-20 Jaccard = 0.74。
- → 换成更严格的全互作组后，结论稳定，定性判断不变。

## 定稿基准（非循环）
- 纯计算分对临床标签：**AUC = 0.40，95%CI [0.15, 0.68]，置换 P = 0.635**（仍不显著、CI 跨 0.5）→ 计算评分不能可靠预测 ACC 临床疗效（结论不变）。
- S_B-neighbor LOMCO Spearman = **−0.111**（类级交叉验证下仍无独立预测力）。

## 最终 Top 20（ADRS v3）
Actinomycin D、Crizotinib、**Mitotane(#3)**、Trilaciclib、Mitomycin、Irinotecan、Daunorubicin、Homoharringtonine、Arsenic trioxide、Ibrutinib、Romidepsin、Zanubrutinib、Cytarabine、Docetaxel、**Cobimetinib(#15)**、**Palbociclib(#16)**、Brigatinib、Axitinib、**Ribociclib(#19)**、Belinostat。

## 关键药物：最终排名（括号为诱导版）
| 药 | 最终 | 诱导版 | 说明 |
|---|---|---|---|
| Mitotane | **3** | 3 | ACC 标准，阳性锚点，稳定 |
| Cobimetinib (MEK) | 15 | 17 | Tier B |
| Palbociclib (CDK4/6) | 16 | 24 | Tier A |
| Ribociclib (CDK4/6) | 19 | 6 | Tier A（全网络下略降）|
| Ixazomib (proteasome) | 21 | 23 | Tier B |
| Olaparib (PARP) | 22 | 11 | Tier B |
| Abemaciclib (CDK4/6) | 25 | 28 | Tier A |
| Doxorubicin | 43 | 42 | 细胞毒混杂被压制，稳定 |
| **Afatinib (EGFR)** | **46** | 7 | 全网络 + 外部证据下大幅下沉，符合 EGFR 临床失败 |
| Cabozantinib | 62 | 71 | 临床有活性但计算低（框架局限的诚实体现）|
| Vinblastine | 105 | 100 | 广谱最强效但 ACC 平庸，正确压制 |

## 定性结论（全部不变，且更稳）
1. 细胞毒混杂被 pan-cancer-adjusted residual 压制（Doxorubicin/Vinblastine 低位）。
2. CDK4/6 仍为 Tier A 证据一致假设（Palbociclib 16 / Ribociclib 19 / Abemaciclib 25），Mitotane #3 锚定。
3. EGFR/Afatinib 经外部证据降至 #46（临床已失败）。
4. 计算评分不能可靠预测临床疗效（AUC 0.40，P=0.635）——方法学警示成立。
5. 全网络 vs 诱导网络 0.972 一致 → 网络构建方式不影响结论。
