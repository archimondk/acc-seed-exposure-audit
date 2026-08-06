# 独立验证：ACC_CellMinerCDB（无湿实验）

数据源：ACC_CellMinerCDB（Lin et al., *Cancer Res Commun* 2024）。用其中**独立于 MIPE 的 NCI-CCR 药敏筛选**（47 药，按 NSC 索引，测于 CU-ACC1/CU-ACC2/H295R/SW13 + 6 个手术样本，-log10 IC50）与 **NCATS**（=我们已有的 MIPE，2665 化合物，同 3 株 ACC 系）交叉比对。

## 验证 1：活性层跨实验室可复现性（PASS，中-良）
同样 3 株 ACC 系上，NCI-CCR vs NCATS 两个独立筛选的 -log10 IC50 相关：
| ACC 系 | n(共测药) | Pearson | Spearman |
|---|---|---|---|
| CU-ACC1 | 21 | 0.61 | **0.83** |
| CU-ACC2 | 20 | 0.41 | 0.58 |
| NCI-H295R | 16 | 0.44 | 0.61 |

结论：两个独立机构、不同 assay 协议下，ACC 药敏排序中-良度一致（Spearman 0.58–0.83）。**直接回应"活性层单一来源（MIPE）"的审稿质疑**——活性信号可复现。（n 小、协议不同，故非完美，如实说明。）

## 验证 2：细胞毒性混杂是"跨数据集普适现象"（强化 v2 选择性校正的必要性）
NCI-CCR 里绝对效力最高的仍是广谱细胞毒药：**Docetaxel（-logIC50 8.9–9.7）、Doxorubicin（5.9–6.8）**，远超靶向药。
→ 广谱细胞毒主导 ACC 体外效力这一现象，在 **MIPE、GDSC(SW13)、NCI-CCR 三个独立数据集中一致重现**。这为 ADRS v2 的"用 NCI-60 基线做选择性校正"提供了强有力的、多数据集的合理性支撑。

## 验证 3：CDK4/6（Tier A）—— 确认"cytostatic 需专门 assay"，非否定
NCI-CCR 中 Palbociclib（-logIC50 ~4.6–4.75，panel 内 17–22 百分位）、Ribociclib（~3.9–4.3，2–7 百分位）**在短程 IC50 下同样显得弱**，且 H295R(RB 缺失) 与其他系差异不明显。
- 这与 MIPE、GDSC 一致：**cytostatic 药在 IC50/AUC 短程 assay 里被系统性低估**，短程 IC50 也无法分辨 RB 依赖的衰老表型。
- 即独立数据**既不能证实也不能证伪** CDK4/6——它证实的是**assay 局限**，正好印证 Tier A 档案的建议：CDK4/6 必须用长程克隆形成 + 衰老标志 + pRB，而非 IC50。

## 附带资源（已下载，待用于生物标志验证）
- `ACCcmdb_NCICCR_rnaseq.txt`（6.4MB，RNA-seq）、`ACCcmdb_NCICCR_cn.txt`（2MB，拷贝数）：含 6 个**真实 ACC 手术样本** + 细胞系。
- 可用于确认 CDKN2A/CDK4/RB1 状态、SOAT1/MGMT/SLFN11/ABCB1 等 seed/biomarker 基因在真实 ACC 肿瘤中的表达/拷贝数 → 把 C_ACC 疾病情境层"落到真实肿瘤"。

## 对论文的意义（诚实版）
- **补上了独立验证**：活性层跨实验室可复现（验证1）；选择性校正的必要性被三数据集共识支撑（验证2）。这两条是实打实、可写进 Results 的独立验证，且不需湿实验。
- **CDK4/6 结论保持审慎**：现有一切公开 ACC 药敏都是短程 IC50/AUC，对 cytostatic 药无判别力；故 CDK4/6 仍是"强机制假设 + 待专门 assay 验证"，不夸大为"已验证有效"。
- **仍存的局限**：ACC 可用模型极少（3 细胞系 + 6 手术样本 + PDX），无大样本药敏基准；ADRS 全排名无法被独立药敏逐一验证，只能验证到层级/机制与可复现性。

产物：`ACCcmdb_NCICCR_focus_drugs.csv`（47 药 × 4 系 -logIC50）。
