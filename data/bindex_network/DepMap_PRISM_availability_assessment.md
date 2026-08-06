# DepMap/PRISM 中 ACC 数据可得性评估（v2 无湿实验验证用）

## 结论先行
**可行且值得做，但 ACC 细胞系极少（约 2 株）**——罕见癌在所有 pan-cancer 数据库里代表性都很低。因此 DepMap 的最大价值不是"大规模 ACC 药敏基准"，而是**基因层功能验证（CRISPR 依赖性）+ 跨数据集活性一致性**，这两条恰好能补上审稿人最可能质疑的"缺独立验证"。

## 已确认的 ACC 细胞系覆盖
| 细胞系 | DepMap | 药敏 | 备注 |
|---|---|---|---|
| **SW-13** | ✓ ACH-001401 | GDSC / PharmacoDB / (PRISM) | NCIt/Orphanet 归为 ACC；SMARCA4 Q164* + TP53 H193Y；RB 阳性 |
| **NCI-H295R** | ✓（CCLE 经典系）| GDSC / PRISM / PharmacoDB | 经典 ACC 系；**RB1 纯合缺失**（CDK4/6i 抵抗生物标志）|
| CU-ACC1 / CU-ACC2 | ✗ | ✗ | 2018 学术系，未进 DepMap/CCLE/GDSC |

- DepMap 数据类型（两株都有）：**CRISPR 基因依赖性（Achilles/Chronos）**、表达、拷贝数、突变；药敏经 GDSC/PRISM。
- PRISM secondary：1448 药 × 499 系；primary：4518 药 × ~578 系——含这两株 ACC 系，与我们 124 FDA 药有部分重叠。

## 三条可行的无湿实验验证（按价值排序）
**1. CRISPR 依赖性验证基因层（最强、最新）**
用 H295R + SW-13 的 Chronos gene effect：检验高 r_ACC / ACC-seed 基因是否在 ACC 系中**优先必需**（比 pan-cancer 更负）。这是对 C_ACC / disease-context 层的**功能性、正交**验证，独立于药敏，最有说服力。
- 附加：H295R（RB1 缺失）对 **CDK4/6 依赖性应低**——正向检验 CDK4/6-RB 生物标志。

**2. 跨数据集活性一致性（验证活性层）**
我们的 MIPE H295R 药敏 vs GDSC/PRISM H295R 药敏做 concordance——证明活性层在独立筛选间可复现，回应"MIPE 单一来源"质疑。

**3. 基因组学佐证生物标志**
用 DepMap 突变/CN 确认 H295R RB1 缺失、SW-13 RB 完整 + SMARCA4 缺失，坐实 CDK4/6 的 RB 状态分层论述（第一手数据，不只引文献）。

## 局限（诚实写入 limitations）
- **仅约 2 株 ACC 系**（MIPE 有 3 株，两边仅 H295R 重叠）→ 无法做稳健的"跨 ACC 系药敏相关"，n 太小。
- 因此定位为**功能/机制正交验证 + 一致性检验**，不是大样本药敏基准。
- SW-13 是否严格 ACC 学界有争议（小细胞型/来源存疑），需注明。

## 下载清单（用户侧，一次性）
- DepMap CRISPR (Chronos) gene effect：`CRISPRGeneEffect.csv`（较大 ~400MB）——或用按基因子集的 API 导出。
- DepMap Model/样本注释：`Model.csv`（定位 ACC 系 ACH ID，很小）。
- 药敏：GDSC2 `GDSC2_fitted_dose_response.xlsx`（较小）或 PRISM secondary（`secondary-screen-dose-response-curve-parameters.csv`）。
- 突变/CN：`OmicsSomaticMutations.csv` / `OmicsCNGene.csv`（大，可选，仅取 2 系）。
来源：https://depmap.org/portal/data_page/ ；GDSC：https://www.cancerrxgene.org/downloads

## 建议
先下**小文件**（Model.csv 定位 ACC 系 + GDSC2 药敏 + CRISPR 仅这 2 系的子集），即可完成验证 1–3，不必下全量大矩阵。完成后这篇的"独立验证"短板基本补齐，可进入写稿。
