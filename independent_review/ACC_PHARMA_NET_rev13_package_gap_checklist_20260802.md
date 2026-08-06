# M5 投稿包缺口清单 — ACC-PHARMA-NET rev13

**生成日期**: 2026-08-02
**检查范围**: `ACC_Supplementary_Materials_list.md` 引用的全部 117 个路径 + rev13 正文引用 + `figures/`、`figure_data/`、`analysis/`、`tests/`、`results/`、`leakage_audit/` 实盘资产
**方法**: 批量存在性检查 + 全项目递归文件名检索 + 通配符展开

---

## 检查结果总览

| 类别 | 引用数 | 存在 | 缺失 | 备注 |
|---|---:|---:|---:|---|
| 补充表 S1–S29 数据文件 | 87 | 45 直查存在,45 经递归找到 | 0(2 个为解析误报) | 索引多数条目只写文件名、未写目录路径 |
| 代码文件 Code S1–S18 | 32 | 32 | 0 | 全部脚本与测试存在 |
| 通配符路径 | 2 | 2(5 个 lambda 文件 + 6 臂 seeds) | 0 | 全部展开存在 |
| 主图资产(figures/revision) | — | 22 个文件 | 0 | 见 B 节命名问题 |
| 补充图 S1–S6 | — | S3/S4/S5 有资产 | **S1/S2 无资产;S6 未完成** | 见 A 节 |

**结论:数据与代码本体全部就位;投稿包缺口集中在(1)补充图资产、(2)投稿元数据、(3)文件命名与路径组织。**

---

## A. 投稿前必须完成(严重缺口)

### A1. 补充图 S6 尚未制作最终资产
- 正文 3.9 引用 "Supplementary Figure S6A–C"(跨平台一致性)+ S6D(生物标志物)
- 图例 rev13 自述:"final asset renumbering deferred until text freeze"
- 现状:仅有旧版 `figures/Fig6_validation.png`(早期六图时代,非三面板规范版);数据源 `figure_data/Fig6abc_crossplatform_concordance.csv`、`Fig6d_biomarker_expression.csv` 存在
- 动作:按 S6A–C/S6D 三面板重新渲染,提供 PNG+PDF+SVG 三格式(MDPI 要求),删除或归档旧 Fig6

### A2. 补充图 S1/S2 无资产文件
- 索引收录 "Figure S1"(四分量权重网格,legacy)与 "Figure S2"(LOMCO 评估,legacy),标注 "excluded from primary robustness claims",但**未列 source 文件,项目内无对应资产**
- 动作:二选一——(a) 重新渲染并提供资产;(b) 从投稿包中移除 S1/S2 并同步修订正文/索引引用(正文未直接引用 S1/S2,选项 b 更省事,但需在图例中保留说明或删除图例条目)

### A3. Data Availability:仓库 URL 与 DOI 未插入
- 正文自述 "The resolving repository URL and DOI will be inserted before submission"
- 动作:创建匿名评审可访问的仓库(如 Zenodo/figshare 预印 DOI 或 GitHub private-invite),将 URL/DOI 写入 Data Availability;Code S17/S18 注明 "release implementation must reproduce the frozen protocol before DOI deposition"——发布版代码需先通过复现门禁

### A4. 全体作者批准语句与 CRediT 未完成
- 稿件内注释 "ALL-AUTHOR CONFIRMATION REQUIRED BEFORE SUBMISSION. After all five authors have approved the roles and final manuscript, append: 'All authors have read and agreed to the published version of the manuscript.'"
- 动作:五位作者批准 CRediT 分工与最终稿后,补上该语句;同时删除稿件中所有 HTML 工作注释

---

## B. 文件命名与路径组织(建议投稿前整改)

### B1. 主图资产文件名与新编号不一致(图例自述)
| 新编号(正文) | 当前文件名 | 需改为 |
|---|---|---|
| Figure 3 | `Fig5_weight_stability_CDK46_primary108.*` | `Fig3_*.pdf/png/svg` |
| Figure 5 | `Fig7_seed_target_leakage_audit.*` | `Fig5_*.pdf/png/svg` |
| Figure S5 | `Fig3_component_correlation_primary108.*` | `FigS5_*.pdf/png/svg` |
| Figure S6 | `Fig6_validation.png`(旧版) | `FigS6_*.pdf/png/svg`(重做) |

- 正文/图例/文件名三者必须一致,否则 MDPI 排版与审稿人核对都会出错
- Figure 1 资产为 SVG + word.png(无 PDF),建议补 PDF 版本
- `Fig7_seed_target_leakage_audit_manifest.json` 随重命名同步更新

### B2. 补充材料索引中 45 个条目缺目录路径
- 例:`ADRS_evidence_informed_108.csv`(实际在 `results/primary_analysis/`)、`verdict.json`(实际在 `leakage_audit/`)、`arm_metrics.json`(实际在 `leakage_audit/arms/*/`)
- 索引的 Source 列当前只写文件名,打包脚本按索引取值会找不到文件
- 动作:把 Source 列全部改为项目根相对路径(与 `results/...`、`leakage_audit/...` 一致),并跑一次打包验证

### B3. 补充索引中 `.pdf`/`.svg` 裸扩展名写法
- 索引中 "Figure S3 的 source" 写作 `` `figures/revision/FigS3_dirichlet_weight_sensitivity.png`; `.pdf`; `.svg` ``——反引号包裹的裸扩展名会被路径工具误解析
- 动作:改为完整路径三连写(与 B1 重命名同步完成)

---

## C. 建议性整改(不阻塞投稿但提升完整性)

1. **README.md 过时**:仍声明 "Current manuscript: rev5"(实际 rev13);更新为 rev13 并同步 "Current manuscript" 小节
2. **旧图资产归档**:`figures/` 根目录的 Fig1–Fig6 旧版(早期分析)建议移入 `figures/archive_legacy/`,避免与 revision 资产混淆;`figure_data/` 根目录旧 CSV 同理
3. **锁文件清理**:根目录存在 `.~lock.ACC_PHARMA_NET_Pharmaceutics_rev4.pdf#`(LibreOffice 残留锁),投稿打包前删除
4. **临时目录清理**:`tmp/`、`__pycache__/`、`.pytest_cache/`、`lu54bxcvj.tmp` 不进投稿包;若仓库公开,建议加入 .gitignore 并确认
5. **机构信息确认**(AUTHOR_METADATA_REQUIRED_rev13.md 已列):机构 2 英文译名(复合医院名)需两家机构确认;通讯作者电话是否加 `+86` 前缀由作者确认
6. **共同第一作者**:目前未登记共同一作;若最终稿需要,须在投稿前由全部作者确认

---

## D. 复现性相关(投稿时可附带)

- 12 个冻结输入、31 个科学输出、12 个图文件的复现门禁均通过(`repro_outputs/major_revision_v6_20260728/`,898.9s,49/49 测试)——已核验该 manifest 存在
- 建议在仓库 README 中给出一键复现命令与运行时间,方便审稿人执行

---

## 优先级排序(若时间有限)

1. **A1(S6 重做)与 A2(S1/S2 取舍)** — 直接影响正文引用完整性,审稿人必查
2. **A3(仓库/DOI)** — MDPI 投稿系统必填,且影响可复现性评价
3. **B1(图重命名)** — 排版与核对必踩
4. **B2(索引补路径)** — 打包验证必踩
5. **A4/C 其余项** — 投稿流程内完成即可
