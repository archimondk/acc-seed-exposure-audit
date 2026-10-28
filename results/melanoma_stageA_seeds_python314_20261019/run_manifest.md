# 阶段A本地复跑与STRING映射核验

最新状态：用户已完成GitHub发布，2026-10-19核验通过。阶段A tag指向af7487a，六文件一致；以下网络失败及tag待办描述保留为历史过程，已由文末发布核验结论取代。阶段B仍未冻结、未执行。

日期：2026-10-19（Asia/Singapore）。协议：melanoma_R1_v0.3及确认的Amendment A1。

## 结果

- 5/5原始输入SHA-256与收到的Python 3.13运行清单一致。
- 三份冻结文件及A1的SHA-256通过协议检查；复跑代码与收到的脚本逐字节一致。
- 37/37输出CSV与收到的Python 3.13结果逐字节一致。比较范围包含队列、突变、RPPA、种子候选和各臂权重表。
- 原生pytest：31 passed，0 failed，0 errors，0 skipped；详见pytest.xml。
- STRING preferred_name精确匹配：7个分支、13个实际生成的分析臂全部通过，所有种子均映射成功；没有少于2个映射种子的臂。
- 主队列：106例E/K（92例E、14例K）和125例参照；主臂M0为ARID2、CDKN2A、PTEN、TP53，M1增加BRAF。主规则RPPA入选0个特征。

机器可读核验：[verification_report.json](verification_report.json)。映射明细：[string_mapping_report.json](string_mapping_report.json)。主臂供后续使用的已检查种子表：[seeds_frozen_string_checked.csv](primary_metastasis/seeds_frozen_string_checked.csv)。

原始seeds_frozen.csv保持复跑时的PENDING_LOCAL_CHECK字段，便于逐字节比较；映射结果写入各分支seeds_frozen_string_checked.csv。全部种子映射成功，因此无须删除成员或重新归一化权重。映射通过不等于阶段B已冻结。

## 环境

| 项目 | 收到的结果清单 | 本地复跑 |
|---|---|---|
| Python | 3.13.16 | 3.14.3 |
| NumPy | 2.5.3 | 2.4.4 |
| SciPy | 1.18.1 | 1.17.1 |
| 系统 | Linux | Windows |
| pytest | 原环境使用替代脚本，据用户说明31项通过 | 原生pytest 9.0.3，31项通过 |

本地版本与项目.python-version及requirements-lock.txt对应版本一致。输出一致结论限于本次脚本及37份CSV；未将本次复跑扩展为全部研究的有效性证明。

## 来源与执行命令

收到的包：`C:/Users/008bu/Downloads/3/melanoma_stageA_package`。

冻结协议保存于`manuscript_revision/plos_one_R1_20261018/04_melanoma_protocol/stageA_v0.3/`，保留原v0.1草案。Python 3.13参考结果保存于`results/melanoma_stageA_reference_python313_20261019/`。

输入来自`data/melanoma_feasibility_20261018/pancan/`；为适配脚本将case_lists/cases_sequenced.txt与其余四份文件放在同一层的预期，建立`stageA_input_view/`只读使用的硬链接视图，未编辑源数据。

以下命令从`D:/03-SCI-ACC药物重定位`运行：

```powershell
python -I analysis/melanoma_seeds.py --data data/melanoma_feasibility_20261018/stageA_input_view --protocol manuscript_revision/plos_one_R1_20261018/04_melanoma_protocol/stageA_v0.3 --out results/melanoma_stageA_seeds_python314_20261019
python -I -m pytest tests/test_melanoma_seeds.py -q -p no:cacheprovider --junitxml=results/melanoma_stageA_seeds_python314_20261019/pytest.xml
python -I scripts/stageA_string_mapping.py 9606.protein.info.v12.0.txt.gz results/melanoma_stageA_seeds_python314_20261019
python -I scripts/verify_melanoma_stageA_reproduction.py --reference results/melanoma_stageA_reference_python313_20261019 --actual results/melanoma_stageA_seeds_python314_20261019 --protocol manuscript_revision/plos_one_R1_20261018/04_melanoma_protocol/stageA_v0.3
```

STRING文件SHA-256：`144de4b0d98c6a7dfde6ddc2591cf88657f27b989eadff4f501450c3ed1f0f1c`，与原流程reproducibility/input_manifest.json一致。未发生STRING下载。

## 各分支映射种子数

| 分支 | M0 | M1 | 判定 |
|---|---:|---:|---|
| primary_metastasis | 4/4 | 5/5 | 通过 |
| S2_E_only | 3/3 | 4/4 | 通过 |
| S3_neoadjuvant_No | 4/4 | 5/5 | 通过 |
| S7_all_samples | 5/5 | 7/7 | 通过 |
| S8_uniform_weights | 4/4 | 未生成此臂 | 通过 |
| S9_multigene_split | 7/7 | 8/8 | 通过 |
| A1sens_strict_braf | 5/5 | 6/6 | 通过 |

S4、S5、S6及传播相关分支不在本次收到的脚本实现内，未将其标记完成。

## 本地提交与tag

项目根目录不是Git仓库。GitHub远程读取首次连接被重置，HTTP/1.1重试仍无法连接443端口。因此建立了仅含六份指定文件的独立本地快照仓库：

`D:/03-SCI-ACC药物重定位/release/melanoma_stageA_git_20261019`

- 本地分支：protocol-stage-a-snapshot。
- 提交：`1177c0a3124c593d667f560897aa67c455ed32d0`。
- annotated tag：`melanoma-protocol-stageA`。
- tag对象：`8c492aa468713606de4dee5a262e9a7118d00a6e`。
- 创建时间：2026-10-19T07:15:24+08:00，使用实际时间。
- 六份提交文件：FREEZE_STAGE_A.txt、FREEZE_STAGE_A_ADDENDUM.txt、AMENDMENT_A1_braf_classification.md、PROJECT_PLAN.md、protocol_config.json、prefreeze_network_check.csv。
- tag中六个文件的字节与冻结原件一致；工作区干净；bundle验证通过。

这是单独的协议快照提交，没有继承远端主分支历史。仅发布该tag可保留远端主分支，不应强推该本地分支替代main。

**GitHub远端tag未创建或核实，第三方时间戳尚未取得。** 本地提交时间晚于收到的种子结果，不能用该tag证明种子生成前已经公开冻结；原冻结和确认时间保留为收到的记录。

可交付文件：

- 六份冻结文件ZIP：`release/melanoma_stageA_frozen_files_20261019.zip`。
- 含提交和tag的Git bundle：`release/melanoma_stageA_20261019.bundle`。

以下为用户网页上传前准备的备用命令，现已停用：用户已将六文件提交到远端main，应在远端提交af7487a上创建tag，不再推送独立本地快照tag，以免产生同名tag冲突。

```powershell
git -C "D:/03-SCI-ACC药物重定位/release/melanoma_stageA_git_20261019" push origin refs/tags/melanoma-protocol-stageA
```

若远端已有同名tag，先比较内容，不覆盖或使用force。发布后核验远端tag指向和六文件哈希，才能把第三方时间戳登记为完成。

## 用户网页提交后的远端核验（2026-10-19）

- 用户已自行网页上传。GitHub API确认提交 `af7487a25ad86a8d5c4e49c880115c4181754c6b`，消息为 `Add melanoma Stage A frozen protocol and A1`；核验时main指向此提交。
- 读取该固定提交下六份文件的Git blob并解码原始字节，逐份计算SHA-256；6/6与本地冻结原件完全一致。
- 核验时远端尚无 `melanoma-protocol-stageA` tag，下一步由用户在网页为此远端提交创建阶段A tag。
- 远端提交与本地独立快照提交不同属正常现象；文件内容相同，不应将本地tag强推覆盖将要创建的远端tag。
- 公开提交现已可核验，但发生在收到的阶段A结果之后，不构成结果生成前公开预注册的证明。
- 提交链接：https://github.com/archimondk/acc-seed-exposure-audit/commit/af7487a25ad86a8d5c4e49c880115c4181754c6b

## 远端tag及Release发布核验完成（2026-10-19）

- Release：https://github.com/archimondk/acc-seed-exposure-audit/releases/tag/melanoma-protocol-stageA
- GitHub API确认draft=false、prerelease=false；发布时间2026-10-19T01:56:55Z，即Asia/Singapore 2026-10-19 09:56:55。
- 远端轻量tag `melanoma-protocol-stageA` 直接指向提交 `af7487a25ad86a8d5c4e49c880115c4181754c6b`，与此前逐份SHA-256核验通过的提交相同。
- 再次以tag为ref取得目录，六份文件Git blob哈希均与本地原件一致。
- latest-release接口仍返回v1.0.3，本次阶段归档未替代正式最新版。
- 本次用户请求的STRING映射、Python 3.14复跑与pytest核验、六文件提交及阶段A tag发布均已完成。阶段B不在本次执行范围，仍未冻结、未执行。
- 保留本地独立快照tag作为历史备份，不推送或强制覆盖现有远端同名tag。
