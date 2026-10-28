# 本地运行：黑色素瘤分支S4–S6（2026-10-19）

S4、S5的种子和S6的描述已在云端由原始datahub文件生成（`results/melanoma_S4_S6/`）。本地需要做三件事：复跑核对、STRING映射、S4/S5的确定性传播。全部命令在仓库根目录执行（与`melanoma_audit`相同的布局：代码放在`analysis/`，`--project`指向返修目录）。

## 0. 放置文件

| 文件 | 放到 |
|---|---|
| `IMPLEMENTATION_NOTE_S4_S6.md` | `<返修目录>/04_melanoma_protocol/` |
| `melanoma_seeds_s4_s6.py`、`melanoma_branches_s4_s5.py` | `analysis/` |
| `test_melanoma_s4_s6.py` | `tests/` |
| `results/melanoma_S4_S6/`（云端输出） | `<返修目录>/results/` |

实施说明是运行前写下的，SHA-256为`50b70e9334ac1cbc0091e61e5f61518e06d4f82be99e30ab43c33ab15249cf71`；请提交并打tag（例如`melanoma-S4-S6-note`），留第三方时间戳。

## 1. 测试

```
pip install statsmodels        # 用于HC3与statsmodels的逐项比对；未安装时该项自动跳过
python -m pytest tests/test_melanoma_s4_s6.py tests/test_melanoma_seeds.py -q
```

预期：S4–S6测试11项通过（未装statsmodels时10项通过、1项跳过），原阶段A测试31项通过。

## 2. 复跑种子并比对

```
python -m analysis.melanoma_seeds_s4_s6 --data <skcm原始文件目录> --protocol <返修目录>/04_melanoma_protocol --out <返修目录>/results/melanoma_S4_S6_local
```

- 两个`seeds_frozen.csv`应与云端逐字节一致：S4 `f6bb914244d83a7bc13cecd19cff19c7dd90d4051d012a907c1e64348530d85b`，S5 `902a19705906d734d0371aedddee5bd4c6cab779134647ff7812219cf08a01a0`。
- 其余表格因NumPy/SciPy版本不同，第12位有效数字可能有差异，请按数值比对。
- **请核对S4识别出的两个样本与任务2登记的两个样本ID是否相同**：TCGA-D9-A6E9-06（V600E）、TCGA-ER-A2NF-06（参照）。不同则停止，不继续第3、4步。

## 3. STRING映射

```
python -I scripts/stageA_string_mapping.py <路径>/9606.protein.info.v12.0.txt.gz <返修目录>/results/melanoma_S4_S6
```

PDCD4是S5新增的种子，必须映射成功；任何臂少于2个映射种子则该臂停止。

## 4. S4/S5确定性传播（无零模型）

```
python -m analysis.melanoma_branches_s4_s5 --project <返修目录> --branches <返修目录>/results/melanoma_S4_S6 --string-dir <STRING目录> --out <返修目录>/results/melanoma_S4_S5_propagation
```

脚本先用同一引擎复算阶段B主队列M0/M1的124药百分位，与`drug_statistics.csv`的最大差超过5×10⁻¹¹就停止；随后按`seed_branches.csv`的口径输出两个分支、两个臂、四种归一化的结果。预计几分钟。

## 5. 请发回

- `results/melanoma_S4_S6/string_mapping_report.json`
- `results/melanoma_S4_S5_propagation/seed_branches_S4_S5.csv`和`run_manifest.json`
- 第2步的比对结论（种子是否逐字节一致；S4两个样本ID是否与任务2一致）
- pytest输出
