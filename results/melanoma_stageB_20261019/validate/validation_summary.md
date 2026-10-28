# 真实 STRING 乳腺癌回归核验

命令：`python -m analysis.melanoma_audit --validate-breast --string-dir . --out results/melanoma_stageB_20261019/validate`

程序标准输出：

```json
{
  "breast_regression": "PASS",
  "max_abs_diff": 4.953260024365136e-13,
  "rows": 432
}
```

冻结阈值为最大绝对差不超过 `5e-11`。逐行对照见 `breast_regression.csv`；核验使用项目根目录的真实 STRING v12 info 和 links 文件。
