# 故障排查

## 流水线在 00_validate_inputs.R 就停了

这是**设计如此**，不是 bug。四类硬门禁：

| 报错关键词 | 含义 | 处理 |
| --- | --- | --- |
| `数据类型不合规` | 数据集是 RNA-seq / ChIP-seq，不是芯片 | 换数据集。`gdstype` 必须含 `array` |
| `物种不合规` | 非人源 | 换数据集 |
| `样本量不合规` | ≥ 10 例 | 换数据集，或按 GSM 筛选子集另存为新的 GSE |
| `分组失败` / `分组歧义` | `group_field` / `group_values` 与数据不匹配 | 用 `node scripts/find_dataset.mjs samples GSEXXXXX` 看样本真实字段值 |

排查分组问题：

```bash
node scripts/find_dataset.mjs check GSE92252   # 会列出识别到的候选组
node scripts/find_dataset.mjs samples GSE92252 # 每个样本的完整 characteristics
```

## GitHub Actions 在 20 分钟被杀掉

`timeout-minutes: 20` 是 spec 的硬约束。超时几乎总是**包在源码编译**：

1. 看 job 日志里 `setup-r-dependencies` 阶段的耗时。
2. 确认 `setup-r@v2` 的 `use-public-rspm: true` 生效 —— 它让 pak 从 Posit Package
   Manager 拉 Linux 二进制包而不是源码。
3. 若 Bioconductor 包仍在编译，把 R 版本与 Bioconductor 版本对齐
   （R 4.3 ↔ Bioc 3.18），版本错配会强制走源码。
4. 最后的兜底：改用 `container: bioconductor/bioconductor_docker:RELEASE_3_18`，
   所有包已预装，job 通常 3–5 分钟跑完。代价是偏离 spec 里显式的 apt + setup-r 步骤。

## KEGG 没有结果

`enrichKEGG` 走 KEGG 在线 REST API，会因限流、网络或授权返回空/报错。
查 `results/enrichment_status.json` 的 `kegg.status` 与 `kegg.reason`。

**报告中不得写"无 KEGG 通路富集"**，只能写"本次未获得 KEGG 结果"。

## PPI 没有出图

查 `results/ppi_status.json`：

- `status: skipped` → 显著 DEG 少于 5 个，无法建网。这是 n=9 的正常表现。
- `status: fallback` → STRINGdb 不可用，已回退为**共表达网络**。
  **这不是物理蛋白互作，不得当作 PPI 证据引用。**
- `status: failed` → STRING 失败且共表达网络中无 `|r| >= 0.9` 的边。

STRINGdb 首次运行要下载约 100 MB 网络文件，超时或网络受限时会失败。

## 热图基因数不是 50

`results/enrichment_status.json` 的 `heatmap_mode` 会写明实际用了什么：

- `top 50 significant DEG by adj.P` → 正常
- `all N significant DEG (< 50 requested)` → 显著基因不足 50，已降级
- `WARNING: no significant DEG, fell back to top 20 by raw P` → 无显著基因

后者在 n=9 设计下是正常的，**不要**通过放宽阈值来"修复"它。

## GO 富集结果为空

先看 `results/deg_summary.json` 的 `n_significant`。若只有几个显著基因，
富集必然为空 —— 这是样本量的限制，不是代码问题。

若显著基因很多但富集仍为空，检查 `enrichment.universe`：`detected` 会把背景限制在
实测基因内，可能过小。

## `impute` 包不可用

spec 要求 KNN 填补。若存在缺失值而 `impute` 未安装，脚本会直接 `stop()` 而不是
静默跳过 —— 静默跳过会让下游统计建立在缺失值上。装上它，或确认数据本身无缺失。

## 本地没有 R

`find_dataset.mjs` 和 `tools/check_r_syntax.mjs` 都不需要 R。完整流水线需要 R 4.3+ 与
Bioconductor，本地跑不了时直接推送到 GitHub Actions。
