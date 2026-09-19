# 故障排查

> **路径约定：** 下文所有 `results/<GSE>/X` 与 `data/<GSE>/X` 里的 `<GSE>` 指
> 当前 config 的 `dataset_id`。产物按数据集分目录，排查时先确认自己看的是哪个数据集 ——
> 拿 GSE64790 的状态文件去解释 GSE42568 的报错是很容易犯的错。

## 流水线在 00_validate_inputs.R 就停了

这是**设计如此**，不是 bug。门禁按 `design_mode` 分两套：

| 报错关键词 | 含义 | 处理 |
| --- | --- | --- |
| `数据类型不合规` | 数据集是 RNA-seq / ChIP-seq，不是芯片 | 换数据集。`gdstype` 必须含 `array` |
| `物种不合规` | 非人源 | 换数据集 |
| `样本量不合规` | `small_sample` 要求 < 10，`cohort` 要求 ≥ 15 | 换数据集，或改 `design_mode`（**不要为了通过而调大上限**） |
| `每组样本数不合规` | `small_sample` 要求 ≥ 3，`cohort` 要求 ≥ 10 | 换数据集，或按 GSM 筛选子集另存为新的 GSE |
| `分组失败` / `分组歧义` | `group_field` / `group_values` 与数据不匹配 | 用 `node scripts/find_dataset.mjs samples GSEXXXXX` 看样本真实字段值 |
| `pairs[[k]] 必须...` | `paired: true` 但 `pairs` 里有样本不属于本数据集 / 重复 / 不成对 | 逐对核对 GSM，每对必须一例 numerator 一例 denominator |
| 样本未配对 | `paired: true` 但有样本没出现在 `pairs` 里 | 补齐 `pairs`，或改 `paired: false` |
| `没有指定配置文件` | 直接 `Rscript scripts/xxx.R` 没带 `--config` | `parse_args()` 故意没有默认配置，见报错里列出的可用配置 |

排查分组问题：

```bash
node scripts/find_dataset.mjs check GSE64790   # 会列出识别到的候选组
node scripts/find_dataset.mjs samples GSE64790 # 每个样本的完整 characteristics
node tools/check_sample_structure.mjs GSE64790 # 分组是否与批次效应混杂
```

## GitHub Actions 被杀掉

`timeout-minutes: 30`。**实测（run 35438543496）**：Install R packages 66s
（增量，`restore-keys` 命中旧缓存）、Run analysis GSE64790 244s / GSE42568 400s，
暖缓存整轮 6m10s / 8m49s。全冷缓存下装包约 720s，整轮约 20 分钟。

超时几乎总是**包在源码编译**：

1. 看 job 日志里 `setup-r-dependencies` 阶段的耗时。
2. 确认 `setup-r@v2` 的 `use-public-rspm: true` 生效 —— 它让 pak 从 Posit Package
   Manager 拉 Linux 二进制包而不是源码。
3. 若 Bioconductor 包仍在编译，把 R 版本与 Bioconductor 版本对齐
   （R 4.3 ↔ Bioc 3.18），版本错配会强制走源码。
4. 最后的兜底：改用 `container: bioconductor/bioconductor_docker:RELEASE_3_18`，
   所有包已预装，job 通常 3–5 分钟跑完。代价是偏离 spec 里显式的 apt + setup-r 步骤。

> 上限之所以写 30 而不是 20：被掐死的 job 存不下缓存，下一轮又是冷缓存，
> 会变成"每次都超时"的死循环。

## WGCNA 报 `unused arguments (weights.x = NULL, weights.y = NULL, cosine = FALSE)`

`blockwiseModules` 内部算 KME 时用 `do.call(corFnc, ...)`，而 `corFnc` 来自包内
常量 `.corFnc = c("cor", "bicor", "cor")` —— 是**字符串**，按名字查找。
不 attach 任何包时它解析到 `stats::cor`，后者没有那三个参数。

**传 `corFnc = WGCNA::cor` 没用**：`blockwiseModules` 形参表里没有 `corFnc`，
参数掉进 `...`，而 KME 那段用包内常量、不看 `...`。

处理：调用期间 `library("WGCNA")`，`on.exit` 立刻 `detach`。见 AGENTS.md 规则 23。

> 注意别写成 `library(WGCNA, character.only = TRUE)` —— 那会去**求值** `WGCNA`
> 这个符号，报 `object 'WGCNA' not found`，看着像"包没装"。

## 步骤 06 / 07 没产物，但 CI 是绿的

它们 `required = FALSE`，失败不让 job 变红。查
`results/<GSE>/wgcna_status.json` 与 `lasso_status.json` 的 `status` 字段：

- `ok` → 真跑了
- `not_applicable` / `not_configured` / `too_few_events` / `too_few_samples`
  / `package_missing` / `empty_signature` → 跳过，`reason` 里写明原因
- **文件不存在，或 `status` 字段缺失** → 步骤崩了。
  `check_acceptance()` 的 `settled()` 会记 FAIL，翻 `state.json` 找 `error`。

## LASSO 的 `cindex_train` 明显低于 `cindex_cv`

不该出现。两者互补（例如 0.121 对 0.793）说明**风险分方向反了** ——
`survival::concordance()` 的公式接口在 `Surv(time, event) ~ risk` 下把预测子
当成生存方向，返回 1 - Harrell C。

本仓库用自己的 `harrell_c()`（`07_lasso.R`），定义写在注释里。
**这类错误不会引起怀疑**：0.121 完全可能是一个真的很差但合理的模型。

## 外部验证 `validation: failed`

看 `validation_error` 的**前缀**，它标明是哪一步：

- `getGEO:` / `exprs:` / `pData:` → 下载或解析
- `fetch_platform_annotation:` → 平台注释抓取
- `map_features_to_symbols:` → 探针映射
- `collapse_to_symbol:` → 基因折叠
- `clinical_table:` → SOFT 临床字段解析
- `探针->基因映射长度 N != 探针数 M` → `map_features_to_symbols()`
  返回的是 **list**，要用 `$symbols` / `$mode`，不能当向量（踩过）

## KEGG 没有结果

`enrichKEGG` 走 KEGG 在线 REST API，会因限流、网络或授权返回空/报错。
查 `results/<GSE>/enrichment_status.json` 的 `kegg.status` 与 `kegg.reason`。

**报告中不得写"无 KEGG 通路富集"**，只能写"本次未获得 KEGG 结果"。

## PPI 没有出图

查 `results/<GSE>/ppi_status.json`：

- `status: skipped` → 显著 DEG 少于 5 个，无法建网。这是 n=9 的正常表现。
- `status: fallback` → STRINGdb 不可用，已回退为**共表达网络**。
  **这不是物理蛋白互作，不得当作 PPI 证据引用。**
- `status: failed` → STRING 失败且共表达网络中无 `|r| >= 0.9` 的边。

STRINGdb 首次运行要下载约 100 MB 网络文件，超时或网络受限时会失败。

## 热图基因数不是 50

`results/<GSE>/enrichment_status.json` 的 `heatmap_mode` 会写明实际用了什么：

- `top 50 significant DEG by adj.P` → 正常
- `all N significant DEG (< 50 requested)` → 显著基因不足 50，已降级
- `WARNING: no significant DEG, fell back to top 20 by raw P` → 无显著基因

后者在 n=9 设计下是正常的，**不要**通过放宽阈值来"修复"它。

## GO 富集结果为空

先看 `results/<GSE>/deg_summary.json` 的 `n_significant`。若只有几个显著基因，
富集必然为空 —— 这是样本量的限制，不是代码问题。

若显著基因很多但富集仍为空，检查 `enrichment.universe`：`detected` 会把背景限制在
实测基因内，可能过小。

## `impute` 包不可用

spec 要求 KNN 填补。若存在缺失值而 `impute` 未安装，脚本会直接 `stop()` 而不是
静默跳过 —— 静默跳过会让下游统计建立在缺失值上。装上它，或确认数据本身无缺失。

## 本地没有 R

`find_dataset.mjs` 和 `tools/check_r_syntax.mjs` 都不需要 R。完整流水线需要 R 4.3+ 与
Bioconductor，本地跑不了时直接推送到 GitHub Actions。
