---
name: geo-brca-microarray-skill
description: Run an end-to-end GEO gene-expression microarray data-mining pipeline for human breast cancer with fewer than 10 samples - download and clean, QC, PCA, sample correlation, limma differential expression, clustered heatmap, GO/KEGG enrichment, and STRING PPI - and ship it as a GitHub Actions workflow that uploads results as an artifact. Use when the user asks for GEO data mining, a GEO/GSE microarray analysis, a breast cancer expression analysis, a small-sample limma DEG workflow, GO/KEGG enrichment on GEO data, a GitHub Actions bioinformatics pipeline, or asks to reproduce or extend GSE92252.
license: MIT
compatibility: R 4.3+ with Bioconductor (GEOquery, limma, clusterProfiler, org.Hs.eg.db, enrichplot, STRINGdb, pheatmap, impute, igraph). Node 18+ for the dataset pre-flight tool. Needs outbound HTTPS to NCBI GEO, STRINGdb and KEGG. Designed to fit GitHub Actions ubuntu-latest in under 20 minutes.
metadata:
  version: "1.0"
  skill-author: built for DSH from a user-supplied agent_spec
---

# /geo-brca-microarray-skill

GEO 乳腺癌小样本芯片数据挖掘流水线。默认数据集 **GSE92252**（人源、NimbleGen 单色芯片、
9 例 = 6 肿瘤 vs 3 正常乳腺组织）。

完整实验设计见 [`EXPERIMENTAL_DESIGN.md`](EXPERIMENTAL_DESIGN.md) —— **改任何东西之前先读它**，
尤其是 §2.2 关于 n=9 的统计功效边界。

## 触发场景

- 「用 GEO 数据做乳腺癌差异表达分析」「GSE 芯片数据挖掘」
- 「小样本（<10）芯片怎么做 limma + GO/KEGG」
- 「搭一个 GitHub Actions 生信流水线，跑完上传 artifact」
- 「复现/扩展 GSE92252」
- 任何要求「清洗 + QC + PCA + 相关性 + 差异基因 + 热图 + 富集 + PPI」全流程的任务

## 工作流

### 1. 先验证数据集，不要相信标题

```bash
node scripts/find_dataset.mjs check GSE92252
node scripts/find_dataset.mjs search --disease "breast cancer" --max-samples 10
```

这一步不是可选的。GEO 的标题与二手描述经常与真实元数据不符 —— 本项目最初提出的
GSE197894 被描述为"表达谱芯片"，实际 `gdstype` 是 RNA-seq，样本 20 例，
照它做流水线会在第一步就终止。**类型必须以 `gdstype` / `GPL` 为准。**

### 2. 改配置，不改代码

所有可调项都在 [`assets/config.yml`](assets/config.yml)：数据集、分组字段、
分组模式、对比方向、阈值、富集背景。切换数据集只需要改这个文件。

### 3. 跑

```bash
# 本地
Rscript scripts/main_analysis.R --config assets/config.yml

# GitHub Actions
gh workflow run geo_analysis.yml
```

`main_analysis.R` 按 00→05 顺序执行，每步独立 `tryCatch`：必需步骤失败即中止，
可选步骤（富集、PPI）失败则记录原因后继续，最后按 spec 的 acceptance criteria
逐项校验产物并写 `results/state.json`。

### 4. 读结果

产物清单见 `EXPERIMENTAL_DESIGN.md` §4。`results/state.json` 是唯一的执行真相来源。

## Gotchas

- **`GES` 不是真实的 GEO 前缀。** 正确的是 `GSE`（GEO Series）。用户说 GES 时按 GSE 处理。
- **样本量 < 10 且要两组对比，几乎排除了所有组织样本数据集。** 人源乳腺癌芯片里
  n<10 的绝大多数是细胞系加药实验，没有"肿瘤 vs 正常"两组。GSE92252 是少数例外
  （6 肿瘤 + 3 正常）。若用户给的数据集是细胞系，分组字段应改为 treatment/control，
  而不是硬套 tumor/normal。
- **20 分钟不是宽裕的预算。** `clusterProfiler` + `STRINGdb` + `enrichplot` +
  `org.Hs.eg.db` 从源码编译要 25–40 分钟。workflow 用
  `r-lib/actions/setup-r-dependencies@v2` 走 Posit Package Manager 二进制包来压时间。
  如果换了包或 R 版本导致退化为源码安装，`timeout-minutes: 20` 会直接杀掉 job。
- **`enrichKEGG` 依赖 KEGG 在线 REST API**，会因限流或授权返回空。按 spec 这不算失败，
  但**报告中不得写"无 KEGG 通路富集"**，只能写"本次未获得 KEGG 结果"，并引用
  `results/enrichment_status.json` 里的原因。
- **STRINGdb 首次运行要下载约 100 MB 网络文件。** 失败时本流水线回退到基于表达谱的
  **共表达网络**（Spearman |r| ≥ 0.9），`ppi_status.json` 会标记 `coexpression_fallback`。
  共表达不是物理互作，**不得当作 PPI 证据引用**。
- **显著 DEG 可能少于 50 个。** 热图会自动降级（全部显著基因 → top 20 by P），
  `enrichment_status.json` 的 `heatmap_mode` 会记录实际用了什么。这是 n=9 的正常表现，
  不是 bug。
- **背景集选择会改变富集结果。** `enrichment.universe: genome`（默认，spec 要求）
  与 `detected`（仅实测基因）结果不同，报告中必须写明用了哪个。
- **不要为 tumor/normal 差异编造机制解释。** 肿瘤 vs 全组织正常的差异主要来自
  组织成分（上皮/脂肪/基质/免疫），不是肿瘤特异驱动事件。见设计文档 §2.2。
- **GSE92252 不做配对分析。** 正常组 N39/N40/N54 中只有 T40/T54 有对应肿瘤样本，
  3 对里只有 2 对成立。`paired: false` 是有依据的，不要"优化"成配对。

## 文件结构

```text
SKILL.md                     本文件
EXPERIMENTAL_DESIGN.md       实验设计（分组依据、统计功效边界、局限）
AGENTS.md                    给其他 agent 的仓库约定
README.md                    人类可读的安装与使用说明
assets/config.yml            全部可调参数
scripts/
├── main_analysis.R          编排器（00→05 + 验收）
├── 00_validate_inputs.R     硬门禁：物种/类型/样本量/分组
├── 01_download_clean.R      下载、探针映射、KNN 填补、quantile 标准化
├── 02_qc_pca_correlation.R  QC 箱线图/密度、PCA、样本相关性
├── 03_deg.R                 limma 差异表达 + 火山图
├── 04_heatmap_enrichment.R  聚类热图 + GO/KEGG
├── 05_ppi.R                 STRING PPI（含共表达回退）
├── find_dataset.mjs         GEO 数据集预检（Node，不需要 R）
└── lib/common.R             配置、日志、状态、SOFT 抓取
references/troubleshooting.md 运行期故障排查
.github/workflows/geo_analysis.yml
```
