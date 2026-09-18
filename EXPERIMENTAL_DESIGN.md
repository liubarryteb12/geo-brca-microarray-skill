# GEO 乳腺癌小样本芯片数据挖掘 — 实验设计

> 数据源：**GSE92252** ｜ 平台：**GPL16025**（Agilent-072363 SurePrint G3 Human GE v3 8×60K）
> 物种：*Homo sapiens* ｜ 类型：Expression profiling by array ｜ 样本量：**9**
> 运行环境：GitHub Actions `ubuntu-latest` ｜ 语言：R / Bioconductor

---

## 1. 为什么是 GSE92252

### 1.1 合规性核验（来自 GEO 元数据，非推测）

| 约束 | 要求 | GSE92252 实测 | 结论 |
| --- | --- | --- | --- |
| 物种 | Homo sapiens | `Homo sapiens` (taxid 9606) | ✅ |
| 数据类型 | 基因芯片 | `Expression profiling by array` / GPL16025 Agilent 双色芯片 | ✅ |
| 疾病 | 乳腺癌 | AR+/ER−/PR− 乳腺癌组织 vs 正常乳腺组织 | ✅ |
| 样本量 | < 10 | **9**（6 肿瘤 + 3 正常） | ✅ |
| 分组可比 | 需两组 | tumor 6 vs normal 3，同一平台同一批次 | ✅ |
| 本地算力 | 轻量 | 9×60K 矩阵，峰值内存 < 1 GB，全程 < 10 min | ✅ |

### 1.2 被否决的候选（这是本设计最关键的一步）

最初提出的四个数据集**全部不合规**，若直接使用，流水线会在第一步 `validate_inputs` 就终止：

| 数据集 | 物种 | 数据类型 | 样本数 | 否决原因 |
| --- | --- | --- | --- | --- |
| **GSE197894** | 人 ✅ | **RNA-seq** ❌ | **20** ❌ | 新辅助化疗前后配对 RNA-seq，非芯片；样本超限；分组不是肿瘤 vs 正常 |
| GSE143334 | **小鼠** ❌ | ChIP-seq ❌ | 6 | 骨骼肌昼夜节律，与乳腺癌无关 |
| GSE161922 | **小鼠** ❌ | RNA-seq ❌ | 6 | B16 黑色素瘤模型 |
| GSE206338 | **小鼠** ❌ | 表达谱 ❌ | — | 小鼠 NASH 肝纤维化模型 |

> 教训：GSE197894 的检索摘要被描述为"表达谱芯片"，但其 `gdstype` 实际是
> `Expression profiling by high throughput sequencing`。**分组与类型必须以 GEO 的
> `gdstype` / `GPL` 字段为准，不能采信标题或二手描述。**

### 1.3 样本构成

| GSM | 标题 | 分组 | 特征 |
| --- | --- | --- | --- |
| GSM2424491 | BreastCancerTissue-HER2negative-B13 | tumor | HER2− |
| GSM2424492 | BreastCancerTissue-HER2negative-B29 | tumor | HER2− |
| GSM2424493 | BreastCancerTissue-HER2negative-B32 | tumor | HER2− |
| GSM2424494 | BreastCancerTissue-HER2positive-T32 | tumor | HER2+ |
| GSM2424495 | BreastCancerTissue-HER2positive-T40 | tumor | HER2+ |
| GSM2424496 | BreastCancerTissue-HER2positive-T54 | tumor | HER2+ |
| GSM2424497 | NormalBreastTissue-Normal-N39 | normal | 正常乳腺 |
| GSM2424498 | NormalBreastTissue-Normal-N40 | normal | 正常乳腺 |
| GSM2424499 | NormalBreastTissue-Normal-N54 | normal | 正常乳腺 |

分组依据 `characteristics_ch1` 的 `tissue:` 字段：肿瘤样本含 `breast tumor`，正常样本含
`Normal breast tissue`。全部 9 例均为女性。

---

## 2. 设计决策与理由

### 2.1 主对比：tumor (6) vs normal (3)，非配对

`contrast = ["tumor", "normal"]`，limma 中 tumor 为分子、normal 为分母，
即 **log2FC > 0 表示肿瘤中上调**。

**为什么不做配对分析**：提交者的实验记录写明正常组为 N39/N40/N54，肿瘤组为
T39/T40/T54，但实际入库的肿瘤样本是 **T32**/T40/T54。仅 T40、T54 能与正常样本按患者
编号对应，N39 无配对肿瘤样本。3 对中只有 2 对成立，配对设计不成立，因此
`paired: false`。

### 2.2 统计功效：这是本设计必须写在最前面的限制

6 vs 3 是一个**极小样本设计**，必须明确它意味着什么：

- limma 的 `eBayes` 经验贝叶斯收缩依赖足够多的基因-样本组合来稳定方差估计。
  n=9 时 moderation 仍然有效（比普通 t 检验强），但**自由度极低**。
- 在 `adj.P < 0.05` 且 `|log2FC| > 1` 下，本设计**只能检出效应量很大的基因**。
  这对"肿瘤 vs 正常乳腺组织"是可行的 —— 两者在组织构成上差异巨大，效应量本就很大。
- **但这也正是最大的解释陷阱**：肿瘤组织与正常乳腺组织的差异，绝大部分来自
  **组织成分差异**（上皮比例、脂肪、基质、浸润免疫细胞），而非肿瘤特异性驱动事件。
  因此 DEG 列表会强烈富集于细胞外基质、免疫应答、脂肪代谢等通路。

> **判读边界**：本流水线产出的是**探索性假设**，不是肿瘤发生机制结论。
> 任何"某基因驱动乳腺癌"的表述都超出了本设计能支持的范围。
> 若要区分"肿瘤特异"与"组织成分"，需要 LCM 显微切割数据或去卷积（如 CIBERSORTx）。

### 2.3 为什么不做 HER2+ vs HER2− 的次级对比

tumor 组内部还有 3 HER2+ vs 3 HER2− 的结构，看起来可以做次级对比。**本设计不纳入**：
每组 n=3 时 limma 无法给出可信的 adj.P，且多重检验校正后几乎不可能有基因通过。
该对比仅在 `results/deg_table.csv` 中保留 `her2` 注释列供人工查看，不产出独立结论。

### 2.4 芯片特有的技术决策

| 问题 | 决策 | 理由 |
| --- | --- | --- |
| 双色 vs 单色 | 按 GPL16025 的实际结构处理，取 log2 ratio | Agilent 双色芯片，GEOquery 返回的已是 log2 ratio |
| 探针 → 基因 | 多探针取**表达方差最大**者 | 比取均值更能保留真实信号，且避免稀释 |
| 缺失值 | KNN 填补（`impute::impute.knn`, k=10） | 芯片常见；记录填补比例 |
| 标准化 | `limma::normalizeBetweenArrays(method="quantile")` | 跨样本可比；QC 保留 before/after 对照 |
| 过滤 | 去除全零/低表达基因（行中位数低于 25 分位） | 降低多重检验负担 |
| 富集背景 | 默认**全基因组**（OrgDb），可切换为实测基因集 | 见 §3.6 |

---

## 3. 分析流程

八个阶段，每阶段产出明确文件；任一步失败写入 `state.json` 并继续执行不依赖该步的后续步骤。

### 3.1 数据获取与校验（`00_validate_inputs.R`）

- 读取 `assets/config.yml`，校验 `dataset_id` / `group_field` / `group_values` / `contrast` 齐全
- 从 GEO SOFT 接口（base R `url()`，不依赖 Bioconductor）拉取 series 与 sample 元数据
- **硬门禁**：物种必须为 `Homo sapiens`；类型必须为 array；样本数必须 < 10；两组样本数均 ≥ 3
- 按 `group_field` 匹配 `group_values` 生成分组；一个样本命中多个组 → 报错退出
- 输出：`data/group.csv`、`data/geo_metadata.json`

### 3.2 下载与清洗（`01_download_clean.R`）

- `GEOquery::getGEO()` 下载表达矩阵与平台注释
- 探针 → 基因 symbol；无 symbol 的探针丢弃
- 过滤全零基因；KNN 填补缺失
- quantile 标准化
- **一致性校验**：表达矩阵列名必须与 `group.csv` 的样本行完全对应，否则报错
- 输出：`data/expr_raw.rds`、`data/expr_clean.rds`、`data/expr_clean.csv`

### 3.3 质控（`02_qc_pca_correlation.R`）

| 图 | 内容 | 判据 |
| --- | --- | --- |
| `boxplot_before_after.pdf` | 标准化前后各样本表达分布 | 标准化后中位数应齐平 |
| `density_plot.pdf` | 标准化前后密度曲线 | 曲线应重合 |
| `pca_plot.pdf` | PCA 散点（PC1/PC2），按分组着色 | 记录 PC1/PC2 方差解释率 |
| `correlation_heatmap.pdf` | 样本间 Pearson 相关热图 | — |
| `correlation_matrix.csv` | Pearson + Spearman 相关矩阵 | **cor < 0.8 的样本标记为离群** |

### 3.4 差异表达（`03_deg.R`）

- 设计矩阵 `~ 0 + group`，对比 `tumor - normal`
- `limma::lmFit` → `contrasts.fit` → `eBayes` → `topTable(n=Inf, adjust="BH")`
- 输出 `deg_table.csv`，**必须包含** `gene, logFC, AveExpr, t, P.Value, adj.P.Val, B`
- 显著标准：`adj.P.Val < 0.05` **且** `|log2FC| > 1`
- 火山图标注 top 基因

### 3.5 聚类热图（`04_heatmap_enrichment.R`）

- 取 top 50 显著 DEG（按 `adj.P.Val` 升序）
- **降级逻辑**：显著基因 < 50 时取全部显著基因；为 0 时取全表 top 20 并标记 `warning`
- 行聚类 euclidean + complete，表达量按行 Z-score
- 输出 `top50_heatmap.pdf`

### 3.6 GO / KEGG 富集（`04_heatmap_enrichment.R`）

- SYMBOL → ENTREZ（`org.Hs.eg.db::bitr`）
- GO BP：`clusterProfiler::enrichGO`，`pAdjustMethod="BH"`，`pvalueCutoff=0.05`
- KEGG：`clusterProfiler::enrichKEGG`，`organism="hsa"`
- dotplot 展示 top 15 条目
- **失败不终止**：KEGG REST API 有速率限制与授权限制，若返回空或报错，写入空表 +
  `results/enrichment_status.json` 记录原因，流程继续

> **背景集选择**：spec 要求"全基因组背景"，本实现遵循该默认。但芯片分析中更严谨的做法是
> 用**实测基因集**作背景（`config.yml` 中 `enrichment.universe: detected` 可切换），
> 因为未在芯片上检出的基因不应计入背景。两种结果会不同，报告中必须写明用了哪种。

### 3.7 差异基因互作（PPI，`05_ppi.R`）

- 显著 DEG symbol → `STRINGdb` 映射（`species=9606`）
- 置信度阈值 `score >= 400`
- 计算节点 degree，**hub 基因 = degree 排名前 10**
- 输出 `PPI_network.png`、`hub_genes.csv`、`ppi_edges.csv`
- **回退**：STRINGdb 需下载 ~100 MB 网络文件，若失败则跳过 PPI 并在
  `results/ppi_status.json` 标记原因，流程继续（不伪造网络图）

### 3.8 汇总（`main_analysis.R`）

按顺序编排 3.1→3.7，逐步 `tryCatch`，写 `results/state.json`：

```json
{"steps":[{"id":"deg","status":"ok","seconds":12.3},...],
 "required_failed":[],"optional_failed":["ppi"]}
```

仅当**必需步骤**（校验/清洗/QC/PCA/相关性/DEG）失败时以非零码退出。

---

## 4. 产出清单

```
results/
├── boxplot_before_after.pdf      QC：标准化前后箱线图
├── density_plot.pdf              QC：密度曲线
├── pca_plot.pdf                  PCA（含方差解释率）
├── correlation_heatmap.pdf       样本相关性热图
├── correlation_matrix.csv        Pearson + Spearman 矩阵 + 离群标记
├── deg_table.csv                 全基因差异分析表
├── volcano_plot.pdf              火山图
├── top50_heatmap.pdf             top DEG 聚类热图（Z-score）
├── GO_dotplot.pdf / GO_table.csv       GO BP 富集
├── KEGG_dotplot.pdf / KEGG_table.csv   KEGG 富集
├── PPI_network.png / hub_genes.csv     STRING PPI 与 hub 基因
├── enrichment_status.json        富集是否为空及原因
├── ppi_status.json               PPI 是否回退及原因
└── state.json                    各步骤执行状态
```

---

## 5. 已知局限（必须随结果一并报告）

1. **样本量**：n=9（6 vs 3）。所有 p 值都不稳健，**不能作为临床或机制结论**。
2. **组织成分混杂**：肿瘤 vs 全组织正常，差异主要反映细胞组成而非肿瘤特异性表达。
3. **无独立验证队列**：本设计不含验证集，DEG 未经任何外部数据复现。
4. **HER2 亚型混杂**：tumor 组内 3 HER2+ / 3 HER2−，增加了组内方差，降低了检出功效。
5. **单平台单批次**：无法评估批次效应，也无法做跨平台一致性检验。
6. **富集分析**：KEGG 依赖在线 API，可能因限流返回空结果；此时结论中不得声称"无 KEGG 通路富集"，
   只能声称"本次未获得 KEGG 结果"。

---

## 6. 复现

```bash
# GitHub Actions：手动触发或 push 触发
gh workflow run geo_analysis.yml

# 本地（需 R 4.3+ 与 Bioconductor）
Rscript scripts/main_analysis.R --config assets/config.yml
```

完整配置见 `assets/config.yml`；运行期故障排查见 `references/troubleshooting.md`。
