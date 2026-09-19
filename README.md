# GEO 乳腺癌小样本芯片数据挖掘流水线

对 GEO 人源乳腺癌基因表达芯片数据（**样本量 < 10**）执行端到端纯生信分析，
产出清洗、QC、PCA、样本相关性、差异基因、聚类热图、GO/KEGG 富集、
差异基因互作网络的**全部图表与表格**，并在 GitHub Actions 上运行后打包为 artifact。

**默认数据集：[GSE64790](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64790)**
人源 · Agilent GPL19612 lncRNA 芯片 · 6 例（3 TNBC vs 3 **配对**正常乳腺组织）

> ✅ **已验证跑通**：[run 35412006459](https://github.com/liubarryteb12/geo-brca-microarray-skill/actions/runs/35412006459)
> —— 暖缓存 3 min 45 s，7 个步骤全部 `ok`，12 项验收全过，24 个产物。
>
> ⚠️ n=6 时**没有任何基因能通过 FDR**（最小 `adj.P` = 0.394）。富集与 PPI 走的是
> 明确标注的 `ranked_fallback` 降级路径，结果只能作假设生成。见
> [`EXPERIMENTAL_DESIGN.md`](EXPERIMENTAL_DESIGN.md) §2.2 与 §2.10。

---

## 快速开始

```bash
# 0. 环境
#    R 4.3+ / Bioconductor；Node 18+（仅用于数据集预检）

# 1. 先验证数据集合规（不要跳过）
node scripts/find_dataset.mjs check GSE64790
node tools/check_sample_structure.mjs GSE64790   # 分组是否与批次混杂

# 2. 改配置（可选）
#    数据集、分组、阈值全部在 assets/config.yml

# 3. 跑
Rscript scripts/main_analysis.R --config assets/config.yml

# 4. 看结果
ls results/
cat results/state.json
```

## GitHub Actions

```bash
gh workflow run geo_analysis.yml
gh run watch
gh run download --name geo-results
```

workflow 手动触发或 push 到 `main` 时运行，上限 20 分钟，产物上传为
artifact `geo-results`。

## 产出

| 文件 | 内容 |
| --- | --- |
| `boxplot_before_after.pdf` / `.png` | 标准化前后表达分布箱线图 |
| `density_plot.pdf` / `.png` | 标准化前后密度曲线 |
| `pca_plot.pdf` / `.png` | PCA（PC1/PC2 方差解释率） |
| `correlation_heatmap.pdf` / `.png` | 样本间 Pearson 相关热图 |
| `correlation_matrix.csv` | Pearson + Spearman 矩阵 + 离群标记 |
| `deg_table.csv` | 全基因 limma 结果（gene/logFC/P.Value/adj.P.Val） |
| `volcano_plot.pdf` / `.png` | 火山图。**颜色 = 方向**（up 红 `#B2182B` / down 蓝 `#2166AC`），**alpha + 大小 = 置信度**（FDR 显著实心大点，名义显著半透明小点）。纵轴统一 raw P |
| `top50_heatmap.pdf` / `.png` | top DEG 聚类热图（行 Z-score，euclidean + complete） |
| `pvalue_histogram.pdf` / `.png` | DE 后 QC：p 值分布（区分"功效不足"与"模型设定错"） |
| `GSEA_GO_dotplot.pdf` / `.png` + `GSEA_GO_table.csv` | **preranked GSEA / GO BP（主力方法）** |
| `GSEA_KEGG_dotplot.pdf` / `.png` + `GSEA_KEGG_table.csv` | preranked GSEA / KEGG |
| `GO_dotplot.pdf` / `.png` + `GO_table.csv` | ORA GO BP（含 `direction` 列，上/下调分开） |
| `KEGG_dotplot.pdf` / `.png` + `KEGG_table.csv` | ORA KEGG（含 `direction` 列） |
| `PPI_network.pdf` / `.png` + `hub_genes.csv` + `ppi_edges.csv` + `ppi_plot_layout.csv` | STRING PPI 网络与 hub 基因。**图做过可读性过滤，并按同心圆环排布**（最大连通分量 → degree 前 200 → 最强 450 条边 → 3 环，内圈 = hub 核心；最大 4 个 Louvain 模块上色）。完整网络见 `ppi_edges.csv`，每个节点的环号/半径/角度见 `ppi_plot_layout.csv` |
| `enrichment_status.json` | 富集模式（`fdr` / `ranked_fallback`）、GSEA 参数、去冗余阈值与原因 |
| `ppi_status.json` | PPI 方法、节点边数、绘图过滤与环参数、回退原因 |
| `state.json` | 每步执行状态 + 验收结果（唯一逐字节不可复现的产物：含耗时与时间戳） |

### 配色

所有图共用 `scripts/lib/common.R` 里的一套语义化色板，来源是 SCI 发表常用色板：

| 角色 | 来源 |
| --- | --- |
| 上调 / tumor、下调 / normal | **ColorBrewer RdBu** 两端 `#B2182B` / `#2166AC` —— 发散热图与方向色同源，热图上的红就是火山图上"上调"的那个红 |
| 显著性（连续） | **viridis**，**截去暗端**（完整 viridis 的 `#365C8D` 距 down 蓝仅 3.1°，深色点会被读成"下调"） |
| 模块 / 多组 | **Okabe-Ito** 变体穷举出的 4 色（绿 / 玫红 / 天蓝 / 橙） |

判据是量化的：**色相相差 15° 以内视为同一个颜色**。按这条算，之前有四处
凭眼睛看不出来的问题（`firebrick` 与 up 红差 0.4°、分类色板在 protanopia 下差 0.002、
magma 中段与 up 红差 2.9°、viridis 暗端与 down 蓝差 3.1°）。

```bash
node tools/check_palette.mjs   # 从 common.R 解析实际色值重算，CI 里是门禁
```

> 每张图都同时出 **PDF（矢量，放大不失真）** 和 **PNG（150 dpi，可直接预览）**。
> artifact 是 zip，PDF 在里面看不了，PNG 是为了打开就能看到。

### 富集为什么有两条路

判据来自 K-Dense `pathway-enrichment` skill：*"a discrete hit list → ORA;
a ranked table with per-gene scores → GSEA"*、*"Never threshold a list and then
feed it to GSEA"*。

| | GSEA（主力） | ORA（辅助） |
| --- | --- | --- |
| 输入 | **完整排序表**（13,948 个基因，不卡阈值） | FDR 显著基因，或降级到 raw P 前 500 |
| 排序/阈值 | limma moderated `t` | `adj.P < 0.05` 且 `\|log2FC\| > 1` |
| 方向 | NES 符号 | **上/下调分开跑**（`direction` 列） |
| 背景 | 不需要 | `detected`（实测 16,487 个基因） |

GSE64790 实测两者差距**两个数量级**：GSEA 拿到 **1,059 条**显著 GO BP 条目
（最好 `adj.P` = 1.0e-8），ORA 在 top-500 上最好只有 3.1e-6。
弱功效、效应弥散的数据正是 GSEA 被设计出来处理的场景。

两条路的结果都做**基因重叠去冗余**（Jaccard ≥ 0.5 折叠为一类），
GO 从 1059 折到 **439** 个代表条目。报告时引用 `representative` 列。

## 配置

全部参数在 [`assets/config.yml`](assets/config.yml)。切换数据集只需改这个文件：

```yaml
dataset_id: GSE64790
group_field: characteristics_ch1
group_values:
  # 判别子串实测：肿瘤是 "tissue: TNBC  tissue"（两个空格），
  # 正常是 "tissue: matched normal breast tissues"
  tumor:  ["tnbc"]
  normal: ["matched normal"]
contrast: ["tumor", "normal"]      # log2FC > 0 表示 tumor 中上调

# 配对关系显式声明，不解析标题后缀；00 步会逐对校验
paired: true
pairs:
  - ["GSM1580581", "GSM1580584"]   # 72y
  - ["GSM1580582", "GSM1580585"]   # 41y
  - ["GSM1580583", "GSM1580586"]   # 52y

thresholds:
  adj_p: 0.05
  log2fc: 1.0
  string_score: 400
analysis:
  ranked_fallback_genes: 500       # 无 FDR 显著基因时富集/PPI 的降级输入（0 = 关闭）
```

## 换数据集

```bash
# 找候选：人源 + 芯片 + 样本数 < 10
node scripts/find_dataset.mjs search --disease "breast cancer" --max-samples 10

# 只列出能凑出两个 >=3 样本组的数据集（慢，会逐个拉样本元数据）
node scripts/find_dataset.mjs search --disease "breast cancer" --max-samples 10 --two-groups

# 看某个数据集的样本明细
node scripts/find_dataset.mjs samples GSE64790

# 完整校验
node scripts/find_dataset.mjs check GSEXXXXX

# 查相关结构：分组是否与批次效应混杂（这一步能省掉一整轮 15 分钟的失败运行）
node tools/check_sample_structure.mjs GSEXXXXX
```

然后把 `dataset_id` 和分组字段写进 `assets/config.yml`。

> **换数据集时最容易漏的两件事**：① 分组与批次混杂（第一版 GSE92252 就栽在这里，
> 三个分组恰好是三个表达批次）；② 平台注释里根本没有基因 symbol
> （GSE112848 形式条件全过，但平台只有 5 列注释，富集和 PPI 会全部跳过）。

## 重要限制

**这是探索性分析，不是机制结论。** n=6（3 vs 3 配对，残差 df = 2）的设计下，
**没有任何基因能通过全基因组 FDR 校正** —— 1,456 个基因在 `raw P < 0.05` 且
`|log2FC| > 1` 水平上显著，但最小 `adj.P` 是 0.394。这是"样本 < 10"这个要求本身的
固有限制，不是分析错误。

因此富集与 PPI 走**明确标注的降级路径**：输入改为按 `raw P` 排序的前 500 个基因，
`enrichment_status.json` / `ppi_status.json` 里的 `deg_mode` 会写成 `ranked_fallback`。
**此时只能说"在最显著的 N 个基因里富集到……"，不能说"显著差异基因富集到……"。**

此外，肿瘤 vs 全组织正常的差异主要来自组织成分（上皮/脂肪/基质/免疫细胞），
而非肿瘤特异驱动事件。

完整的统计功效边界、分组依据与被否决的数据集，见
[`EXPERIMENTAL_DESIGN.md`](EXPERIMENTAL_DESIGN.md)。

## 故障排查

见 [`references/troubleshooting.md`](references/troubleshooting.md)。

## 许可

MIT
