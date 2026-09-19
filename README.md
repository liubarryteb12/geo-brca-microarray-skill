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
| `boxplot_before_after.pdf` | 标准化前后表达分布箱线图 |
| `density_plot.pdf` | 标准化前后密度曲线 |
| `pca_plot.pdf` | PCA（PC1/PC2 方差解释率） |
| `correlation_heatmap.pdf` | 样本间 Pearson 相关热图 |
| `correlation_matrix.csv` | Pearson + Spearman 矩阵 + 离群标记 |
| `deg_table.csv` | 全基因 limma 结果（gene/logFC/P.Value/adj.P.Val） |
| `volcano_plot.pdf` | 火山图（标注 top 基因） |
| `top50_heatmap.pdf` | top DEG 聚类热图（行 Z-score，euclidean + complete） |
| `GO_dotplot.pdf` / `GO_table.csv` | GO BP 富集 |
| `KEGG_dotplot.pdf` / `KEGG_table.csv` | KEGG 通路富集 |
| `PPI_network.png` / `hub_genes.csv` / `ppi_edges.csv` | STRING PPI 网络与 hub 基因 |
| `enrichment_status.json` | 富集模式（`fdr` / `ranked_fallback`）与原因 |
| `ppi_status.json` | PPI 方法、节点边数与回退原因 |
| `state.json` | 每步执行状态 + 验收结果 |

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
