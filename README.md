# GEO 乳腺癌小样本芯片数据挖掘流水线

对 GEO 人源乳腺癌基因表达芯片数据（**样本量 < 10**）执行端到端纯生信分析，
产出清洗、QC、PCA、样本相关性、差异基因、聚类热图、GO/KEGG 富集、
差异基因互作网络的**全部图表与表格**，并在 GitHub Actions 上运行后打包为 artifact。

**默认数据集：[GSE92252](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92252)**
人源 · Agilent GPL16025 芯片 · 9 例（6 肿瘤 vs 3 正常乳腺组织）

---

## 快速开始

```bash
# 0. 环境
#    R 4.3+ / Bioconductor；Node 18+（仅用于数据集预检）

# 1. 先验证数据集合规（不要跳过）
node scripts/find_dataset.mjs check GSE92252

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
| `PPI_network.png` / `hub_genes.csv` | STRING PPI 网络与 hub 基因 |
| `enrichment_status.json` | 富集为空时的原因 |
| `ppi_status.json` | PPI 回退原因 |
| `state.json` | 每步执行状态 + 验收结果 |

## 配置

全部参数在 [`assets/config.yml`](assets/config.yml)。切换数据集只需改这个文件：

```yaml
dataset_id: GSE92252
group_field: characteristics_ch1
group_values:
  tumor:  ["breast tumor"]
  normal: ["normal breast tissue"]
contrast: ["tumor", "normal"]      # log2FC > 0 表示 tumor 中上调
thresholds:
  adj_p: 0.05
  log2fc: 1.0
  string_score: 400
```

## 换数据集

```bash
# 找候选：人源 + 芯片 + 样本数 < 10
node scripts/find_dataset.mjs search --disease "breast cancer" --max-samples 10

# 只列出能凑出两个 >=3 样本组的数据集（慢，会逐个拉样本元数据）
node scripts/find_dataset.mjs search --disease "breast cancer" --max-samples 10 --two-groups

# 看某个数据集的样本明细
node scripts/find_dataset.mjs samples GSE92252

# 完整校验
node scripts/find_dataset.mjs check GSEXXXXX
```

然后把 `dataset_id` 和分组字段写进 `assets/config.yml`。

## 重要限制

**这是探索性分析，不是机制结论。** n=9（6 vs 3）的设计只能检出效应量很大的基因，
且肿瘤 vs 全组织正常的差异主要来自组织成分（上皮/脂肪/基质/免疫细胞），
而非肿瘤特异驱动事件。

完整的统计功效边界、分组依据与被否决的数据集，见
[`EXPERIMENTAL_DESIGN.md`](EXPERIMENTAL_DESIGN.md)。

## 故障排查

见 [`references/troubleshooting.md`](references/troubleshooting.md)。

## 许可

MIT
