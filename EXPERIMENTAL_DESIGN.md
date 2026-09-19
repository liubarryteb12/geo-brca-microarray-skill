# GEO 乳腺癌小样本芯片数据挖掘 — 实验设计

> 数据源：**GSE92252** ｜ 平台：**GPL16025**（NimbleGen Homo sapiens Expression Array [100718_HG18_opt_expr]，45,033 探针）
> 物种：*Homo sapiens* ｜ 类型：Expression profiling by array ｜ 样本量：**9**
> 运行环境：GitHub Actions `ubuntu-latest` ｜ 语言：R / Bioconductor

---

## 1. 为什么是 GSE92252

### 1.1 合规性核验（来自 GEO 元数据，非推测）

| 约束 | 要求 | GSE92252 实测 | 结论 |
| --- | --- | --- | --- |
| 物种 | Homo sapiens | `Homo sapiens` (taxid 9606) | ✅ |
| 数据类型 | 基因芯片 | `Expression profiling by array` / GPL16025 单色 NimbleGen 芯片 | ✅ |
| 疾病 | 乳腺癌 | AR+/ER−/PR− 乳腺癌组织 vs 正常乳腺组织 | ✅ |
| 样本量 | < 10 | **9**（6 肿瘤 + 3 正常） | ✅ |
| 分组可比 | 需两组 | tumor 6 vs normal 3，同一平台同一批次 | ✅ |
| 本地算力 | 轻量 | 9×45K 矩阵，峰值内存 < 1 GB，全程 < 10 min | ✅ |

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
| 单色 vs 双色 | 按单色处理：表达值中位数 > 50 时才补做 log2 | GPL16025 是 NimbleGen 单色芯片，`exprs()` 返回的是 log2 强度而非 ratio |
| 探针 → 基因 | **三级降级映射**，见 §2.5 | 该平台注释只有 ID / GB_ACC / DESCRIPTION，没有 symbol 列 |
| 多探针同基因 | 取**表达方差最大**者 | 比取均值更能保留真实信号，且避免稀释 |
| 缺失值 | KNN 填补（`impute::impute.knn`, k=10） | 芯片常见；记录填补比例 |
| 标准化 | `limma::normalizeBetweenArrays(method="quantile")` | 跨样本可比；QC 保留 before/after 对照 |
| 过滤 | 去除全 NA / 无变异特征（有效值 < 2 或标准差为 0） | 这类特征对任何下游统计都无贡献 |
| 富集背景 | 默认**全基因组**（OrgDb），可切换为实测基因集 | 见 §3.6 |

### 2.5 GPL16025 的注释现实与三级映射

**这是本项目最容易踩空的地方。** GPL16025 的注释表只有三列：

```text
ID            GB_ACC        DESCRIPTION
AB000409      AB000409      MAP kinase interacting serine/threonine kinase 1
```

没有 `GENE_SYMBOL`，也没有 GEO curated 注释（`GPL16025.annot.gz` 返回 404）。
只写「取 symbol 列」的流水线会在下载完成后直接失败。因此 `01_download_clean.R`
实现多级降级，逐级实测覆盖率后取最好的一条：

| 级别 | 途径 | 说明 |
| --- | --- | --- |
| 1 | 平台注释的 symbol 列 | 通用路径，GPL16025 上不可用 |
| 2 | `GB_ACC` → `org.Hs.eg.db` 的 `ACCNUM` | GenBank accession，需剥掉 `.1` 之类的版本后缀 |
| 3 | `GB_ACC` 中的 RefSeq 子集 → `REFSEQ` | `NM_` / `NR_` / `XM_` / `XR_` 开头的记录 |
| 4 | `DESCRIPTION` → `GENENAME`（精确） | 注释里存的是**基因全名**而非 symbol，正好对应 GENENAME |
| 5 | `DESCRIPTION` → `GENENAME`（归一化） | 见下 |

**第 5 级是必需的，不是锦上添花。** 该平台的设计年代是 2007 年前后，DESCRIPTION
用的是当时的基因名，与今天的 GENENAME 经常只差标点或一个括号补充：

| DESCRIPTION（平台） | GENENAME（当前） | 精确 | 归一化 |
| --- | --- | --- | --- |
| `SH3-domain binding protein 2` | `SH3 domain binding protein 2` | ✗ | ✓ |
| `Rap guanine nucleotide exchange factor (GEF) 2` | `Rap guanine nucleotide exchange factor 2` | ✗ | ✓ |
| `solute carrier family 15 (oligopeptide transporter), member 1` | `solute carrier family 15 member 1` | ✗ | ✓ |

归一化 = 去括号内容 → 去所有非字母数字 → 转小写。所有级别的覆盖率都会打印到日志，
最终采用哪条、覆盖率多少，记录在 `data/clean_stats.json` 与 `data/feature_mode.json`。

**另外，不要用 GEOquery 的 `getGPL=TRUE` 取这个平台的注释。** 那条路会下载
`GPL16025_family.soft.gz`（**182 MB**，含该平台上千个 GSM 的完整记录）并在 R 里解析；
而 GEO 的 CGI `view=full` 返回**同样完整的 45,033 行**注释表，只有 **2.6 MB**。
流水线走后者，并缓存为 RDS。

全部途径都达不到 **50% 覆盖率**时，流水线**不报错**，而是退回探针层面：
QC / PCA / 相关性 / limma DEG / 热图全部照常产出（这些不依赖基因身份），
但 GO/KEGG 会被跳过，原因写入 `results/enrichment_status.json`；
STRING 查询也跳过，PPI 直接走共表达回退，原因写入 `results/ppi_status.json`。

**报告结论前必须先看这两个文件** —— 「做了 GO 富集」和「因为映射不到 symbol 所以没做」
是两个完全不同的结论。

### 2.6 为什么不用 spec 里写的 R 4.3.0

spec 指定 `r-version: '4.3.0'`。首次实跑证明这个组合在 `ubuntu-latest` 上装不上包，
job 在 **Install R packages** 一步就失败，`pak::repo_status()` 显示五个 Bioconductor
仓库全部 `ok=FALSE`。根因有两条，都不是代码问题：

| 问题 | 实测 | 后果 |
| --- | --- | --- |
| `ubuntu-latest` 现在是 noble (24.04)，而 P3M 的 Linux 二进制按**当前 R** 构建 | R 4.3 请求的二进制与 noble 上实际提供的版本对不上 | 全部退化为源码编译 |
| P3M 的 Bioconductor 镜像已下线 3.18 | `bioconductor.posit.co/packages/3.18/bioc` → 404；3.19–3.24 正常 | Bioconductor 仓库不可用，安装直接失败 |

修正：`r-version: 'release'`。这样 CRAN 依赖走 P3M 的 Linux 二进制，只有 Bioconductor
包需要源码编译，而它们绝大多数是纯 R（只有 limma / impute / GOSemSim 带少量 C/Fortran），
编译量可以忽略。改完后 **Install R packages 从「77 秒即失败」变成「730 秒成功」**。

> 附带结论：P3M 的 Bioconductor 镜像**只提供源码**（路径是 `src/contrib`，没有
> `__linux__` 段），所以「Bioconductor 全二进制」这条路在 Posit 侧并不存在。
> 想再快只能换成预装包的镜像，但 `bioconductor/bioconductor_docker` 官方镜像
> 按自己的描述只装**系统依赖**、不含 R 包，换过去并不能省时间。

### 2.7 运行时间：20 分钟上限是够的

spec 要求单次运行 < 20 分钟。**冷缓存**实测（run 35409119340，成功）：

| 步骤 | 耗时 |
| --- | --- |
| Set up job + Checkout | 3 s |
| Install system dependencies | 13 s |
| Setup R | 24 s |
| Cache R library | 1 s |
| **Install R packages** | **722 s（12 min 2 s）** |
| Record environment | 8 s |
| **Run analysis pipeline** | **115 s（1 min 55 s）** |
| Upload + Summarise + Post | 6 s |
| **整轮** | **14 min 58 s** |

所以 `timeout-minutes: 20` 保持不变，**冷缓存也在预算内**。

`setup-r-dependencies` 自带的缓存对本仓库**不生效**（没有 `DESCRIPTION` / lockfile，
它的 post 步骤实测 0 s 直接跳过），所以额外加了一步显式 `actions/cache` 缓存
`R_LIBS_USER`。缓存的 key 用固定的 `rlib-Linux-bioc-v1`，**不用 workflow 文件哈希** ——
否则每次改 workflow（哪怕只改超时或注释）都会让 12 分钟的源码编译重来一遍。
改动包列表时手动把 `bioc-v1` 递增。

> 为什么当初把超时放宽到 30 又改回来：第一次成功前只测到「装包 730 s 失败」，
> 分析耗时未知，按最坏情况估了 30 分钟。实测分析只要 115 s，20 分钟足够，
> 就按 spec 收回来了。

### 2.8 分组模式的坑（实跑踩到）

首次跑通流水线时，第 00 步报「分组失败」，3 个样本未命中任何组。原因是分组模式写成了
`"breast tumor"`，而 `tissue` 字段的实际取值并不一致：

| 组 | `characteristics_ch1` 中的 tissue 值 | 含 "breast tumor"？ |
| --- | --- | --- |
| HER2− 肿瘤（3 例） | `AR-positive, ER-/PR-negative, HER2-negative breast tumor` | ✅ |
| HER2+ 肿瘤（3 例） | `AR-positive, ER-/PR-negative, HER2-positive tumor` | ❌ **没有 "breast"** |
| 正常组织（3 例） | `Normal breast tissue` | ❌ |

正确的判别子串是 `"tumor"`（正常样本里不含该词）。**换任何数据集之前，先跑
`node scripts/find_dataset.mjs samples GSEXXXXX` 看真实的 characteristics 取值，
不要根据数据集标题推断。**

顺带一提，这个错误被第 00 步拦下了，而不是让 6 个肿瘤样本里只有 3 个进入分析 ——
这正是硬门禁存在的意义。

### 2.9 实跑发现：分组与一个全局表达位移高度混杂（**本设计最重要的限制**）

首轮成功的运行（run 35409119340）里，QC 把 **9 个样本全部**判为相关性离群，
`median_min_pearson` 只有 0.122。查相关矩阵后，结构非常清楚：

| 对比 | Pearson r |
| --- | --- |
| HER2− 三个肿瘤内部（B13/B29/B32） | 0.976 – 0.984 |
| HER2+ 三个肿瘤内部（T32/T40/T54） | 0.926 – 0.931 |
| 三个正常内部（N39/N40/N54） | 0.955 – 0.963 |
| HER2+ ↔ 正常 | 0.82 – 0.83 |
| HER2− ↔ HER2+ | **0.31 – 0.38** |
| HER2− ↔ 正常 | **0.11 – 0.19** |

**9 个样本分成三个紧致簇，而这三个簇恰好就是三个分组。** 簇内 r > 0.92，
簇间低到 0.11。

**这不是流水线引入的。** 直接对 GEO 下载的原始强度做同样的计算：

| 处理 | 簇内平均 r | 簇间平均 r |
| --- | --- | --- |
| 原始强度（不做任何处理） | 0.986 | 0.602 |
| log2(x+1) | 0.966 | 0.564 |
| log2 + quantile 标准化 | 0.966 | 0.567 |

结构在**原始沉积值**里就存在，标准化只是原样保留。也就是说这不是
log2 判断错误、不是 quantile 用错、不是探针折叠引入的。

> **用 `tools/check_sample_structure.mjs` 可以在跑 R 之前就发现这件事**（不需要 R）：
>
> ```bash
> node tools/check_sample_structure.mjs GSE92252
> ```
>
> 它输出相关矩阵、在指定阈值处切分出的相关簇，以及簇内/簇间平均相关。
> 该工具在 26,528 个探针上计算，本流水线在折叠后的 13,261 个 symbol 上计算，
> 所以**绝对数值会不同**（工具给出的 HER2− ↔ 正常是 0.30–0.37，
> 流水线是 0.11–0.19）—— quantile 标准化的参考分布取决于参与计算的行集合。
> 但两者给出的结构完全一致：三个紧致簇、簇间显著偏低。**看结构，不要纠结小数位。**

**为什么这很严重：** HER2− 与 HER2+ 都是乳腺肿瘤，两者 r 只有 0.35 在生物学上
说不通（同组织不同亚型通常 r > 0.95）。这更像是**批次效应与分组完全混杂**，
或者提交者对这批数据做了某种分组相关的处理。无论哪种，后果都一样：

> tumor-vs-normal 的差异基因表里，**分不清多少来自恶性转化、多少来自这个全局位移**。
> 1423 个显著基因不能当作肿瘤特异事件清单使用。

**旁证：** top DEG 是 `IGHG3` / `IGHG1` / `IGHG2` / `IGHV4-31`（免疫球蛋白重链）
和 `SPP1`。免疫球蛋白基因在肿瘤 vs 全组织正常里排最前，是**浆细胞浸润造成的
组织成分差异**的典型特征，与 §2.2 的警告一致，而不是肿瘤细胞内在的改变。

**流水线的处理：** 不删样本、不做批次校正（n=9、3 簇、无重复批次，
任何校正都会把分组本身一起扣掉）。改为把事实记录进 `results/qc_summary.json`：

- `mean_within_group_pearson` / `mean_between_group_pearson`
- `group_confounded_with_global_shift`（组间平均相关低于离群阈值时为 `true`）

同时 QC 会在日志里直接警告"分组与全局表达位移高度混杂"。

**结论口径：** 本流水线在 GSE92252 上产出的是**流程演示与假设生成**，
不是可引用的乳腺癌差异表达结论。任何下游解读都必须先处理这个混杂。

---

## 3. 分析流程

八个阶段，每阶段产出明确文件；任一步失败写入 `state.json` 并继续执行不依赖该步的后续步骤。

### 3.1 数据获取与校验（`00_validate_inputs.R`）

- 读取 `assets/config.yml`，校验 `dataset_id` / `group_field` / `group_values` / `contrast` 齐全
- 从 GEO SOFT 接口（base R `url()`，不依赖 Bioconductor）拉取 series 与 sample 元数据
- **硬门禁**：物种必须为 `Homo sapiens`；类型必须为 array；样本数必须 < 10；两组样本数均 ≥ 3
- 按 `group_field` 匹配 `group_values` 生成分组；一个样本命中多个组 → 报错退出
- 输出：`data/group.csv`、`data/meta.csv`、`data/platform.txt`、`data/geo_metadata.json`

### 3.2 下载与清洗（`01_download_clean.R`）

- `GEOquery::getGEO(..., getGPL = FALSE)` 只下载表达矩阵；平台注释单独用轻量 CGI
  接口取（2.6 MB，而非 182 MB 的 family 文件），详见 §2.5
- 原始强度自动识别：中位数 > 50 时执行 `log2(x + 1)`
- 探针 → 基因 symbol：多级降级映射（§2.5）；全部途径覆盖率 < 50% 时退回探针层面
- 同一 symbol 的多探针取方差最大者
- 过滤全 NA / 无变异特征（有效值 < 2 或标准差为 0）
- KNN 填补缺失（`impute::impute.knn`, k=10）
- quantile 标准化
- **一致性校验**：表达矩阵列名必须与 `group.csv` 的样本行完全对应，否则报错
- 输出：`data/expr_raw.rds`、`data/expr_clean.rds`、`data/expr_clean.csv`、
  `data/clean_stats.json`、`data/feature_mode.json`

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

- SYMBOL → ENTREZ（`clusterProfiler::bitr` + `org.Hs.eg.db`）
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
