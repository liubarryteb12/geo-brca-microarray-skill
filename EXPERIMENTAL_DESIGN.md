# GEO 乳腺癌芯片数据挖掘 — 实验设计

> **本文档的路径约定：** 文中出现的 `results/X` 与 `data/X` 一律指
> `results/<GSE>/X` 与 `data/<GSE>/X` —— 产物**按数据集分目录**，
> 跑第二个数据集不会覆盖第一个。配置也一个数据集一个文件
> （`assets/config.<GSE>.yml`）。目录名由 `dataset_id` 在 `load_config()` 里派生。

> 数据源：**GSE64790** ｜ 平台：**GPL19612**（Agilent-062918 OE Human lncRNA Microarray V4.0，111,088 探针）
> 物种：*Homo sapiens* ｜ 类型：Expression profiling by array ｜ 样本量：**6**（3 例 TNBC + 3 例配对正常）
> `design_mode: small_sample` ｜ 运行环境：GitHub Actions `ubuntu-latest` ｜ 语言：R / Bioconductor
>
> 另有队列级设计：**GSE42568**（121 例，`design_mode: cohort`），见 §6。

---

## 1. 为什么是 GSE64790

### 1.1 合规性核验（来自 GEO 元数据，非推测）

| 约束 | 要求 | GSE64790 实测 | 结论 |
| --- | --- | --- | --- |
| 物种 | Homo sapiens | `Homo sapiens` (taxid 9606) | ✅ |
| 数据类型 | 基因芯片 | `Expression profiling by array` / GPL19612 | ✅ |
| 疾病 | 乳腺癌 | 三阴性乳腺癌（TNBC）组织 vs 配对正常乳腺组织 | ✅ |
| 样本量 | < 10 | **6**（3 肿瘤 + 3 正常） | ✅ |
| 每组样本数 | ≥ 3 | 3 vs 3 | ✅ |
| 分组可比 | 需两组 | 同平台、同批次，全部 6 例两两 r = 0.914–0.967 | ✅ |
| 本地算力 | 轻量 | 6×65K 矩阵，全程约 106 s | ✅ |

### 1.2 候选筛选过程（这是本设计最关键的一步）

**"样本 < 10 + 乳腺癌 + 人源 + 芯片"这四个条件本身很容易满足，难的是找到一个
分组没有被批次效应混杂的数据集。** 第一版选的 GSE92252 满足全部形式条件，
实跑后才发现三个分组恰好是三个表达批次（簇间 r 低到 0.11），
DEG 结果无法解释 —— 详见 §2.9。

所以这次先写筛选器，对 GEO 中所有"人源 + 乳腺癌 + 芯片 + 4 ≤ n ≤ 9"的 series
逐个拉 GSM 元数据，要求：恰好两组、每组 ≥ 3 例、组织样本（排除细胞系）、
无处理/转染。约 200 个候选里只剩 **4 个**是真正的组织样本两组设计：

| 数据集 | 分组 | 平台 | 实测结论 |
| --- | --- | --- | --- |
| **GSE64790** | TNBC 3 vs 配对正常 3 | GPL19612 | ✅ **采用**。6 例两两 r = 0.914–0.967，单一相关簇 |
| GSE112848 | DCIS 3 vs 良性组织 3 | GPL16956 | ❌ 平台是 Arraystar lncRNA V3，注释**只有 5 列**（ID/类型/BUILD/序列/SPOT_ID），**没有任何基因注释**，无法做富集 |
| GSE73613 | 浸润性癌 2 vs 正常 2 | GPL570 | ❌ 每组只有 2 例，达不到"每组 ≥ 3"的硬门禁 |
| GSE207304 | 乳腺癌细胞外泌体 3 vs 正常 3 | GPL26963 | ❌ 外泌体，非组织 |

> 教训一：**形式合规 ≠ 可用。** 类型、物种、样本量都能过门禁，但如果分组与批次
> 共线，产出的 DEG 表就是不可解释的。选数据集时必须先查相关结构，
> `tools/check_sample_structure.mjs` 就是干这个的。
>
> 教训二：**平台注释要先查。** GSE112848 看着完美（3 vs 3、组织、配对），
> 但 GPL16956 根本没有基因注释列，流水线会退回探针模式，GO/KEGG/STRING 全部跳过。

### 1.3 样本构成

| GSM | 标题 | 分组 | 患者 | 年龄 |
| --- | --- | --- | --- | --- |
| GSM1580581 | TNBC tissue 1 | tumor | pair01 | 72y |
| GSM1580582 | TNBC tissue 2 | tumor | pair02 | 41y |
| GSM1580583 | TNBC tissue 3 | tumor | pair03 | 52y |
| GSM1580584 | matched normal breast tissues 1 | normal | pair01 | 72y |
| GSM1580585 | matched normal breast tissues 2 | normal | pair02 | 41y |
| GSM1580586 | matched normal breast tissues 3 | normal | pair03 | 52y |

配对依据是**年龄一一对应**（72/72、41/41、52/52），且提交者在 `Series_overall_design`
里明确写了 "matched histological normal breast tissues"。
配对关系在 `assets/config.<GSE>.yml` 里**显式声明**，不靠解析标题后缀。

分组依据 `characteristics_ch1` 的 `tissue:` 字段：肿瘤含 `TNBC  tissue`，
正常含 `matched normal breast tissues`。

### 1.4 被否决的初始候选

最初提出的四个数据集**全部不合规**，若直接使用，流水线会在第一步
`validate_inputs` 就终止：

| 数据集 | 物种 | 数据类型 | 样本数 | 否决原因 |
| --- | --- | --- | --- | --- |
| **GSE197894** | 人 ✅ | **RNA-seq** ❌ | **20** ❌ | 新辅助化疗前后配对 RNA-seq，非芯片；样本超限；分组不是肿瘤 vs 正常 |
| GSE143334 | **小鼠** ❌ | ChIP-seq ❌ | 6 | 骨骼肌昼夜节律，与乳腺癌无关 |
| GSE161922 | **小鼠** ❌ | RNA-seq ❌ | 6 | B16 黑色素瘤模型 |
| GSE206338 | **小鼠** ❌ | 表达谱 ❌ | — | 小鼠 NASH 肝纤维化模型 |

> 教训：GSE197894 的检索摘要被描述为"表达谱芯片"，但其 `gdstype` 实际是
> `Expression profiling by high throughput sequencing`。**分组与类型必须以 GEO 的
> `gdstype` / `GPL` 字段为准，不能采信标题或二手描述。**

---

## 2. 设计决策与理由

### 2.1 主对比：tumor (3) vs normal (3)，**配对**

`contrast = ["tumor", "normal"]`，limma 中 tumor 为分子、normal 为分母，
即 **log2FC > 0 表示肿瘤中上调**。

该数据集是"同一患者的肿瘤 + 正常组织"配对设计，所以做配对分析：
设计矩阵 `~ 0 + groups + patient`，患者作为阻断因子。

**为什么必须配对**：患者间差异往往比肿瘤/正常差异还大。不阻断的话，这部分方差
全部落进残差，真正的信号会被埋掉。配对后残差 df = 6 − 4 = **2**。

配对关系显式写在 `assets/config.<GSE>.yml` 的 `pairs` 里，`00_validate_inputs.R` 会校验：
每个 GSM 只出现一次、每对必须一例 tumor 一例 normal、不允许有样本落单。
`03_deg.R` 若发现配对设计不可用（秩不足等）会**自动退回非配对**并把原因写进
`deg_summary.json` 的 `paired_fallback_reason`，绝不让它拖垮整条流水线。

### 2.2 统计功效：这是本设计必须写在最前面的限制

3 vs 3 配对（残差 df = 2）是一个**极小样本设计**，实跑结果把它的边界暴露得非常清楚：

| 指标 | GSE64790 实测 |
| --- | --- |
| 检验基因数 | 16,487 |
| `raw P < 0.05` 且 `\|log2FC\| > 1` | **1,456** 个 |
| `raw P < 0.001` | 36 个 |
| 最小 `raw P` | 6.0e-5（KRT14） |
| 最小 `adj.P` | **0.394** |
| `adj.P < 0.05` 且 `\|log2FC\| > 1` | **0 个** |

**信号是真实的**：1,456 个基因在 raw P < 0.05 水平上显著，远超随机预期的 5%；
top 命中（KRT14、SPARCL1、TAGLN、SDPR、PPARG、PGR）也全是教科书级的
肿瘤 vs 正常乳腺组织基因。**但全基因组 BH 校正过不去** ——
在 16,487 个基因上做 BH，最小的 raw P 需要达到约 3e-6 才能得到 adj.P < 0.05，
而实测最小值是 6e-5，差了 20 倍。

> **这是"样本 < 10"这个要求本身的固有限制，不是分析错误。**
> 任何 n < 10 的乳腺癌全基因组芯片数据集都会撞上同一堵墙。

#### 2.2.1 在既定 n 下的最小可检测效应（MDE）

上面是**多重检验**的论证。还缺**功效**的论证：这个设计**能够**检出多大的效应？
用 K-Dense `statistical-power` skill 的 `mde()` / `power()` 计算
（`t_paired`，3 对，df=2，α=0.05 双侧）：

| 目标功效 | MDE（Cohen's d_z） |
| --- | --- |
| 80% | **3.26** |
| 90% | 3.92 |

即：只有配对差值的效应量达到 **3.26 个标准差**，本设计才有 80% 把握检出。
n=3 时各真实效应量对应的实际功效：

| 真实 d_z | 0.5 | 1.0 | 1.5 | 2.0 | 2.5 | 3.0 | 4.0 | 5.0 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 功效 | 0.08 | 0.18 | 0.32 | 0.47 | 0.62 | 0.75 | 0.91 | 0.98 |

要达到 80% 功效：d_z=2.0 需 **5 对**，d_z=1.5 需 **6 对**，d_z=1.0 需 **10 对**，
d_z=0.5 需 **34 对**。GSE64790 只有 3 对，连 d_z=1.0 这种"中等偏大"的效应也只有
18% 把握检出。

> **这不是 post-hoc observed power。** 用观测到的效应反算功效是循环论证 ——
> 它是 p 值的确定性函数，不提供任何新信息，是审稿人常见的批评点。
> 这里报的是**敏感性分析**：在固定 n 下，本设计**能够**检出多大的效应。

> **关于配对的实际收益**：n=3 时配对检验的 MDE 是 3.26，而同样 3 vs 3 **不配对**
> 是 3.07 —— **配对在这个样本量下没有换来功效优势**。原因是 df 从 4 降到 2，
> t 临界值从 2.78 升到 4.30，这个代价超过了配对降低方差的收益。
> 配对仍然是**正确的**（它控制患者个体基线，且配对差值本身方差更小，
> 同样的生物学差异会映射成更大的 d_z），但**不要指望 3 对配对能提升检出能力**。

因此下游（富集、PPI）走**明确标注的降级路径**，见 §2.10。

**解释陷阱**：肿瘤组织与正常乳腺组织的差异，绝大部分来自
**组织成分差异**（上皮比例、脂肪、基质、浸润免疫细胞），而非肿瘤特异性驱动事件。

> **判读边界**：本流水线产出的是**探索性假设**，不是肿瘤发生机制结论。
> 任何"某基因驱动乳腺癌"的表述都超出了本设计能支持的范围。
> 若要区分"肿瘤特异"与"组织成分"，需要 LCM 显微切割数据或去卷积（如 CIBERSORTx）。

### 2.3 为什么不做亚型次级对比

TNBC 组内部只有 3 例，没有可分的次级结构。即使有，每组 n=3 时 limma 也无法给出
可信的 `adj.P`。次级对比不纳入本设计。

### 2.4 芯片特有的技术决策

| 问题 | 决策 | 理由 |
| --- | --- | --- |
| 单色 vs 双色 | 按单色处理：表达值中位数 > 50 时才补做 log2 | GPL19612 矩阵已是 log2 尺度（中位数约 5.4，最大 19），实测**不会**触发 log2；若误取 log2 会把信号压平 |
| 探针 → 基因 | 优先平台注释的 `GeneSymbol` 列，见 §2.5 | GPL19612 有 23 列注释，含 `GeneSymbol` / `GenbankAccession` / `GeneName` |
| 多探针同基因 | 取**表达方差最大**者 | 比取均值更能保留真实信号，且避免稀释 |
| 缺失值 | KNN 填补（`impute::impute.knn`, k=10） | 实测 0 个缺失值，仍保留该步 |
| 标准化 | `limma::normalizeBetweenArrays(method="quantile")` | 跨样本可比；QC 保留 before/after 对照 |
| 过滤 | 去除全 NA / 无变异特征（有效值 < 2 或标准差为 0） | 这类特征对任何下游统计都无贡献 |
| 富集背景 | 默认**全基因组**（OrgDb），可切换为实测基因集 | 见 §3.6 |

### 2.5 GPL19612 的注释与映射判据

GPL19612（Agilent-062918 OE Human lncRNA V4.0）有 111,088 行注释、23 列，
其中 `GeneSymbol`、`GenbankAccession`、`GB_ACC`、`GeneName` 都可用于映射。
矩阵里 65,531 个探针，逐列实测覆盖率：

| 途径 | 覆盖探针数 | 覆盖率 |
| --- | --- | --- |
| `GeneSymbol` | 21,812 | 33.3% |
| `GenbankAccession` | 22,492 | 34.3% |
| `GB_ACC` | 19,675 | 30.0% |
| `GeneName` | 21,650 | 33.0% |

`GeneSymbol` 被选中，最终得到 **16,487 个唯一基因**。

**这里有个判据陷阱，值得单独记下来。** 早期版本只按"注释覆盖率 ≥ 50%"决定
是否走 symbol 模式。这在 lncRNA / 外显子芯片上必然误判：这类芯片**大部分探针
本来就是非编码的**，覆盖率天然上不去。GSE64790 只有 33.3%，
按老判据会被降级成探针模式，GO/KEGG/STRING 全部跳过 —— 而它实际有
16,487 个带 symbol 的基因，做基因层面分析绰绰有余。

现在的判据是**覆盖率与基因数取「或」**：

```r
if (best$coverage < MIN_SYMBOL_COVERAGE && n_genes < MIN_SYMBOL_GENES)  # 50% 且 5000 个
```

覆盖率低但基因数够 → 照常做基因层面分析；两者都不够才退回探针模式。

### 2.6 多级降级映射

`01_download_clean.R` 实现多级降级，逐级实测覆盖率后取最好的一条：

| 级别 | 途径 | 说明 |
| --- | --- | --- |
| 1 | 平台注释的 symbol 列 | 通用路径；GPL19612 走的就是这条（`GeneSymbol`） |
| 2 | `GB_ACC` → `org.Hs.eg.db` 的 `ACCNUM` | GenBank accession，需剥掉 `.1` 之类的版本后缀 |
| 3 | `GB_ACC` 中的 RefSeq 子集 → `REFSEQ` | `NM_` / `NR_` / `XM_` / `XR_` 开头的记录 |
| 4 | `DESCRIPTION` → `GENENAME`（精确） | 注释里存的是**基因全名**而非 symbol，正好对应 GENENAME |
| 5 | `DESCRIPTION` → `GENENAME`（归一化） | 见下 |

**这套降级不是为 GPL19612 写的，而是为 GPL16025 那类"注释里没有 symbol"的平台准备的。**
该平台（NimbleGen 100718_HG18_opt_expr）的注释表只有 `ID / GB_ACC / DESCRIPTION` 三列，
GEO curated 注释也返回 404，只写「取 symbol 列」的流水线会在下载完成后直接失败。
它的 DESCRIPTION 用的是 2007 年前后的基因名，与今天的 GENENAME 经常只差标点或括号补充：

| DESCRIPTION（平台） | GENENAME（当前） | 精确 | 归一化 |
| --- | --- | --- | --- |
| `SH3-domain binding protein 2` | `SH3 domain binding protein 2` | ✗ | ✓ |
| `Rap guanine nucleotide exchange factor (GEF) 2` | `Rap guanine nucleotide exchange factor 2` | ✗ | ✓ |
| `solute carrier family 15 (oligopeptide transporter), member 1` | `solute carrier family 15 member 1` | ✗ | ✓ |

归一化 = 去括号内容 → 去所有非字母数字 → 转小写。所有级别的覆盖率与**唯一 symbol 数**
都会打印到日志，最终采用哪条、覆盖率多少，记录在 `data/clean_stats.json` 与
`data/feature_mode.json`。

**另外，不要用 GEOquery 的 `getGPL=TRUE` 取平台注释。** 对 GPL16025 那条路会下载
`GPL16025_family.soft.gz`（**182 MB**，含该平台上千个 GSM 的完整记录）并在 R 里解析；
而 GEO 的 CGI `view=full` 返回**同样完整的 45,033 行**注释表，只有 **2.6 MB**。
流水线走后者，并缓存为 RDS。

**什么时候退回探针层面**：覆盖率与唯一基因数**都不够**时（`< 50%` **且** `< 5000 个`）。
此时流水线**不报错**，而是退回探针层面：
QC / PCA / 相关性 / limma DEG / 热图全部照常产出（这些不依赖基因身份），
但 GO/KEGG 会被跳过，原因写入 `results/enrichment_status.json`；
STRING 查询也跳过，PPI 直接走共表达回退，原因写入 `results/ppi_status.json`。

> **为什么是「或」而不是只看覆盖率**：lncRNA / 外显子芯片上大部分探针本来就是非编码的，
> 覆盖率天然上不去。GPL19612 只有 33.3%，但对应 **16,487 个基因**，
> 做基因层面分析绰绰有余。只卡覆盖率会把这类平台误判成探针模式，下游全部跳过。

**报告结论前必须先看这两个文件** —— 「做了 GO 富集」和「因为映射不到 symbol 所以没做」
是两个完全不同的结论。

### 2.7 为什么不用 spec 里写的 R 4.3.0

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

| 步骤 | 冷缓存 | 暖缓存（run 35412006459） |
| --- | --- | --- |
| Set up job + Checkout | 3 s | 4 s |
| Install system dependencies | 13 s | 17 s |
| Setup R | 24 s | 29 s |
| Cache R library | 1 s | 4 s |
| **Install R packages** | **722 s（12 min 2 s）** | **42 s** |
| Record environment | 8 s | 8 s |
| **Run analysis pipeline** | **115 s** | **114 s** |
| Upload + Summarise + Post | 6 s | 3 s |
| **整轮** | **14 min 58 s** | **3 min 45 s** |

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

> 换到 GSE64790 时同样要小心：肿瘤是 `tissue: TNBC  tissue`（注意**两个空格**），
> 正常是 `tissue: matched normal breast tissues`。判别子串用 `"tnbc"` 和
> `"matched normal"`，已用 `find_dataset.mjs check` 实测 6/6 命中且无歧义。

### 2.9 为什么否决 GSE92252（**换数据集的原因**）

第一版设计用的是 GSE92252（9 例：6 肿瘤 + 3 正常）。它满足全部形式条件 ——
人源、芯片、乳腺癌、n < 10、两组 —— 但实跑后发现**三个分组恰好是三个表达批次**：

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
tumor-vs-normal 的差异基因表里，**分不清多少来自恶性转化、多少来自这个全局位移**。

**旁证：** top DEG 是 `IGHG3` / `IGHG1` / `IGHG2` / `IGHV4-31`（免疫球蛋白重链）
和 `SPP1`。免疫球蛋白基因在肿瘤 vs 全组织正常里排最前，是**浆细胞浸润造成的
组织成分差异**的典型特征，与 §2.2 的警告一致，而不是肿瘤细胞内在的改变。

**所以这个数据集被换掉了。** 不删样本、不做批次校正（n=9、3 簇、无重复批次，
任何校正都会把分组本身一起扣掉）—— 这些补救手段都救不了一个分组与批次共线的设计。
改为在筛选阶段就用 `tools/check_sample_structure.mjs` 把关，见 §1.2。

**作为通用防线保留：** 流水线仍会把相关结构记进 `results/qc_summary.json`：

- `mean_within_group_pearson` / `mean_between_group_pearson`
- `group_confounded_with_global_shift`（组间平均相关低于离群阈值时为 `true`）

同时 QC 会在日志里直接警告"分组与全局表达位移高度混杂"。
在 GSE64790 上这三项分别是 0.9082 / 0.894 / **false** —— 干净的对照。

### 2.10 无 FDR 显著基因时的降级路径（GSE64790 实际走的就是这条）

§2.2 说明：n=6 时全基因组 BH 校正几乎不可能有基因通过。
但 spec 要求产出 GO/KEGG 富集和 PPI 网络。**如果因为没有显著基因就把这两步跳过，
交付物就只剩一半。** 所以下游走一条**明确标注的降级路径**：

| 情形 | 下游输入 | `deg_mode` |
| --- | --- | --- |
| FDR 显著基因 ≥ 5 个 | `adj.P < 0.05` 且 `\|log2FC\| > 1` 的基因 | `fdr` |
| FDR 显著基因 < 5 个 | 按 `raw P` 排序、`\|log2FC\| > 1` 的前 500 个基因 | `ranked_fallback` |

`ranked_fallback_genes` 在 `assets/config.<GSE>.yml` 里可调；设为 `0` 即关闭降级
（富集/PPI 直接跳过并记录原因）。

**降级必须被标注，这是硬规则。** `enrichment_status.json` 和 `ppi_status.json`
都会写入完整的 `deg_mode` 与 `deg_reason`，例如：

```json
"deg_mode": "ranked_fallback",
"deg_reason": "FDR 显著基因仅 0 个（需 >= 5）。样本数 6 下对 16487 个基因做
               BH 校正过严，最小的 adj.P 为 0.394。退回按 raw P 排序、
               |log2FC| > 1 的前 500 个基因。**这是假设生成，不是显著差异基因清单。**"
```

> **结论口径**：`ranked_fallback` 模式下产出的富集条目与 PPI hub 基因，
> 只能表述为"在最显著的 500 个基因里富集到……"，**不得**写成
> "显著差异基因富集到……"。这与 AGENTS.md 的硬规则 2 是同一条要求。

**GSE64790 实测结果**（run 35412006459）：

- 富集输入 500 个基因 → GO 44 条（`mitotic cell cycle phase transition`、
  `DNA replication`、`meiotic spindle assembly`…），KEGG 14 条
  （`Integrin signaling`、`PI3K-Akt signaling pathway`、`Mismatch repair`…）
- PPI 输入 500 个基因，STRING v12.0 映射成功 485 个 → 390 节点 / 3,284 边，
  hub 基因为 `GAPDH`、`CD34`、`IGF1`、`CHEK1`、`CCNB2`、`CXCL12`、`EXO1`、
  `CENPA`、`MAD2L1`

这些结果在生物学上自洽：TNBC 相对正常乳腺组织，增殖/细胞周期通路上调
（mitotic、DNA replication、CHEK1、CCNB2、MAD2L1、CENPA、EXO1），
基质与血管相关基因下调（CD34、IGF1、CXCL12、KRT14、TAGLN、SDPR、PPARG）。
**方向合理，但按上面的口径，这仍只是假设生成。**

---

## 3. 分析流程

八个阶段，每阶段产出明确文件；任一步失败写入 `state.json` 并继续执行不依赖该步的后续步骤。

### 3.1 数据获取与校验（`00_validate_inputs.R`）

- 读取 `assets/config.<GSE>.yml`，校验 `dataset_id` / `group_field` / `group_values` / `contrast` 齐全
- 从 GEO SOFT 接口（base R `url()`，不依赖 Bioconductor）拉取 series 与 sample 元数据
- **硬门禁**：物种必须为 `Homo sapiens`；类型必须为 array；样本数必须 < 10；两组样本数均 ≥ 3
- 按 `group_field` 匹配 `group_values` 生成分组；一个样本命中多个组 → 报错退出
- `paired: true` 时校验 `pairs`：每个 GSM 只出现一次、每对必须一例 numerator 一例
  denominator、不允许有样本落单（§2.1）
- 输出：`data/group.csv`（含 `patient` 列）、`data/meta.csv`、`data/platform.txt`、
  `data/geo_metadata.json`

### 3.2 下载与清洗（`01_download_clean.R`）

- `GEOquery::getGEO(..., getGPL = FALSE)` 只下载表达矩阵；平台注释单独用轻量 CGI
  接口取（远小于 family 文件），详见 §2.5
- 原始强度自动识别：中位数 > 50 时执行 `log2(x + 1)`
- 探针 → 基因 symbol：多级降级映射（§2.6）；覆盖率与基因数**都不够**时才退回探针层面
- 同一 symbol 的多探针取方差最大者
- 过滤全 NA / 无变异特征（有效值 < 2 或标准差为 0）
- KNN 填补缺失（`impute::impute.knn`, k=10），**调用前 `set.seed(cfg$analysis$seed)`**
- quantile 标准化
- **一致性校验**：表达矩阵列名必须与 `group.csv` 的样本行完全对应，否则报错
- 输出：`data/expr_raw.rds`、`data/expr_clean.rds`、`data/expr_clean.csv`、
  `data/clean_stats.json`、`data/feature_mode.json`

#### 3.2.1 可复现性：五个随机源，加一个不是随机源的问题

这条是被实测逼出来的。连续两轮 CI 的 GSEA 显著 GO 条目数不同（1059 / 1110），
看起来像蒙特卡洛的随机波动，实际上根因不在 GSEA。

**四个"没设种子"的随机源：**

| 随机源 | 症状 | 处理 |
| --- | --- | --- |
| `impute::impute.knn` 用 `sample()` 处理并列近邻 | 插补值每次不同 → t 统计量 → GSEA 排序 → 所有 p 值全部分叉 | 调用前 `set.seed()` |
| `fgsea` 的 `fgseaMultilevel` 自适应采样 | NES / p 值在第 4 位有效数字分叉，条目集合都不同（1072 vs 1095） | 调用前 `set.seed()` |
| `igraph::layout_with_fr()` 随机初始位置 | `PPI_network.png` 节点摆位不同（边完全相同） | 调用前 `set.seed()` |
| `ggrepel::geom_text_repel()` 用环境 RNG 做标签排布 | 12 张图里 11 张一致，只有 `volcano_plot.png` 不一致 | 传 `seed=` |

> **`gseGO(seed = 123)` 这个参数不足以保证复现。** 它被接受了（日志里没有退回警告），
> 但没有真正转发到 fgsea 的采样器。**不要依赖它** —— 直接在调用前设 RNG 状态。
> 初版就是踩了这个坑：以为传了 `seed=` 就完事了。
>
> **`ggrepel` 同理**：`seed` 的默认值是 `NA`（不是 `NULL`），不传就用当前 RNG 状态。
> 而 `save_pdf` 会把同一个绘图表达式求值两次（PDF 一次、PNG 一次），
> 第二次接着第一次消耗过的状态跑 —— 所以两次的标签位置本来就不同。

**第五个不是随机源，设种子也没用：**

| 来源 | 症状 | 处理 |
| --- | --- | --- |
| 多线程 BLAS 的归约顺序 | `deg_table.csv` 浮点末位分叉：`6.00193941779545e-05` vs `...546e-05`。t 统计量在 15 位有效数字上看着一样（实际不是），极小的 P 值把这个差异放大到可见 | 把 `OMP_NUM_THREADS` / `OPENBLAS_NUM_THREADS` / `MKL_NUM_THREADS` 钉为 1 |

> 这一类**不能靠设种子解决**，因为没有随机抽样可设 —— 只有线程调度。
> 它也让"之前几轮报逐字节一致"这个说法变得不严谨：比对是真的、结果确实一致，
> 但没有任何东西在固定算术，换台机器或换个负载就可能不一致。
> **结论比支撑它的机制更强**，这是要避免的。

另外 `enrichGO` 写出的 `geneID` 列**内部基因顺序不稳定**（实测 234 行里 230 行的
排列不同，而 p 值、计数完全一致）。这对科学结论没有影响，但让产物无法逐字节比对，
"两轮是否一致"这种验证就做不了。`normalise_gene_lists()` 在落盘前排序消除该差异。

**验证方式**：连续两轮 CI（`push` + `workflow_dispatch`），比对全部产物的 SHA256。
当前状态：**13 个结果表 + 12 张 PNG 全部逐字节一致**；
只有 `state.json` 不同，那是各步骤耗时与时间戳，属预期。

`assets/config.<GSE>.yml` 的 `analysis.seed`（默认 `20260919`）**不可删除**；
所有种子都紧挨着各自的随机调用设置，与上游消耗了多少随机数无关。

### 3.3 质控（`02_qc_pca_correlation.R`）

| 图 | 内容 | 判据 |
| --- | --- | --- |
| `boxplot_before_after.pdf` | 标准化前后各样本表达分布 | 标准化后中位数应齐平 |
| `density_plot.pdf` | 标准化前后密度曲线 | 曲线应重合 |
| `pca_plot.pdf` | PCA 散点（PC1/PC2），按分组着色 + 形状 | 记录 PC1/PC2 方差解释率；组内 95% 椭圆 |
| `pca_ellipse.csv` | 椭圆的坐标（每组 120 个点 + 半径系数） | 椭圆跨度 / 数据跨度比值记进日志 |
| `correlation_heatmap.pdf` | 样本间 Pearson 相关热图 | — |
| `correlation_matrix.csv` | Pearson + Spearman 相关矩阵 | **cor < 0.8 的样本标记为离群** |

额外计算**分组混杂指标**（§2.9）：`mean_within_group_pearson`、
`mean_between_group_pearson`、`group_confounded_with_global_shift`。
组间平均相关低于离群阈值时日志直接警告。

#### 3.3.1 组内置信椭圆：为什么半径是 2.45 SD，为什么坐标要落盘

**半径用 `sqrt(qchisq(0.95, 2))` = 2.45 SD，不用 ggplot 默认的口径。**
ggplot 的 `stat_ellipse` 默认 `type = "t"`，半径是 `sqrt(2 * qf(0.95, 2, n-2))`：

| n/组 | 半径系数 | 后果 |
| ---: | ---: | --- |
| 3 | **6.16 SD** | 椭圆是数据跨度的 6 倍，点被压成中心一小团 |
| 5 | 4.37 SD | 尚可 |
| 10 | 2.99 SD | 接近 `"norm"` |

本设计每组 3 个样本，所以用 `"norm"`。代价是它把协方差当作**已知**，
而 n=3 时协方差只由 **2 个自由度**估计 —— 椭圆**低估**了真实离散程度。
这不是一个严格的 95% 区间，副标题里写明了"a visual aid, not a test"。

**坐标落盘成 `pca_ellipse.csv`，因为"椭圆画没画"不能靠看图判断。**
实测用 `stat_ellipse()` 时它每组 3 个样本**产出空数据**：ggplot 不报错、
坐标范围也没被椭圆撑大，所以图上只有点、没有椭圆，而图注却写着有椭圆 ——
**"少画了一层"和"画对了"在图上看起来一样**。改成 `ellipse_points()` 自己算坐标后，
`pca_ellipse.csv` 有 240 行，日志给出跨度比值 2.39（椭圆 305.6 / 数据 127.8）。
点在面板里约占 42% 宽度 —— 这是 95% 椭圆应有的样子，但值得知道。

> **这是本仓库第三类"静默失败"**：第一类是空白图（设备开了又关），
> 第二类是被裁掉的副标题（§3.10），第三类是**少画的图层**。
> 三者共同点是**产物存在、大小正常、看起来没问题**。
> 所以判据只能是"算一遍"而不是"看一眼"。

### 3.4 差异表达（`03_deg.R`）

- 非配对：设计矩阵 `~ 0 + group`；配对：`~ 0 + group + patient`（§2.1）
- 对比 `tumor - normal`；配对设计若秩不足自动退回非配对并记录原因
- `limma::lmFit` → `contrasts.fit` → `eBayes(trend=TRUE, robust=TRUE)` →
  `topTable(n=Inf, adjust="BH")`（robust 不可用时退回标准 `eBayes`）
- 输出 `deg_table.csv`，**必须包含** `gene, logFC, AveExpr, t, P.Value, adj.P.Val, B`
- 显著标准：`adj.P.Val < 0.05` **且** `|log2FC| > thresholds.log2fc`（当前 0.5）
- 摘要写 `paired` / `residual_df` / `n_significant` / `n_nominal_only`，火山图标注 top 基因

> **把 `log2fc` 从 1.0 调到 0.5 不会让火山图出现显著基因。**
> 判据是 `adj.P < 0.05` **且** `|log2FC| > log2fc`，AND 关系，调小倍数阈值
> 只会让集合更小。GSE64790 的最小 `adj.P` 是 0.394，全基因组 BH 校正下
> **没有任何基因**能过第一关。限制在样本量（n=6），不在倍数阈值，见 §2.2。

**火山图用两个正交通道编码两个变量**（`03_deg.R`）：

| 通道 | 变量 | 取值 |
| --- | --- | --- |
| **颜色** | 方向 | up = 红 `#B2182B`、down = 蓝 `#2166AC`、不显著 = 灰 |
| **alpha + 大小** | 置信度 | FDR 显著 = 实心大点、名义显著 = 半透明小点、其余 = 极小淡点 |

**为什么不是把名义显著单独染成一种颜色。** 早先的做法是给这批基因一个独立的
橙色 —— 结果图上**看不到任何上调/下调**，因为红蓝两档是空的（FDR 无一通过）。
读者拿到一张只有橙点和灰点的火山图，读不出方向，而方向正是火山图存在的理由。

现在上调红、下调蓝一眼可见，而"未过 FDR"这个事实由**半透明**承载，
并写进图例（`nominal P < 0.05, NOT FDR-significant (n)`）与副标题。
两个通道编码两个变量，不是一色两义。GSE64790 实测 1909 个名义显著基因，
上调下调都有（像素核对：红 4865 / 蓝 4497）。

纵轴一律用 **raw P**。曾经想让名义显著的点用 raw P、其余用 `adj.P`，但那样
名义显著的点会画在 `-log10(0.05) = 1.3` 这条线**之上**，看起来像通过了 FDR ——
正好是这张图要避免的误导。而全部用 `adj.P` 也不行：GSE64790 的 `adj.P` 全在
0.394~1 之间，纵轴范围只有 0~0.4，整个图压成一张饼。
**FDR 阈值无法用一条横线表示** —— BH 校正是逐基因的，没有常数截断，这点写在副标题里。

#### 3.4.1 p 值分布诊断（DE 之后的 QC 关卡）

来自 K-Dense `bulk-rnaseq` skill 的 QC 清单：

> "A well-behaved test gives a roughly **uniform** histogram with a **peak near 0**
> (the true positives). A peak near 1, or a U-shape, signals a problem:
> misspecified design, unmodeled batch, or filtering issues.
> **Fix the design rather than trusting the gene list.**"

输出 `pvalue_histogram.pdf/.png`（直方图 + 均匀分布的期望线），并在
`deg_summary.json` 写入 `n_p_lt_0p001` / `frac_p_lt_0p05` / `pvalue_diagnosis`。

**这张图对本设计尤其关键**：n=6 时几乎不可能有基因通过 FDR，
光看"0 个显著基因"**分不清是功效不足还是模型设定错了**。
判据：P<0.05 的基因占比超过原假设期望（5%）2 倍以上 → 有真实信号；
低于 0.5 倍 → 基本没信号。

GSE64790 实测 **11.6%**（超出 2.3 倍），36 个基因 P<0.001 →
`signal_present_but_underpowered`，日志同时给出 WARN 说明这是样本量限制而非设计错误。

### 3.5 聚类热图（`04_heatmap_enrichment.R`）

- 取 top 50 显著 DEG（按 `adj.P.Val` 升序）
- **降级逻辑**：显著基因 < 50 时取全部显著基因；为 0 时取全表 top 20 并标记 `warning`
- 行聚类 euclidean + complete，表达量按行 Z-score
- 输出 `top50_heatmap.pdf`

### 3.6 preranked GSEA + GO / KEGG 富集（`04_heatmap_enrichment.R`）

富集走**两条路**，判据来自 K-Dense `pathway-enrichment` skill：

> "a discrete hit list → ORA; a ranked table with per-gene scores → GSEA"
> "Never threshold a list and then feed it to GSEA — that discards the ranking GSEA depends on."
> "Better when effects are broad/subtle or when a hit list would be very short or very long."
> "A very long one (> 2000) loses specificity — prefer GSEA in both extremes."

#### A. preranked GSEA（主力）

- **输入是完整的排序表，不卡任何阈值**。GSE64790 实测 13,948 个基因进入排序
  （16,487 个基因中能映射到 ENTREZ 的部分）
- **排序指标 = limma 的 moderated t 统计量**，不是 log2FC。理由同 skill：
  *"Rank by the test statistic (sign = direction, magnitude = evidence). This is
  more stable than ranking by log2FoldChange, which is noisy for low-count genes."*
- GO BP：`clusterProfiler::gseGO`；KEGG：`clusterProfiler::gseKEGG`（引擎均为 `fgsea`）
- 基因集大小限制 `minGSSize=15` / `maxGSSize=500`：过小的集合靠几个基因就能显著，
  过大的泛化集合（"metabolic process"）没有信息量
- `seed=123` 固定置换随机性，保证 p 值可复现
- **NES 的符号即方向**：`NES > 0` 表示该基因集在 tumor 一侧富集

#### B. ORA（辅助，**按上/下调分开跑**）

- 输入基因由 `select_degs()` 决定：FDR 显著基因，或降级到 raw P 前 N 个（§2.10）
- **上调、下调各自独立跑一遍 `enrichGO` / `enrichKEGG`**，结果表带 `direction` 列
- SYMBOL → ENTREZ（`clusterProfiler::bitr` + `org.Hs.eg.db`）
- `pAdjustMethod="BH"`，`pvalueCutoff=0.05`；KEGG 用 `organism="hsa"`
- dotplot 的 **x 轴是方向**，一眼看出条目由哪一侧驱动

> **为什么必须拆方向**（K-Dense `pathway-enrichment`）：
> *"ORA is direction-agnostic unless you split up/down lists; GSEA NES sign gives direction."*
> 对肿瘤 vs 正常组织这是致命的：上调的是增殖，下调的是基质/脂肪/血管。
> 混在一起跑会得到"两条方向相反的通路同时富集"这种无法解释的结果。
>
> **实测证明了这一点**：拆分前 `Integrin signaling`、`PI3K-Akt signaling pathway`
> 被报为"富集"；拆分后看清它们**全部来自下调基因**（属血管/基质簇）。
> "PI3K-Akt 下调"和"PI3K-Akt 上调"是完全不同的生物学陈述，合并分析无法区分。

#### C. 条目去冗余（两条路都做）

GO 会返回大量近义条目 —— 实测下调簇的前 4 名就是
`vasculature development` / `blood vessel morphogenesis` / `blood vessel development` /
`angiogenesis`，这是 **1 个发现重复了 4 次，不是 4 个发现**。

按**基因重叠 Jaccard 单链接聚类**折叠（`reduce_terms_by_overlap()`，阈值 0.5），
每类保留 `adj.P` 最小的那个作代表，并记录 `representative` / `cluster_size` 两列。
不引入 `GOSemSim` 这类重依赖。

实测效果：GSEA GO **1059 → 439** 个代表条目；ORA GO **234 → 72**。

#### D. 失败处理

**失败不终止**：KEGG REST API 有速率限制与授权限制，若返回空或报错，写入空表 +
`results/enrichment_status.json` 记录原因，流程继续。

> **背景集**：默认 `detected`（实测基因集，GSE64790 为 16,487 个）。
> K-Dense `pathway-enrichment` 把过大的背景列为 ORA 结果误导人的头号来源：
> *"Using too large a background makes ordinary housekeeping categories look
> significant — the most common way ORA results mislead."*
> 背景应当是"本实验**可能**检出的基因"，未在芯片上检出的基因不应计入。
> spec 里写的是"全基因组背景"；需要按 spec 口径复现时把
> `assets/config.<GSE>.yml` 的 `enrichment.universe` 改回 `genome` 即可。
> **两种结果不同，报告中必须写明用了哪种。**

### 3.7 差异基因互作（PPI，`05_ppi.R`）

- 输入基因由 `select_degs()` 决定：FDR 显著基因，或降级到 raw P 前 N 个（§2.10）
- symbol → `STRINGdb` 映射（`species=9606`）
- 置信度阈值 `score >= 400`
- **直接读 STRING 的 `protein.links` 文件建边**，不用 `STRINGdb::get_interactions()` ——
  后者实测在 GSE92252 上静默返回 0 行，导致"STRING 未返回任何达到阈值的互作"的假象
- 计算节点 degree，**hub 基因 = degree 排名前 10**
- 输出 `PPI_network.png`、`hub_genes.csv`、`ppi_edges.csv`
- **回退顺序**：STRING links 文件 → 共表达网络（标注 `method=coexpression`）→ 跳过
- **状态先落盘再画图**：`ppi_status.json` 在绘图**之前**写出，绘图单独 `tryCatch`。
  这样即使画图失败（如 `layout_with_fr` 拒绝负权重），节点/边数与原因仍然留存
- **回退不伪造网络图**：任何回退都在 `results/ppi_status.json` 标记 `method` 与原因

#### 3.7.1 网络图的可读性过滤（**只影响图，不影响数据**）

390 节点 / 3284 条边全画出来是一团毛线，读者只能看到"中间密、边上稀"。
过滤三步：

1. **最大连通分量**（零散小碎片对"互作网络"没有信息量）
2. **degree 前 `ppi_plot_max_nodes` 个节点**
3. **最强的 `ppi_plot_max_edges` 条边**，再丢掉因此孤立的节点

> **第 3 步才是密度控制，第 2 步不是。** 只按 degree 截节点是**反效果**：
> 实测取 degree 前 200 个节点后，边数反而从 3284 涨到 **4906** ——
> 高 degree 的节点彼此高度互联，取它们等于取网络最密的核。
> 密度是边的属性，不是节点的属性。

着色用 **Louvain 社区**（不是 degree 深浅 —— 连续渐变看不出任何结构）。
实测切出 20 个模块，而经计算验证的色盲安全色板只有 4 色，
所以**只给最大的 4 个上色**，其余归入灰色 `other`。
给 20 个模块上色必然走插值降级，插出来的颜色没验证过、彼此也分不开 ——
那样的图看着花花绿绿，实际读不出结构。完整网络仍在 `ppi_edges.csv`。

#### 3.7.2 布局：同心圆环，不是力导向

力导向（Fruchterman-Reingold）在这张图上是失败的：几百个节点挤成一团，
中间密到看不出结构、边上又空着，而且**每次运行布局都不一样**（一个随机源）。

同心圆环把"谁是 hub"直接编码成**半径**：内圈 = degree 最高的核心，
外圈 = 边缘基因。读者不用去找中心，一眼就知道层次。
实测 146 节点 / 450 边 → 3 环 **28 / 49 / 69** 个节点。

三条实现要点：

* **每环节点数按半径成比例**（周长 ∝ 半径），否则内圈挤成一坨、外圈稀稀拉拉
* **环内按社区排序**（同一模块占同一角度扇区）—— 这样模块内的边是短弦，
  模块间的边才跨圆心，边交叉大幅减少。按 degree 排会让每条边都横穿全图
* 半径等距递增，配一圈很淡的参考圆，让"几圈"这件事看得见

实测核对（`ppi_plot_layout.csv`，两轮 CI 逐字节一致）：

| 环 | 节点 | 半径 | degree 范围 | 角间距（理论） |
| --- | --- | --- | --- | --- |
| 1（核心） | 28 | 1.00 | **10–26** | 12.8–12.9°（12.9°） |
| 2 | 49 | 1.72 | **4–10** | 7.3–7.4°（7.3°） |
| 3 | 69 | 2.44 | **2–4** | 5.2–5.3°（5.2°） |

degree 区间严格不重叠，半径零误差，角度均匀。

> **图上的 hub 是子网络的 hub，不是全网络的 hub。** 环序按**过滤后**子网络的
> degree 排，而最强的 450 条边集中在一个增殖模块上 ——
> 核心环是 MAD2L1(26)、CENPA(24)、CCNB2(24)、MCM6(22)、MELK(22)、
> EXO1(20)、CHEK1(20)、BIRC5(20)。而 `hub_genes.csv` 排的是**全网络** 388 个节点，
> 前列是 GAPDH、CD34、IGF1、KIT、CXCL12。
>
> 两个排名都对，回答的是不同问题。**只给图不给这句话，读者会把核心环
> 当成"hub 基因"的答案，而那个文件给的是另一批基因。** 所以副标题里写明
> 环序来自哪个网络。

> **绘图必须与网络构建隔离。** 出图代码的错误如果逃逸，会被调用方的 `tryCatch`
> 当成"STRING 失败"接住，整条路径静默退化成共表达网络 ——
> **一个画图的 bug 换掉了分析方法**，而 `ppi_status.json` 里"有图、有 hub 基因"，
> 一切看着正常。实测连续两轮 CI 就是这么把共表达网络当成 PPI 结果发出去的。
> 根因是 STRING 分支建图时只传了两列，图里没有 `weight` 边属性。

> **布局也要落盘。** 光看 PNG 分不出"第 3 环是不是真的在外圈"，
> 从像素反推是猜。`ppi_plot_layout.csv` 写下每个节点的环号、半径、角度、
> 坐标、degree 与模块 —— 环结构就成了可核对的数据，而不是视觉印象。

### 3.8 配色：SCI 常用色板 + 算出来的判据

全部色值集中在 `scripts/lib/common.R` 的 `PAL` 与三个色板函数里，
**来源是 SCI 发表中通用的几套色板**，不是自己调的：

| 角色 | 来源 | 为什么是它 |
| --- | --- | --- |
| 上调 / tumor、下调 / normal | **ColorBrewer RdBu** 的两个端点 `#B2182B` / `#2166AC` | RdBu 是基因组学论文里最常用的发散尺；方向色与发散尺**同源**，热图上的红就是火山图上"上调"的那个红 |
| 显著性（连续） | **viridis**（截去暗端） | Nature Methods 推荐，感知均匀、灰度单调、色盲安全 |
| 模块 / 多组 | **Okabe-Ito** 及其加深变体，穷举筛出的 4 色 | 色盲安全定性色板的通行选择 |

判据取自 `better-colors` skill，其中一条给出了量化标准：

> **色相相差 15° 以内视为同一个颜色**，一个颜色只能有一个含义。

按这条算，本仓库先后暴露了四个**凭眼睛看不出来**的问题：

| 问题 | 实测数字 |
| --- | --- |
| p 值直方图的 `firebrick` 参考线 vs up 红 | 色相只差 **0.4°** —— 按判据就是同一个颜色 |
| 旧分类色板的棕 vs up 红（protanopia） | OKLab 距离 **0.002**，几乎完全重合 |
| magma 序列色中段 vs up 红 | 色相只差 **2.9°** —— 序列色的中段就是"上调红" |
| **viridis 暗端 `#365C8D` vs down 蓝** | 色相只差 **3.1°** —— 深色点会被读成"下调" |

最后一条是**门禁先发现、图还没画出来就拦下的**：完整 viridis 的暗端三档色相是
318 / 292 / **255**，而 down 蓝是 252.4。所以 viridis **必须截去暗端**，
从 `#277F8E`（h=211.3，距 down 41.1°）起用，保留感知均匀性与灰度单调性
（亮度 0.552 → 0.918）。

最终语义映射：

| 角色 | 色值 | 白底对比度 |
| --- | --- | --- |
| 上调 / tumor | `#B2182B` | 6.87:1 |
| 下调 / normal | `#2166AC` | 5.90:1 |
| 不显著 / other | `#BDBDBD` | 1.88:1（中性灰，不承载色相语义） |
| 显著性序列 | viridis 截断版（青 → 亮黄） | 亮度单调递增 |
| 模块分类 | `#2E6B3E` `#B5527A` `#56B4E9` `#E69F00` | 6.39 / 4.72 / 2.31 / 2.25 :1 |

模块四色是**穷举**出来的，不是贪心：贪心会选出两个绿色（h=149.9 / 165.5，
只差 15.6°）。约束三条 —— 与 up / down / 序列色色相拉开 15° 以上、
在三种色盲下 OKLab 距离 ≥0.05、白底对比度 ≥2.0。
最优组合的最小两两距离是 **0.246**。

`tools/check_palette.mjs` 从 `common.R` 解析出实际色值重新算一遍，接在 CI 里。
**这道门禁验过能失败**：把 viridis 的暗端加回去，它立刻报
`序列色板有一档 #365C8D 与 down 蓝 色相只差 3.1° < 15°，该档会被读成方向色`。

另外**颜色不是唯一的语义载体**（`better-colors` 的另一条）：
PCA 的分组同时用颜色和**形状**编码，密度图用颜色和**线型**，
火山图的显著性用 **alpha + 大小**（见 §3.4），黑白打印或色盲读者仍能分辨。

### 3.9 汇总（`main_analysis.R`）

按顺序编排 3.1→3.7，逐步 `tryCatch`，写 `results/state.json`：

```json
{"steps":[{"id":"deg","status":"ok","seconds":12.3},...],
 "required_failed":[],"optional_failed":["ppi"]}
```

仅当**必需步骤**（校验/清洗/QC/PCA/相关性/DEG）失败时以非零码退出。
一旦某必需步骤失败，**后续所有步骤都记为 `skipped`**，不再继续跑。

`check_acceptance()` 读取 `enrichment_status.json` / `ppi_status.json` 判断
"富集/PPI 是否有结果**或**有记录在案的原因"。注意这两个文件的形状不同：
前者是 `{"go": {...}, "kegg": {...}}`，后者把 `status` 放在顶层。

### 3.10 出图的三条硬约束（都是踩出来的）

**A. 副标题必须折行。** ggplot 的副标题**不换行**，超出画布宽度的部分被
**静默裁掉** —— 不是显示成省略号，是根本没有。所以"字没显示全"从图上
完全看不出来：一张被切掉尾巴的图和一张本来就短的图长得一样。
实测 PPI 副标题 455 字符、火山图 275 字符，尾巴都进了裁纸刀。

`wrap_subtitle()` 按画布宽度和字号算每行字符数（副标题 8.5pt、无衬线体
平均字宽约 0.52 em，留 15% 余量），全部副标题都走它。
`tools/check_r_syntax.mjs` 挡住超过 100 字符又没折行的字面量 —— 门禁本身
也验过会响（临时塞一个 136 字符的副标题，它报错）。

> 门禁第一版有假阳性：`plot.subtitle = element_text(...)` 是主题设置，
> 它把后面的 `"bottom"` / `"horizontal"` 当成副标题文本算了进去。
> 现在排除 `element_text` 和带 `.` 前缀的匹配。

**B. 图例默认在底部横排。** 原来全部在右侧，直接吃掉图的宽度。
PPI 那张图本来就是方的（`coord_fixed`），却为一个右侧图例多留了 2.5 英寸。
改到 `theme_paper()` 默认底部横排后，12 张图的面积降到 69%–94%：

| 图 | 之前 | 之后 | 面积 |
| --- | --- | --- | ---: |
| PPI_network | 1500×1275 | 1125×1170 | 69% |
| boxplot_before_after | 1500×900 | 1200×825 | 73% |
| density_plot | 1200×1050 | 1050×900 | 75% |
| GO / KEGG dotplot | 1350×1200 | 1200×1050 | 78% |
| top50_heatmap | 1200×900 | 1050×825 | 80% |
| GSEA dotplot | 1350×1050 | 1200×975 | 83% |
| correlation_heatmap | 1200×1050 | 1050×975 | 81% |
| pca / pvalue_histogram | 1050×900 | 975×900 | 93% |
| volcano_plot | 1275×975 | 1200×975 | 94% |

pheatmap 的色条固定在右侧、没有位置参数，它是细长条、占宽有限，保留。

**C. 图例键的填充色要显式给。** PPI 节点的 `degree` 图例是**空的**：
节点是 `shape = 21` + `colour = "white"`，而 `fill` 映射在另一个 scale 上，
size 图例的键继承不到 `fill`，于是画成"白描边 + 无填充"——
**在白底上完全隐形**，图例只剩一个标题。用
`override.aes = list(fill = PAL$muted, ...)` 补回来。
实测修复后底部区域有 1885 个 `#666666` 像素、5 个横向分离的色块，
宽度 14/19/23/25/28 px —— 正好对应 degree 5/10/15/20/25 五个递增的圆。

`pretty()` 还会给出数据范围外的断点（degree 从 2 起却标了 0），一并过滤。

---

## 4. 产出清单

**每张图都同时出 PDF 和 PNG**：PDF 是矢量图（放大不失真），PNG 是为了能直接看 ——
artifact 是打包成 zip 下载的，PDF 在里面不能预览，得先解压再找阅读器。
PNG 由 `save_pdf()` 在 150 dpi 下渲染，PNG 失败只记 warning，不影响 PDF 与流程。

```
results/
├── boxplot_before_after.pdf/.png  QC：标准化前后箱线图
├── density_plot.pdf/.png          QC：密度曲线
├── pca_plot.pdf/.png              PCA（含方差解释率）
├── correlation_heatmap.pdf/.png   样本相关性热图
├── correlation_matrix.csv         Pearson + Spearman 矩阵 + 离群标记
├── deg_table.csv                  全基因差异分析表
├── volcano_plot.pdf/.png          火山图
├── pvalue_histogram.pdf/.png      DE 后 QC：p 值分布（均匀性 + 0 附近是否有峰）
├── top50_heatmap.pdf/.png         top DEG 聚类热图（Z-score）
├── GSEA_GO_dotplot.pdf/.png + GSEA_GO_table.csv     preranked GSEA / GO BP（主力）
├── GSEA_KEGG_dotplot.pdf/.png + GSEA_KEGG_table.csv preranked GSEA / KEGG
├── GO_dotplot.pdf/.png + GO_table.csv       ORA GO BP（含 direction 列）
├── KEGG_dotplot.pdf/.png + KEGG_table.csv   ORA KEGG（含 direction 列）
├── PPI_network.png + hub_genes.csv + ppi_edges.csv   STRING PPI 与 hub 基因
├── enrichment_status.json         富集模式（fdr / ranked_fallback）、GSEA 参数、去冗余阈值及原因
├── ppi_status.json                PPI 方法、节点边数及回退原因
└── state.json                     各步骤执行状态 + 验收结果
```

---

## 5. 已知局限（必须随结果一并报告）

1. **样本量**：n=6（3 vs 3 配对，残差 df = 2）。**没有任何基因能通过 FDR**
   （最小 `adj.P` = 0.394），下游富集与 PPI 走的是 `ranked_fallback` 降级路径（§2.10）。
   所有结果只能作为**假设生成**，不能作为临床或机制结论。
   **本设计的 MDE 为 d_z = 3.26（80% 功效）**，见 §2.2.1 —— 中等效应基本检不出。
2. **组织成分混杂**：肿瘤 vs 全组织正常，差异主要反映细胞组成而非肿瘤特异性表达。
   本数据集尤为明显 —— top 基因是 KRT14、SPARCL1、TAGLN、SDPR、PPARG 等
   基质/脂肪/上皮比例相关基因。GSEA 的下调簇（血管发育、肌肉系统、循环系统）
   与上调簇（有丝分裂、染色体分离）正是这一组成的两个侧面。
3. **配对自由度极低**：3 对只剩 2 个残差自由度。实测配对并没有换来功效优势
   （MDE 3.26 vs 不配对的 3.07，§2.2.1）—— df 从 4 降到 2 的代价超过了方差缩减的收益。
   若能拿到更多配对样本，功效会显著改善：d_z=1.5 需 6 对，d_z=1.0 需 10 对。
4. **无独立验证队列**：本设计不含验证集，结果未经任何外部数据复现。
5. **lncRNA 芯片的平台局限**：GPL19612 的 65,531 个探针里只有 33.3% 带 `GeneSymbol`，
   覆盖 16,487 个基因。非编码转录本的信息在 symbol 折叠时被丢弃了。
   另有 2,539 个基因无法映射到 ENTREZ，未进入 GSEA 排序表（13,948 个）。
6. **单平台单批次**：无法评估批次效应，也无法做跨平台一致性检验。
7. **富集分析**：
   - KEGG 依赖在线 API，可能因限流返回空结果；此时结论中不得声称"无 KEGG 通路富集"，
     只能声称"本次未获得 KEGG 结果"。
   - **GSEA 与 ORA 的结论强度不同**：GSEA 不卡阈值、用完整排序表，对弥散效应敏感；
     ORA 的输入是 `ranked_fallback` 的前 500 个基因，**只能表述为
     "在最显著的 N 个基因里富集到……"**，不能表述为"显著差异基因富集到……"。
   - 两条路都做了基因重叠去冗余（Jaccard ≥ 0.5 折叠为一类），
     报告时应引用 `representative` 列，不要罗列同一簇的多个近义条目。

---

## 6. 队列级设计：GSE42568

上面 §1–§5 是 `design_mode: small_sample` 的设计（GSE64790，n=6）。
GSE42568 是**同一个对比的放大版**，用途是把 GSE64790 上"信号存在但功效不足"
的假设放到有功效的队列上验证。

### 6.1 合规性核验

| 约束 | 要求（`cohort`） | GSE42568 实测 | 结论 |
| --- | --- | --- | --- |
| 物种 | Homo sapiens | `Homo sapiens` (taxid 9606) | ✅ |
| 数据类型 | 基因芯片 | `Expression profiling by array` / GPL570 | ✅ |
| 疾病 | 乳腺癌 | 104 例癌 vs 17 例正常乳腺 | ✅ |
| 样本量 | ≥ 15 | **121** | ✅ |
| 每组样本数 | ≥ 10 | 104 vs 17 | ✅ |
| 基因注释 | 需要 symbol | GPL570 有 `Gene Symbol` 列 → 23,500 个基因 | ✅ |

### 6.2 实际结果（run 35434644032）

**DEG 是实打实的，没有走降级路径**：`deg_mode: "fdr"`，3,795 个基因通过
`adj.P < 0.05` 且 `|log2FC| > 1`（上调 1,845 / 下调 1,950），残差自由度 119。
这与 GSE64790 的 0 个形成直接对照 —— **同一套代码、同一个对比，差别只在样本量**。

| | GSE64790 (n=6) | GSE42568 (n=121) |
|---|---|---|
| FDR 显著基因 | **0** | **3,795** |
| 最小 adj.P | 0.394 | — |
| 富集输入 | `ranked_fallback` 前 500 | `fdr` 3,795 |
| GO 条目（去冗余后） | — | 799 → 278 |
| PPI | 共表达回退（STRING 映射 485/500） | STRING 3,254 节点 / 92,872 边 |

### 6.3 临床终点与 EPV 上限

`tools/check_clinical_endpoints.mjs` 实测（**这是选数据集时最容易漏掉的一步**）：

| 终点 | 事件数 / 风险集 | Cox 签名上限（EPV ≥ 10） |
| --- | --- | --- |
| `overall survival event`（配 `overall survival time_days`） | **35 / 104** | **3 个基因** |
| `relapse free survival event`（配 `relapse free survival time_days`） | **48 / 104** | 4 个基因 |
| `er_status` | 67 vs 34（少数类 34） | 3 个基因（logistic） |
| `lymph node status` | 59 vs 45（少数类 45） | 4 个基因（logistic） |

**35 个事件意味着 OS 模型超过 3 个基因就开始过拟合。** 文献里常见的
"8 基因预后签名"配 35 个事件是过拟合，只是很少被查。

### 6.4 配套验证队列：GSE20685

| | GSE42568 | GSE20685 |
|---|---|---|
| 平台 | **GPL570** | **GPL570** |
| n | 121（104 癌 + 17 正常） | 327（全为原发癌） |
| 终点 | OS 35 事件 / RFS 48 事件 | 死亡 83 事件 / 转移 83 事件 |
| 其他临床 | ER / grade / 淋巴结 / 年龄 / 大小 | 分期 / 亚型(I–VI) / 化疗方案 / 随访年数 |
| 角色 | 发现集 | **独立验证集** |

**两者同平台**，注释和表达量尺度直接可比，省掉跨平台校正这一整块风险。

> ⚠️ GSE20685 的事件列与时间列**名称对不上**（`event_death` 配
> `follow_up_duration (years)`），工具会如实报"需人工确认配对关系"而不是猜。

### 6.5 WGCNA 与 LASSO-Cox（步骤 06 / 07，已实现）

两个都是**可选步骤**（`STEPS` 里 `required = FALSE`）：不适用时写状态文件说明原因，
而不是静默跳过。**"失败"和"不适用"是两回事** —— 它们都让 job 保持绿色，
所以验收项查的是状态文件里的 `status` 字段（`settled()`），
文件不存在或字段缺失一律记 FAIL。实测第一次跑时 WGCNA 崩了而 CI 全绿，
就是这个验收项把它揪出来的。

#### 步骤 06：WGCNA 共表达网络

- **只用 104 例癌，不带那 17 例正常。** 带上正常样本的话第一个模块必然是
  "肿瘤 vs 正常"轴，反映的是组织成分 —— 而 DEG 已经答过那件事了。
  WGCNA 要回答的是癌组织**内部**的异质性。
- 门禁：`design_mode == "cohort"`、n ≥ 15、基因数 ≥ 2000。
  **GSE64790（n=6）必然跳过**，状态是 `not_applicable` 并写明原因。
- 输入取方差最大的 5000 个基因；`goodSamplesGenes` 过滤。
- 软阈值：取**最小**的 power 使 scale-free R² ≥ 0.8；达不到就记录最大 R²
  并标 `reached: false`，**不假装通过**。实测 GSE42568 得 power=8（R²=0.858）。
- `blockwiseModules(networkType="signed", TOMType="signed", mergeCutHeight=0.25,
  randomSeed=<seed>)`。`randomSeed` 必须显式传（规则 11）。
  **`corFnc` 传函数对象没用** —— 见 AGENTS.md 规则 23，包内按字符串查函数，
  必须临时挂载 WGCNA。
- 模块-性状关联：**Pearson + BH 校正**。几十上百次检验不校正的话，
  p<0.05 的格子会有一堆是偶然的，而它们在热图上和真信号长得一样。
- 产物：`wgcna_modules.csv`、`wgcna_module_sizes.csv`、`wgcna_soft_power.csv`、
  `wgcna_module_trait.csv`、`wgcna_soft_power.pdf`、
  `wgcna_module_trait_heatmap.pdf`、`wgcna_status.json`。

#### 步骤 07：LASSO-Cox 预后签名

- **终点必须由 config 显式指定**（`survival.time_column` / `event_column` /
  `event_value`），不做关键词自动配对。实测 GSE20685 是 `event_death` 配
  `follow_up_duration (years)` —— 名字里没有共同词，自动配对一定会配错，
  而**配错不报错**，只会算出错的 C-index。
- EPV 门禁（规则 7）：GSE42568 的 OS 只有 35 个事件 → 上限 3 个基因。
- 候选基因 = FDR 显著 DEG；没有显著基因时退到 `ranked_fallback_genes`
  并按原始 P 排序，**状态里标明是降级**。
- `glmnet(family="cox", alpha=1)` + 重复 CV（显式 `foldid`）+ `lambda.1se`。
- **C-index 用自己实现的 `harrell_c()`，不用 `survival::concordance()` 的公式接口** ——
  实测后者在 `Surv(time, event) ~ risk` 下返回 1 - Harrell C（训练集 0.121
  vs 交叉验证 0.793，正好互补），而 0.121 看着像个正常数字，不会引起怀疑。
- 报**三个** C-index：训练集、交叉验证、外部验证（GSE20685）。
  验证集用**训练集的均值方差**标准化，KM 切点用**训练集的中位数** ——
  在验证集里重新取中位数会让两个队列的"高风险"不是同一个定义。
- 报重复 CV 选出基因数的**分布**（实测 [16, 3, 3, 22, 3]，只换 foldid 就差 7 倍）
  和 `signature_stable`。只报"最终签名 N 个基因"是把不稳定性藏起来。
- 当 `lambda.1se` 超出 EPV 上限时，**额外给一个显式满足 EPV 的版本**
  （`lasso_coefficients_epv.csv`，取"非零系数 ≤ epv_cap"的最大 lambda），
  两个都在外部队列上打分，让读者看到少要基因的代价。
- 产物：`lasso_coefficients.csv`、`lasso_coefficients_epv.csv`、
  `lasso_risk_scores.csv`、`lasso_stability.csv`、
  `lasso_selection_frequency.csv`、`lasso_cv_curve.csv`、`lasso_km.pdf`、
  `lasso_status.json`。

> **GSE64790 做不了 LASSO**：没有随访终点，且 `min(n, p) = 6`、
> 最大事件数 6 → EPV 上限 0。状态是 `not_applicable`。

### 6.6 尚未做的（本设计范围之外）

---

## 7. 复现

```bash
# GitHub Actions：指定数据集（push 时两个数据集各跑一个 job）
gh workflow run geo_analysis.yml -f dataset=GSE42568
gh run download --name geo-results-GSE42568

# 本地（需 R 4.3+ 与 Bioconductor）。没有默认配置，必须显式指定。
Rscript scripts/main_analysis.R --config assets/config.GSE42568.yml
```

**已验证的成功运行**：[run 35434644032](https://github.com/liubarryteb12/geo-normal-pipeline-skill/actions/runs/35434644032)
—— 两个 job 全绿（GSE64790 5m47s、GSE42568 6m25s），7 个步骤全部 `ok`。

> **关于"逐字节可复现"：不要把浮点末位也当成已验证的性质。**
>
> 同一 commit 连跑多轮会得到**两组**不同的浮点末位。实测数据（GSE64790）：
>
> | 轮次 | Azure 区域 | 数据文件指纹 |
> |---|---|---|
> | c370bd0 ×3 | westus3 / eastus2 / westus2 | `04d8295b` |
> | c370bd0 ×2 | eastus / eastus | `04d8295b` |
> | **c370bd0 ×1** | **westus3** | **`6839eb19`** |
> | 3f2a373 ×2 | centralus | `6839eb19` |
>
> 关键那一行：**同一个 commit、同一个区域（westus3）跑出了另一组末位。**
> 所以差异既不是 commit 引起的（我曾据此误判过一次），也不是区域引起的 ——
> 是 GitHub 托管 runner 的 CPU 型号在**同一 Azure 区域内也不单一**。
>
> 差异幅度：`deg_table.csv` 最大绝对差 **9.9e-14**（logFC 量级 ~1）、
> `pca_ellipse.csv` 最大 **1.0e-12**（坐标量级 ~100），都是**约 1 个 ULP**。
> 而同一批运行的 `data/` 输入文件（`group.csv` / `clean_stats.json` /
> `feature_mode.json`）**逐字节一致**，所有**报告出来的**统计量也完全一致
> （16487 个基因、PC1=46.1%、1909 个名义显著、GSEA 1080/101、
> STRING 485/500、17 项验收全过）。
>
> 机制是 OpenBLAS 在**运行期**按 CPU 型号分发 SIMD 内核，向量宽度不同则
> 归约顺序不同。`OMP_NUM_THREADS=1` 只解决了**线程调度**那一半，
> 管不到**内核分发**这一半 —— 所以 workflow 里另外钉了 `OPENBLAS_CORETYPE`。
>
> **可以放心声称的是结构与量级可复现；不要把"浮点末位一致"写进结论。**
> 验证方式是同一 commit 连跑多轮比对 SHA256（AGENTS.md 规则 12）。

完整配置见 `assets/config.<GSE>.yml`；运行期故障排查见 `references/troubleshooting.md`。
