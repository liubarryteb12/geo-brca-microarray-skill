---
name: geo-brca-microarray-skill
description: Run an end-to-end GEO gene-expression microarray data-mining pipeline for human breast cancer with fewer than 10 samples - download and clean, QC, PCA, sample correlation, limma differential expression, clustered heatmap, GO/KEGG enrichment, and STRING PPI - and ship it as a GitHub Actions workflow that uploads results as an artifact. Use when the user asks for GEO data mining, a GEO/GSE microarray analysis, a breast cancer expression analysis, a small-sample limma DEG workflow, GO/KEGG enrichment on GEO data, a GitHub Actions bioinformatics pipeline, or asks to reproduce or extend GSE64790.
license: MIT
compatibility: R 4.3+ with Bioconductor (GEOquery, limma, clusterProfiler, org.Hs.eg.db, enrichplot, STRINGdb, pheatmap, impute, igraph). Node 18+ for the dataset pre-flight tool. Needs outbound HTTPS to NCBI GEO, STRINGdb and KEGG. Designed to fit GitHub Actions ubuntu-latest in under 20 minutes.
metadata:
  version: "1.0"
  skill-author: built for DSH from a user-supplied agent_spec
---

# /geo-brca-microarray-skill

GEO 乳腺癌小样本芯片数据挖掘流水线。默认数据集 **GSE64790**（人源、Agilent lncRNA 芯片
GPL19612、6 例 = 3 TNBC vs 3 配对正常乳腺组织）。

完整实验设计见 [`EXPERIMENTAL_DESIGN.md`](EXPERIMENTAL_DESIGN.md) —— **改任何东西之前先读它**，
尤其是 §2.2 关于 n=6 的统计功效边界（**没有任何基因能通过 FDR**）和 §2.10 的降级路径。

## 触发场景

- 「用 GEO 数据做乳腺癌差异表达分析」「GSE 芯片数据挖掘」
- 「小样本（<10）芯片怎么做 limma + GO/KEGG」
- 「搭一个 GitHub Actions 生信流水线，跑完上传 artifact」
- 「复现/扩展 GSE64790」
- 任何要求「清洗 + QC + PCA + 相关性 + 差异基因 + 热图 + 富集 + PPI」全流程的任务

## 工作流

### 1. 先验证数据集，不要相信标题

```bash
node scripts/find_dataset.mjs check GSE64790
node scripts/find_dataset.mjs search --disease "breast cancer" --max-samples 10
```

这一步不是可选的。GEO 的标题与二手描述经常与真实元数据不符 —— 本项目最初提出的
GSE197894 被描述为"表达谱芯片"，实际 `gdstype` 是 RNA-seq，样本 20 例，
照它做流水线会在第一步就终止。**类型必须以 `gdstype` / `GPL` 为准。**

然后查样本相关结构，确认分组没有被批次效应混杂：

```bash
node tools/check_sample_structure.mjs GSE64790
```

这一步同样是"先花 10 秒，省掉 15 分钟"：等 R 流水线跑完才发现分组与一个全局表达位移
共线，DEG 结果就已经不可用了。第一版用的 GSE92252 正是在这里被查出问题并换掉的
（见 Gotchas 第一条）。

**还要查平台注释列有没有基因 symbol。** 形式条件全过的数据集可能根本没有基因注释
（GSE112848 的平台只有 5 列，无任何 symbol/accession），那样富集和 PPI 会全部跳过。

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

- **先查分组是否与全局表达位移混杂，再解读 DEG。** GSE92252 实测：9 个样本分成
  三个紧致簇，簇内 Pearson r > 0.92，而 **HER2− ↔ HER2+ 只有 0.31–0.38、
  HER2− ↔ 正常只有 0.11–0.19**，三个簇恰好就是三个分组。这个结构在**原始沉积值**
  里就存在（原始强度簇内 0.986 / 簇间 0.602），不是 log2 或 quantile 引入的。
  后果：tumor-vs-normal 的差异分不清多少来自恶性转化、多少来自这个位移。
  `results/qc_summary.json` 的 `group_confounded_with_global_shift` 为 `true` 时，
  结论只能按假设生成写，不能按肿瘤特异事件写。**换数据集时用
  `tools/check_sample_structure.mjs` 提前拦掉。** 换到 GSE64790 后这三项是
  0.9082 / 0.894 / **false**。
- **`min Pearson < 阈值 即离群` 这条规则会整体失效。** 它假设所有样本是同一组织的
  技术重复；一旦分组自带全局位移，**每个**样本都会与另一组的样本低相关，
  于是报出"9 个样本全是离群"。这时要看 `mean_within_group_pearson` vs
  `mean_between_group_pearson`，而不是离群样本个数。
- **注释覆盖率低不等于不能做基因层面分析。** lncRNA / 外显子芯片上大部分探针
  **本来就是非编码的**，覆盖率天然上不去：GPL19612 只有 33.3% 的探针带 `GeneSymbol`，
  但那是 **16,487 个基因**，做 GO/KEGG/STRING 绰绰有余。判据要**覆盖率与基因数取「或」**
  （`< 50%` **且** `< 5000 个`才退回探针模式），只卡覆盖率会把这类平台误判成探针模式，
  下游全部跳过。
- **不要用 `STRINGdb::get_interactions()`。** 实测 v12.0 上 1321/1423 个基因映射成功、
  links 文件也下好了，但它**返回 0 行且不报错**，流程被静默推进到共表达回退。
  本流水线直接读 `protein.links` 文件自己取子图。
- **`igraph::layout_with_fr()` 不接受负权重。** Spearman r 可以为负，直接把带符号的 r
  当权重会报 `Weights must be positive for Fruchterman-Reingold layout` 并中断步骤。
  布局权重取 `|r|`，有符号的 r 另存在边表里。
- **先落状态 JSON，再出图。** 出图是最后一步，画不出来时如果状态还没写，
  就会同时丢掉图和状态文件（本仓库踩过：`PPI_network.png` 和 `ppi_status.json`
  一起消失，只剩 `hub_genes.csv`）。
- **验收检查要能处理两种状态文件形状。** `enrichment_status.json` 是
  `{"go": {...}, "kegg": {...}}`，而 `ppi_status.json` 把 `status` 放在**顶层**。
  无条件写 `s[[key]]$status` 会在后者上抛 `$ operator is invalid for atomic vectors`，
  而且是在**所有分析步骤都成功之后**才崩，看起来像流程失败。
- **`<<-` 在 `tryCatch` 表达式里会跳过当前帧。** 表达式在调用函数的帧里求值，
  `<<-` 从**外层**环境开始找，于是本地变量没被赋值、外层环境被悄悄创建了一个同名变量。
  本仓库踩过：配对设计实际生效（残差 df = 2），摘要却报 `paired: false`。
  要在 `tryCatch` 里改本地变量，就把赋值挪到 `tryCatch` 外面。
- **`force(expr)` 两次不等于跑两遍。** R 的 promise 有记忆，第二次 force 直接返回
  缓存值。想用同一段绘图代码出两种格式，必须 `substitute()` 抓住**未求值**的表达式
  再 `eval()` 两次。本仓库踩过：`save_pdf()` 里 `force(expr); force(expr)`
  让 PNG 设备开了又关、什么都没画，产出 8 张空白 PNG，而 PDF 正常、日志无警告、
  验收全过 —— 因为验收只判断文件**存在**。所以 CI 里另有一道
  `node tools/check_figures.mjs results`，独立解码 PNG 像素判空。
- **每张图同时出 PDF 和 PNG。** PDF 是矢量图，PNG 是为了能直接看 ——
  artifact 是打包成 zip 的，PDF 在里面不能预览。PNG 失败只记 warning，不中断流程。
- **`GES` 不是真实的 GEO 前缀。** 正确的是 `GSE`（GEO Series）。用户说 GES 时按 GSE 处理。
- **分组模式不要凭直觉写。** GSE92252 的 `tissue` 字段里，HER2− 肿瘤写的是
  `... HER2-negative breast tumor`，而 HER2+ 肿瘤写的是 `... HER2-positive tumor`
  —— **没有 "breast"**。用 `"breast tumor"` 作模式会漏掉 3 个样本，第 00 步直接报
  「分组失败」。判别子串要用 `"tumor"`。GSE64790 则是 `TNBC  tissue`（**两个空格**）
  与 `matched normal breast tissues`。改任何数据集前先跑
  `node scripts/find_dataset.mjs samples GSEXXXXX` 看真实取值。
- **样本量 < 10 且要两组对比，几乎排除了所有组织样本数据集。** 人源乳腺癌芯片里
  n<10 的绝大多数是细胞系加药实验，没有"肿瘤 vs 正常"两组。约 200 个候选里只有
  4 个是真正的组织两组设计。若用户给的数据集是细胞系，分组字段应改为 treatment/control，
  而不是硬套 tumor/normal。
- **20 分钟不是宽裕的预算。** `clusterProfiler` + `STRINGdb` + `enrichplot` +
  `org.Hs.eg.db` 从源码编译要 25–40 分钟。workflow 用
  `r-lib/actions/setup-r-dependencies@v2` 走 Posit Package Manager 二进制包来压时间。
  如果换了包或 R 版本导致退化为源码安装，`timeout-minutes: 20` 会直接杀掉 job。
- **`enrichKEGG` 依赖 KEGG 在线 REST API**，会因限流或授权返回空。按 spec 这不算失败，
  但**报告中不得写"无 KEGG 通路富集"**，只能写"本次未获得 KEGG 结果"，并引用
  `results/enrichment_status.json` 里的原因。
- **有完整排序表时用 GSEA，不要卡阈值跑 ORA。** 这是本仓库最容易犯的方法学错误，
  初版就犯了。K-Dense `pathway-enrichment`：*"a discrete hit list → ORA;
  a ranked table with per-gene scores → GSEA"*、*"Never threshold a list and then
  feed it to GSEA"*。GSE64790 实测差距是**两个数量级**：GSEA 在完整排序表上拿到
  **1,059 条**显著 GO BP 条目（最好 adj.P = 1.0e-8），而 ORA 在 top-500 上最好只有
  3.1e-6。**弱功效、效应弥散的数据正是 GSEA 被设计出来处理的场景。**
  排序指标用 limma 的 moderated `t`，不要用 log2FC（低表达基因的 logFC 噪声极大）。
- **ORA 必须按上/下调分开跑。** ORA 本身方向无关，混在一起跑就分不清某条通路是被
  上调还是下调基因驱动的。本仓库实测：合并时 `Integrin signaling`、`PI3K-Akt
  signaling pathway` 被报为"富集"；拆开后看清它们**全部来自下调基因**（血管/基质簇）。
  "PI3K-Akt 上调"和"PI3K-Akt 下调"是完全不同的生物学陈述。
- **背景集默认用实测基因集（`detected`），不要用全基因组。** K-Dense
  `pathway-enrichment` 把过大的背景列为 ORA 误导人的头号来源：*"Using too large a
  background makes ordinary housekeeping categories look significant."*
  背景应当是"本实验**可能**检出的基因"。spec 写的是全基因组，改
  `enrichment.universe: genome` 可切回。
- **GO 条目必须去冗余后再报告。** 实测下调簇的前 4 名是
  `vasculature development` / `blood vessel morphogenesis` / `blood vessel development` /
  `angiogenesis` —— 这是 1 个发现重复 4 次。按基因重叠 Jaccard 单链接聚类折叠
  （`reduce_terms_by_overlap()`），GO 从 1059 折到 439。报告时引用 `representative` 列。
- **`n < 10` 时一定要看 p 值直方图。** K-Dense `bulk-rnaseq` 的 QC 关卡：
  分布应接近均匀且在 0 附近有峰。**光看"0 个显著基因"分不清是功效不足还是模型设定错了**，
  这张图能分开。GSE64790 实测 P<0.05 占 11.6%（期望 5%）、36 个 P<0.001 →
  `signal_present_but_underpowered`。峰在 1 或 U 形说明设计有问题，那时候连排序表都不能用。
- **不要报 post-hoc observed power。** 用观测到的效应反算功效是循环论证
  （它是 p 值的确定性函数）。要报就报**敏感性分析**：固定 n 下的 MDE。
  GSE64790 的 MDE 是 d_z = 3.26（80% 功效）—— 中等效应根本检不出。
- **配对在 n=3 时不换来功效。** 实测配对 MDE 3.26 vs 不配对 3.07：df 从 4 降到 2、
  t 临界值从 2.78 升到 4.30，代价超过了方差缩减的收益。配对仍然**正确**（控制个体基线），
  但别指望 3 对配对能提升检出能力。d_z=1.5 需要 6 对，d_z=1.0 需要 10 对。
- **有三个随机源必须先设种子，否则整条流水线不可复现。** 实测连续两轮 CI 的 GSEA
  显著 GO 条目数不同（1059 / 1110），看着像蒙特卡洛波动，其实源都不在 GSEA 本身：
  1. `impute.knn` 用 `sample()` 处理并列近邻 → 插补值变 → t 统计量变 → GSEA 排序变 → 全变
  2. `fgsea` 的 `fgseaMultilevel` 自适应采样 → NES / p 值第 4 位有效数字分叉
  3. `igraph::cluster_louvain()` 的初始顺序 → 模块划分变 → PPI 图上色变
     （**力导向布局 `layout_with_fr` 也曾经是第四个源，换成同心圆环布局后就没了**）

  **`gseGO(seed = 123)` 这个参数不足以保证复现** —— 它被接受（日志无退回警告）
  但没转发到 fgsea 的采样器，别依赖它，直接在调用前 `set.seed()`。
  种子要紧挨着各自的随机调用设置。验证方式：连跑两轮比对 SHA256，
  当前 14 个结果表 + 12 张 PNG 全部逐字节一致。
- **第四个随机源是 `ggrepel`，而且它连"同一张图的 PDF 和 PNG 不一致"都会造成。**
  `geom_text_repel` 的 `seed` 默认是 `NA`（**不是 `NULL`**），不传就用环境 RNG
  做标签排布的模拟退火。而 `save_pdf` 会把同一个绘图表达式求值两次（PDF 一次、
  PNG 一次），第二次接着第一次消耗过的状态跑 —— 所以两次的标签位置本来就不同，
  跨运行更会随上游随机数消耗量漂移。实测 12 张图里 11 张一致，只有火山图不一致。
- **设种子解决不了多线程 BLAS。** 实测 `deg_table.csv` 出现
  `6.00193941779545e-05` vs `6.00193941779546e-05`：t 统计量在 15 位有效数字上
  看着一样（实际不是），残差平方和走并行归约，线程调度不同就末位不同，
  极小的 P 值把这个差异放大到可见。**只能把 `OMP_NUM_THREADS` /
  `OPENBLAS_NUM_THREADS` / `MKL_NUM_THREADS` 钉为 1。**
  这条的教训是：报"逐字节一致"之前先确认**是什么机制在保证它** ——
  之前几轮报的一致有一半是运气。
- **画图代码的错误会被当成"分析方法失败"。** 实测踩过：STRING 分支建图时只传了
  两列，图里没有 `weight` 边属性；出图代码按权重筛边拿到 `NULL`，报
  `argument 1 is not a vector`；这个错误逃逸到 STRING 分支的 `tryCatch`，
  被当成"STRINGdb 失败"，整条路径静默退化成共表达网络。
  **一个画图 bug 换掉了分析方法**，而 `ppi_status.json` 里"有图、有 hub 基因"，
  一切看着正常 —— 连续两轮 CI 把共表达网络当 PPI 结果发了出去。
  绘图必须自己 `tryCatch` 兜住，方法本身如实记录。
- **网络图按节点筛是反效果，密度要按边筛。** 388 节点 / 3284 边全画是毛线球；
  但只取 degree 前 200 个节点后边数反而涨到 **4906** —— 高 degree 的节点彼此高度互联，
  取它们等于取最密的核。正确顺序：最大连通分量 → degree 前 N 个节点 → **最强的 M 条边**。
  另外按 degree 深浅上色看不出任何结构，要按 **Louvain 社区**上色。
- **力导向布局在几百个节点上就是一团糊，用同心圆环。** 半径直接编码 degree：
  内圈 = hub 核心，外圈 = 边缘基因，读者不用去找中心。三条要点 ——
  每环节点数**按半径成比例**（周长 ∝ 半径，否则内圈挤成一坨）；
  环内**按社区排序**（同一模块占同一角度扇区，边交叉大幅减少）；
  配很淡的参考圆让"几圈"看得见。实测 146 节点 → 3 环 28/49/69，
  degree 区间 10-26 / 4-10 / 2-4 严格不重叠。
  **顺带干掉了力导向布局那个随机源。**
- **图上的 hub 是子网络的 hub。** 环序按**过滤后**子网络的 degree 排，核心环是
  MAD2L1、CENPA、CCNB2、MCM6、MELK 这一整个增殖模块；而 `hub_genes.csv` 排的是
  **全网络**，前列是 GAPDH、CD34、IGF1、KIT。两个排名都对，回答的是不同问题 ——
  **副标题必须写明环序来自哪个网络**，否则读者会把核心环当成 hub 基因的答案。
- **布局要落盘。** 光看 PNG 分不出"第 3 环是不是真的在外圈"，从像素反推是猜。
  `ppi_plot_layout.csv` 写下每个节点的环号、半径、角度、坐标、degree 与模块。
- **`enrichGO` 的 `geneID` 列基因顺序不稳定。** 实测 234 行里 230 行的排列不同，
  而 p 值、计数完全一致 —— 科学上可复现，但字节不可复现，导致"两轮是否一致"
  这种验证做不了。`normalise_gene_lists()` 落盘前排序消除该差异。
- **颜色判据是量化的：色相相差 15° 以内视为同一个颜色。** 按这条算，
  本仓库先后暴露了四个凭眼睛看不出来的问题：p 值直方图的 `firebrick` 与 up 红
  只差 **0.4°**（就是同一个颜色）；旧分类色板的棕在 protanopia 下与 up 红距离 **0.002**
  （几乎重合）；magma 序列色中段与 up 红只差 **2.9°**；**viridis 暗端 `#365C8D`
  与 down 蓝只差 3.1°**（深色点会被读成"下调"）。
  色板来源固定：方向 = **ColorBrewer RdBu** 两端，连续 = **viridis 截去暗端**，
  分类 = Okabe-Ito 变体穷举出的 4 色。改色值前跑 `node tools/check_palette.mjs`。
- **颜色不能是唯一的语义载体。** PCA 的分组同时用颜色和形状，密度图用颜色和线型，
  火山图的显著性用 alpha + 大小（见下一条）。
- **火山图"没有差异基因"通常不是倍数阈值的问题。** 判据是
  `adj.P < 0.05` **且** `|log2FC| > 阈值`，是 AND —— 调小倍数阈值只会让集合更小。
  n=6 时全基因组 BH 校正下没有任何基因能过第一关（GSE64790 最小 adj.P = 0.394）。
  这时**不要**把名义显著的基因统一染成一个独立颜色（那样图上就没有上调下调了，
  而方向正是火山图存在的理由）。正确做法是**两个正交通道**：
  颜色编码方向（up 红 / down 蓝），alpha + 大小编码置信度（FDR 显著 = 实心大点，
  名义显著 = 半透明小点），图例写明 `NOT FDR-significant`。
  纵轴统一用 raw P：混用 raw P 和 adj.P 会让名义显著的点画在 FDR 线之上，
  看着像通过了 FDR。
- **JSON 键名不要带点号。** `frac_p_below_0.05` 这种键在 JS 里无法用属性访问
  （`summary.frac_p_below_0.05` 会解析成 `summary.frac_p_below_0` 再加 `.05`，直接语法错误），
  必须写成 `summary['frac_p_below_0.05']`。用 `frac_p_lt_0p05` 这类键名省事。
  同理，`deg_table.csv` 的 `P.Value`、`p.adjust` 列在 JS 里也都要用方括号访问。
- **STRINGdb 首次运行要下载约 100 MB 网络文件。** 失败时本流水线回退到基于表达谱的
  **共表达网络**（Spearman |r| ≥ 0.9），`ppi_status.json` 会标记 `coexpression_fallback`。
  共表达不是物理互作，**不得当作 PPI 证据引用**。
- **n<10 时几乎不可能有基因通过 FDR，这不是 bug。** GSE64790（3 vs 3 配对）实测：
  1,456 个基因 `raw P < 0.05` 且 `|log2FC| > 1`，但最小 `adj.P` 是 **0.394** ——
  在 16,487 个基因上做 BH，最小 raw P 要到 ~3e-6 才够，实测最好只有 6e-5。
  富集与 PPI 因此走 `ranked_fallback`（raw P 前 500 个），`deg_mode` 会写进两个
  status JSON。**此时只能说"在最显著的 N 个基因里富集到……"**，
  不能说"显著差异基因富集到……"。这是 spec 里"样本 < 10"的固有限制。
- **背景集选择会改变富集结果。** `enrichment.universe: genome`（默认，spec 要求）
  与 `detected`（仅实测基因）结果不同，报告中必须写明用了哪个。
- **不要为 tumor/normal 差异编造机制解释。** 肿瘤 vs 全组织正常的差异主要来自
  组织成分（上皮/脂肪/基质/免疫），不是肿瘤特异驱动事件。见设计文档 §2.2。
- **配对关系要显式声明，不要解析标题后缀。** GSE64790 是"同一患者肿瘤 + 正常组织"
  配对设计（年龄 72/72、41/41、52/52），`paired: true` 且 `pairs` 在 config 里逐对写死。
  猜后缀的做法换数据集就失效，而且**猜错了不会报错**，只会让配对分析静默变错。
  第 00 步会校验：每个 GSM 只出现一次、每对一例 numerator 一例 denominator、无落单样本。

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
tools/
├── check_r_syntax.mjs       R 静态检查（括号配平、配置键一致性）
├── check_sample_structure.mjs  样本相关结构预检（Node，不需要 R）
└── check_figures.mjs        图不是空白的（独立解码 PNG 像素，不需要 R）
references/troubleshooting.md 运行期故障排查
.github/workflows/geo_analysis.yml
```
