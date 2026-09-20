# GEO 常规表达谱流水线（bulk 芯片）

对 GEO 人源基因表达**芯片**数据执行端到端纯生信分析，产出清洗、QC、PCA、
样本相关性、差异基因、聚类热图、GO/KEGG 富集、差异基因互作网络的**全部图表与表格**，
并在 GitHub Actions 上运行后打包为 artifact。

> **疾病无关。** 两个已验证数据集恰好是乳腺癌（GSE64790 n=6、GSE42568 n=121），
> 但 `scripts/` 里没有任何一处依赖疾病 —— 换 `assets/config.<GSE>.yml` 就换分析。
> 名字里的 `normal` 指**常规 bulk 芯片**，用来和姊妹项目
> [`scrna-pipeline-skill`](https://github.com/liubarryteb12/scrna-pipeline-skill)（单细胞）、
> [`spatial-pipeline-skill`](https://github.com/liubarryteb12/spatial-pipeline-skill)（空间转录组）区分。

## 这是一个**框架**，不是一条焊死的流水线

本仓库提供的是**生信分析的骨架与判据**：数据门禁、方法学约定、验收项、
产物清单。**具体跑什么由输入数据和配置文件决定**，步骤本身可增删 ——
加一步、换一种方法、关掉某个可选步骤，都是预期用法，不是"改坏了"。

所以「这个仓库能做什么」的答案在 `assets/config.*.yml` 和验收项里，
**不在目录结构里**。`scripts/` 中没有任何一处硬编码某个疾病或某个平台。

后续会有一个独立的「流水线编排模块」，让使用者挑选分析模块并串起来，
再与本仓库对接。**那部分不在本仓库职责范围内** —— 本仓库只负责把每一步做对。

## 怎么拿到它：云端仓库是唯一真源

本 skill **不需要"安装"**，也不依赖任何一台机器上的目录。真源是 GitHub 仓库：

    https://github.com/liubarryteb12/geo-normal-pipeline-skill

要用的时候从云端拉下来：

```bash
./use.sh                          # 拉取/更新到 ~/.cache/dsh-skills/，打印路径
./use.sh --register               # 需要本机 agent 直接发现它时才加
./use.sh --ref v1.0               # 钉住某一版
SKILL=$(./use.sh --print-path)    # 只取路径，便于脚本里用
./use.sh --clean                  # 清掉缓存副本并撤销注册
```

**拉取后会校验 `SKILL.md` 存在。** 远端改名或换结构时会明确报错，
而不是安静地给一个空目录 —— 实测过：`git clone` 失败时后面的步骤照样会跑，
最后就是靠这道校验拦住的。

> `use.sh` 里每个可能失败的步骤都显式 `|| die`，**不依赖 `set -e`**。
> 实测（bash 5.3）在 `resolved="$(pull)"` 这种「函数在命令替换里」的结构下，
> 函数内部的失败不一定会中止外层脚本。出错的路径必须自己说出来。

**两套设计模式**，由配置里的 `design_mode` 选择（不是一个门禁换个阈值，见
[`AGENTS.md`](AGENTS.md) 规则 1）：

| | `small_sample` | `cohort` |
|---|---|---|
| 规模 | 总样本 < 10，每组 ≥ 3 | 总样本 ≥ 15，每组 ≥ 10 |
| 用途 | 探索性、假设生成 | 队列级验证 |
| 数据集 | [GSE64790](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64790) · 6 例（3 TNBC vs 3 **配对**正常） · GPL19612 | [GSE42568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568) · 121 例（104 癌 vs 17 正常） · GPL570 |

> ⚠️ GSE64790 在 n=6 时**没有任何基因能通过 FDR**（最小 `adj.P` = 0.394）。富集与 PPI
> 走的是明确标注的 `ranked_fallback` 降级路径，结果只能作假设生成。见
> [`EXPERIMENTAL_DESIGN.md`](EXPERIMENTAL_DESIGN.md) §2.2 与 §2.10。
>
> GSE42568 的用途正是把这个假设放到**有功效的队列**上验证。它还有完整的随访
> （OS 35 个事件、RFS 48 个事件），是同平台（GPL570）验证队列 GSE20685 的配套发现集。

---

## 快速开始

```bash
# 0. 环境
#    R 4.3+ / Bioconductor；Node 18+（仅用于数据集预检）

# 1. 先验证数据集合规（不要跳过）
node scripts/find_dataset.mjs check GSE42568
node tools/check_sample_structure.mjs GSE42568        # 分组是否与批次混杂
node tools/check_clinical_endpoints.mjs GSE42568      # 有没有随访终点、多少个事件

# 2. 配置：一个数据集一个文件，assets/config.<GSE>.yml
#    数据集、分组、阈值、design_mode 都在里面

# 3. 跑（没有默认配置，必须显式指定）
Rscript scripts/main_analysis.R --config assets/config.GSE42568.yml

# 4. 看结果（产物按数据集分目录，跑第二个数据集不会覆盖第一个）
ls results/GSE42568/
cat results/GSE42568/state.json
```

## GitHub Actions

```bash
gh workflow run geo_analysis.yml -f dataset=GSE42568
gh run watch
gh run download --name geo-results-GSE42568
```

手动触发时指定数据集；push 到 `main` 时**两个数据集各跑一个 job**。
上限 30 分钟，artifact 名带数据集，下载下来不会混。

**实测（run 35438543496）**：Install R packages 66s（增量）、Run analysis
GSE64790 244s / GSE42568 400s，暖缓存整轮 6m10s / 8m49s。全冷缓存装包约 720s、
整轮约 20 分钟 —— 上限写 30 而不是 20，因为被掐死的 job 存不下缓存，
下一轮又是冷缓存，会变成"每次都超时"的死循环。

## 产出

所有产物落在 `results/<GSE>/` 下。

| 文件 | 内容 |
| --- | --- |
| `01-02-01-unit1-boxplot-before-after.pdf` / `.png` | 标准化前后表达分布箱线图 |
| `01-02-02-unit1-density-plot.pdf` / `.png` | 标准化前后密度曲线 |
| `01-02-03-unit1-pca-plot.pdf` / `.png` + `pca_ellipse.csv` | PCA 散点，分组用**颜色 + 形状**双重编码；组内 95% 正态椭圆（半径 2.45 SD，坐标落在 `pca_ellipse.csv` 里可核对） |
| `01-02-04-unit1-correlation-heatmap.pdf` / `.png` | 样本间 Pearson 相关热图（序列色，相关性无负值故不用发散色） |
| `correlation_matrix.csv` | Pearson + Spearman 矩阵 + 离群标记 |
| `deg_table.csv` | 全基因 limma 结果（gene/logFC/P.Value/adj.P.Val） |
| `01-03-01-unit1-volcano-plot.pdf` / `.png` | 火山图。**颜色 = 方向**（up 红 `#B2182B` / down 蓝 `#2166AC`），**alpha + 大小 = 置信度**（FDR 显著实心大点，名义显著半透明小点）。纵轴统一 raw P |
| `01-04-01-unit1-top50-heatmap.pdf` / `.png` + `top50_heatmap_genes.csv` | top DEG 聚类热图（行 Z-score，euclidean + complete）。**行名放不放得下是算出来的**：一行标签要 `字号 + 2.5pt`，画布能给 `高(in) x 72 x 0.82` 点。放不下就整张不显示行名（GSE42568 实测 50 行 > 容量 45 → 隐藏），此时靠 `top50_heatmap_genes.csv` 对照，且该表是**图上的显示顺序**（行聚类自己算再传给 pheatmap，两边同一棵树） |
| `label_decisions.csv` | **每张图的标签决策落盘**：`figure / label / n_labels / capacity / height_in / fontsize / panel_frac / min_gap / shown`。日志里有同样的算式，但 CI 日志会滚掉，文件不会 —— 想回答"这张图为什么把行名藏了"直接查这张表 |
| `01-03-02-unit1-pvalue-histogram.pdf` / `.png` | DE 后 QC：p 值分布（区分"功效不足"与"模型设定错"） |
| `01-04-02-unit1-gsea-go-dotplot.pdf` / `.png` + `GSEA_GO_table.csv` | **preranked GSEA / GO BP（主力方法）** |
| `01-04-03-unit1-gsea-kegg-dotplot.pdf` / `.png` + `GSEA_KEGG_table.csv` | preranked GSEA / KEGG |
| `01-04-04-unit1-go-ora-dotplot.pdf` / `.png` + `GO_table.csv` | ORA GO BP。**按方向分面**（两个面板各有标题和条目数）—— 初版两个 x 轴标签都写成 "up in tumor"，整张图上 "down" 一次都没出现，而数据里下调比上调还多 |
| `01-04-05-unit1-kegg-ora-dotplot.pdf` / `.png` + `KEGG_table.csv` | ORA KEGG，同样按方向分面 |
| `01-05-01-unit1-ppi-network.pdf` / `.png` + `hub_genes.csv` + `ppi_edges.csv` + `ppi_plot_layout.csv` + `PPI_network_caption.txt` | STRING PPI 网络与 hub 基因。**图做过可读性过滤，并按同心圆环排布**（最大连通分量 → degree 前 200 → 最强 450 条边 → 3 环，内圈 = hub 核心；最大 4 个 Louvain 模块上色）。完整网络见 `ppi_edges.csv`，每个节点的环号/半径/角度见 `ppi_plot_layout.csv`，图注七段式（方法/输入/展示范围/视觉编码/环/陷阱/文件）另存为可检索的 txt |
| `wgcna_modules.csv` + `wgcna_module_sizes.csv` + `wgcna_soft_power.csv` + `wgcna_module_trait.csv` + `01-06-01-unit1-wgcna-scale-free-fit.pdf` + `01-06-02-unit1-wgcna-module-trait-heatmap.pdf` + `wgcna_status.json` | **WGCNA 共表达模块（可选，仅 `cohort`）**。只用肿瘤组（带上正常样本的话第一个模块必然是"肿瘤 vs 正常"轴，而 DEG 已经答过那件事）。软阈值取最小的 R²≥0.8 的 power，达不到就如实记录不假装通过；模块-性状做 **BH 校正**。实测 GSE42568：15 个模块（14 非 grey）、power=8（R²=0.858）、126 次检验中 30 个校正后显著 |
| `lasso_coefficients.csv` + `lasso_coefficients_epv.csv` + `lasso_risk_scores.csv` + `lasso_stability.csv` + `lasso_selection_frequency.csv` + `lasso_cv_curve.csv` + `01-07-01-unit1-lasso-km-training.pdf`（有外部队列时另出 `01-07-01-unit2-lasso-km-validation.pdf`）+ `lasso_status.json` | **LASSO-Cox 预后签名（可选，需随访终点）**。终点由 config 显式指定（自动配对在字段名不规整时一定配错，而配错不报错）。报**三个** C-index：训练集 / 交叉验证 / 外部验证。实测 GSE42568：16 基因签名训练集 0.879、CV 0.793、**外部验证 0.672**；EPV 合规的 3 基因版本外部 0.642。重复 CV 选出 [16,3,3,22,3] → `signature_stable: false` |
| `enrichment_status.json` | 富集模式（`fdr` / `ranked_fallback`）、GSEA 参数、去冗余阈值与原因 |
| `ppi_status.json` | PPI 方法、节点边数、绘图过滤与环参数、回退原因 |
| `tf_regulon_enrichment.csv` + `tf_activity_by_sample.csv` + `tf_activity_group_test.csv` + `01-08-01-unit1-tf-regulon-enrichment.pdf` + `01-08-02-unit1-tf-activity-group-difference.pdf` + `tf_status.json` | **TF 调控（可选）**。调控子来自 `dorothea`（置信度 A/B/C，271 个 TF、13223 条关系）。两件事分开做：**调控子富集**用 Fisher 精确检验（背景集 = 检测到的基因，不依赖样本量）；**调控子活性**是靶基因 z-score 的 mor 加权均值，再做组间比较。实测 GSE42568：87 个 TF 富集 BH<0.05（top `E2F4` p=1.4e-15、`TEAD1`、`E2F1`、`KLF5` —— 全是乳腺癌经典的增殖驱动 TF）；190 个 TF 活性组间差异（top `KDM5B` d=3.38、`FOXO1` d=−3.26、`PPARA`）。n=6 的 GSE64790 走 `ranked_fallback`，top `ESR1` p=5.1e-06，**活性组间 0 个显著**（每组 3 个样本，Wilcoxon 没有功效 —— 如实报 0，不制造显著性） |
| `state.json` | 每步执行状态 + 验收结果（唯一逐字节不可复现的产物：含耗时与时间戳） |

> **可选步骤失败不会让 CI 变红。** 06/07/08 是 `required = FALSE`，实测第一次跑时
> WGCNA 崩了而两个 job 全绿。所以这些脚本在**每一条退出路径**上都要写状态文件，
> 验收项 `settled()` 查的就是里面的 `status` 字段 ——
> 文件不存在或字段缺失一律记 FAIL。`not_done`（明确说了做不了 + 为什么）也算已定论。

> **但"状态文件在"不等于"结果是对的"。** TF 那一步第一次跑时
> `tf_status.json` 正常产出、验收全绿，而里面的 `figure_written` 其实是
> `false` —— 图一张没出（`ggsave` 收到 `dpi = NULL` 报错，被绘图层的
> `tryCatch` 接住），同时 121 个样本因为列名写成了 `sample` 而非 `gsm`
> **全部匹配失败**，组间比较静默变成空操作。所以验收里多了一条
> **读 status 的真实字段并核对文件在磁盘上**，而不只是检查文件存在。

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

> **每张图还会多一份带版本号的副本**，形如 `PPI_network__2180651.png`，
> 后缀是产出它的 commit 短 SHA。规范文件名（`01-05-01-unit1-ppi-network.png`）保留给验收和下游脚本；
> 副本是为了解决"同名文件在不同 commit 上内容不同、光看名字分不出是哪一版"的问题 ——
> 实测就是这么把两轮 CI 的图当成同一张的。
> `results/state.json` 的 `build` 字段也记了 commit / run id / 分支，
> 从 artifact 里解出来的结果可以自证出处。**用命名区分，不改图片内容。**

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

**一个数据集一个配置文件**：`assets/config.<GSE>.yml`。切换数据集就是换文件，
`--config` 必须显式指定（故意没有默认值 —— 多数据集下静默默认到其中某一个，
正是"跑错数据集"的来源）。

```yaml
design_mode: small_sample   # 或 cohort。决定走哪套门禁，见 AGENTS.md 规则 1
dataset_id: GSE64790        # 产物目录 results/<dataset_id>/ 由此派生

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
# 找候选：人源 + 芯片
node scripts/find_dataset.mjs search --disease "breast cancer" --max-samples 10

# 只列出能凑出两个 >=3 样本组的数据集（慢，会逐个拉样本元数据）
node scripts/find_dataset.mjs search --disease "breast cancer" --max-samples 10 --two-groups

# 看某个数据集的样本明细
node scripts/find_dataset.mjs samples GSE42568

# 有没有随访终点、多少个事件、EPV 换算的签名基因数上限
node tools/check_clinical_endpoints.mjs GSE42568

# 完整校验
node scripts/find_dataset.mjs check GSEXXXXX

# 查相关结构：分组是否与批次效应混杂（这一步能省掉一整轮 15 分钟的失败运行）
node tools/check_sample_structure.mjs GSEXXXXX
```

然后把 `dataset_id` 和分组字段写进 `assets/config.<GSE>.yml`
（一个数据集一个配置文件，见 [`AGENTS.md`](AGENTS.md) 规则 2）。

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
