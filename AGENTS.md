# AGENTS.md — 仓库约定

本仓库是一个 **GEO 常规表达谱（bulk 芯片）分析流水线**，同时是一个 agent skill。
疾病无关：`scripts/` 里没有一处依赖疾病或平台，换 `assets/config.<GSE>.yml` 就换分析
（已验证的两个数据集恰好是乳腺癌 GSE64790 / GSE42568）。
改代码前先读 [`EXPERIMENTAL_DESIGN.md`](EXPERIMENTAL_DESIGN.md)。

> **这是一个框架，不是一条焊死的流水线。** 步骤、参数、产物、验收项都随分析需求
> 变化；加一步、换一种方法、关掉某个可选步骤都是预期用法。见
> [`README.md`](README.md) 开头与 [`references/module0.md`](references/module0.md)。

## 硬性规则

1. **不要绕过 `00_validate_inputs.R`。** 物种 / 数据类型 / 样本量 / 分组四项门禁是
   这个流水线唯一防止"用错数据得出结论"的机制。
   **门禁分两套，由 config 的 `design_mode` 选择，不是一个门禁换个阈值：**

   | | `small_sample` | `cohort` |
   |---|---|---|
   | 总样本 | < 10 | >= 15（WGCNA 通行下限） |
   | 每组 | >= 3 | >= 10 |
   | FDR 显著基因 | 可能为 0，走 `ranked_fallback` 降级 | 应有 |
   | 措辞约束 | "在最显著的 N 个基因里富集到……" | 可直接说"差异基因富集到……" |

   **不要为了让某个数据集通过而调大 `small_sample` 的上限。** 两种设计的降级路径
   和措辞约束都不同，混成一个门禁会让契约变含糊。规模不够就跑 `small_sample`，
   够就显式写 `cohort`。
2. **一个数据集一个配置文件、一个产物目录。** 配置是 `assets/config.<GSE>.yml`，
   产物落在 `results/<GSE>/` 与 `data/<GSE>/`（目录由 `dataset_id` 在
   `load_config()` 里派生，不要在脚本或配置里写死）。
   **`parse_args()` 故意没有默认配置** —— 多数据集下静默默认到其中某一个，
   正是"跑错数据集"的来源。
3. **不要把"没有结果"写成"没有富集"。** KEGG 空结果、STRING 回退都必须通过
   `results/<GSE>/enrichment_status.json` / `ppi_status.json` 记录原因，
   结论中引用该原因，不得升级为生物学结论。**`deg_mode: ranked_fallback` 时尤其注意**：
   只能说"在最显著的 N 个基因里富集到……"，不能说"显著差异基因富集到……"（§2.10）。
4. **不要为 n<10 的结果编造机制解释。** 见设计文档 §2.2。肿瘤 vs 全组织正常的差异
   主要反映组织成分，不是肿瘤特异事件。
5. **不要凭空开关配对分析。** `paired: true` 必须在 `pairs` 里**显式声明**配对关系，
   且能拿出依据（年龄/患者编号/`Series_overall_design`），不能靠解析样本标题后缀。
   见设计文档 §2.1。
6. **换数据集前先查三件事**：`tools/check_sample_structure.mjs`（分组是否与批次混杂）、
   平台注释列（有没有 `GeneSymbol` 之类的基因注释）、以及
   `tools/check_clinical_endpoints.mjs`（有没有随访终点、有多少个事件）。
   前两者任一不合格数据集就不可用（见设计文档 §1.2）；第三个决定能不能做预后模型。
7. **做 LASSO / Cox 之前先算 EPV。** 通行判据是每个入选变量至少 10 个事件。
   实测 GSE42568 的 OS 只有 **35 个事件** → 签名超过 **3 个基因**就开始过拟合；
   GSE20685 的死亡有 83 个事件 → 上限 8 个。文献里常见的"8 基因预后签名"配
   35 个事件是过拟合。`tools/check_clinical_endpoints.mjs` 会直接把上限打出来。
8. **新增分析步骤要同时改三处**：脚本、`main_analysis.R` 的 `STEPS`、
   `check_acceptance()` 的验收项。漏掉后两处会让步骤静默不执行。
9. **富集分析的方法学判据来自 K-Dense `pathway-enrichment` skill，不要凭直觉改**：
   - 有**完整排序表**就用 preranked GSEA，**不要卡阈值跑 ORA**（初版犯过这个错，
     灵敏度差两个数量级）。排序指标用 limma 的 moderated `t`，不用 log2FC。
   - ORA 必须**按上/下调分开跑**，合并会丢掉方向（实测 `PI3K-Akt` 其实全是下调的）。
   - 背景集默认 `detected`（实测基因集），不用全基因组。
   - GO 条目必须**基因重叠去冗余**后再报告，不要罗列同一簇的近义条目。
10. **不要报 post-hoc observed power。** 要报就报固定 n 下的 MDE（敏感性分析）。
    `n < 10` 时必须看 `pvalue_histogram`：峰在 1 或 U 形说明设计有问题，
    那时候连排序表都不能用。判据来自 K-Dense `bulk-rnaseq` / `statistical-power`。
11. **引入任何随机调用都必须紧挨着它 `set.seed(cfg$analysis$seed)`。** 已知四个源：
    `impute.knn`、`fgsea`（`gseGO(seed=)` **不可靠**，必须自己设 RNG）、
    `layout_with_fr`、**`ggrepel::geom_text_repel`**（`seed` 默认是 `NA` 不是 `NULL`）。
    `analysis.seed` 不可删。**验证方式是连跑两轮比对 SHA256**，
    不是看一眼日志说"应该没问题"。
12. **不要靠设种子解决一切，也不要以为钉了线程数就逐字节可复现。**
    浮点末位分叉有**两个独立**的来源，必须分别处理：

    | 来源 | 机制 | 处理 |
    |---|---|---|
    | 线程调度 | 多线程归约的求和顺序随调度变化 | `OMP_NUM_THREADS` / `OPENBLAS_NUM_THREADS` / `MKL_NUM_THREADS` = 1 |
    | 内核分发 | OpenBLAS 在**运行期**按检测到的 CPU 型号选 SIMD 内核，向量宽度不同则归约顺序不同 | `OPENBLAS_CORETYPE: Haswell` |

    **只钉线程数是不够的 —— 这一点是实测出来的，不是推理。** 同一个 commit
    `c370bd0` 连跑 6 轮：5 轮给出一组末位，1 轮给出另一组。那一轮和其中两轮
    **同在 westus3**，所以差异既不是 commit 引起的（我曾据此误判过一次，
    还白改了一轮代码），也不是区域引起的 —— GitHub 托管 runner 的 CPU 型号在
    **同一 Azure 区域内也不单一**。

    实测差异幅度：`deg_table.csv` 最大绝对差 9.9e-14（logFC 量级 ~1）、
    `pca_ellipse.csv` 最大 1.0e-12（坐标量级 ~100），都是**约 1 个 ULP**；
    而同一批运行的 `data/` 输入文件（`group.csv` / `clean_stats.json` /
    `feature_mode.json`）逐字节一致。实测 `deg_table.csv` 还出现过
    `6.00193941779545e-05` vs `...546e-05`。

    **加上 `OPENBLAS_CORETYPE` 之后**，同一 commit 连跑 4 轮、跨 4 个 Azure 区域
    （centralus / westus / eastus2 / eastus），结果文件指纹**全部相同**。
    同一区域的 centralus 在加之前给另一组末位、加之后给这一组 ——
    说明这个变量确实生效，不是空转。

    **所以：报"逐字节一致"之前先确认是什么机制在保证它。**
    可以放心声称的是**结构与量级可复现**（16487 个基因、PC1=46.1%、
    1909 个名义显著、GSEA 1080/101、STRING 485/500、17 项验收全过 ——
    这些在两种末位下都相同）；**不要**把"浮点末位也一致"当作已验证的性质。
    验证方式是**同一 commit 连跑多轮比对 SHA256**，不是看一眼日志说"应该没问题"。
13. **颜色只能有一个含义，且判据是量化的：色相相差 15° 以内视为同一颜色。**
    改任何色值前先跑 `node tools/check_palette.mjs`，它从 `common.R` 解析实际值重算。
    不要在某个脚本里就地写 `"#C1443C"` 之类的字面量 —— 一律走 `PAL$*`。
    色板来源固定：方向 = **ColorBrewer RdBu** 两端，连续 = **viridis**，
    分类 = Okabe-Ito 变体。**viridis 必须截去暗端** —— `#365C8D` 距 down 蓝仅 3.1°，
    深色点会被读成"下调"。门禁就是为这条设的。
14. **出图代码的错误不得逃逸到方法级的 `tryCatch`。** 实测踩过：画图代码因为
    图缺 `weight` 边属性而报错，被 STRING 分支的 `tryCatch` 当成"STRING 失败"接住，
    **一个画图 bug 静默换掉了分析方法**，而状态 JSON 里看着一切正常。
    绘图要单独兜住，方法本身如实记录。
15. **图上的 hub 是子网络的 hub。** PPI 图按**过滤后**子网络的 degree 排环序，
    `hub_genes.csv` 排的是**全网络**，两者前列基因不同（核心环是增殖模块，
    全网络前列是 GAPDH / CD34 / IGF1）。副标题必须写明环序来自哪个网络 ——
    否则读者会把核心环当成"hub 基因"的答案，而那个文件给的是另一批基因。
16. **布局要落盘。** 从 PNG 反推"第 3 环是不是真的在外圈"是猜。
    `ppi_plot_layout.csv` 记录每个节点的环号、半径、角度、坐标、degree 与模块，
    环结构因此是可核对的数据而不是视觉印象。
17. **副标题必须走 `wrap_subtitle()`。** ggplot 的副标题**不换行** ——
    超出画布宽度的部分被**静默裁掉**，不是显示成省略号，所以"字没显示全"
    从图上完全看不出来。实测 PPI 副标题 455 字符、火山图 275 字符，尾巴都被切了。
    `tools/check_r_syntax.mjs` 会挡住超过 100 字符又没折行的副标题。
    **注意 `plot.subtitle = element_text(...)` 是主题设置，不算副标题文本。**
18. **图例默认在底部横排。** `theme_paper()` 已经设好。右侧图例直接吃掉图宽，
    底部横排在同样信息量下几乎不增加图幅。pheatmap 的色条固定在右侧、无位置参数，
    它是细长条、占宽有限，保留即可。
19. **不要用 `stat_ellipse()`，自己算椭圆坐标。** 实测每组 3 个样本时它**产出空数据**，
    ggplot 不报错、坐标范围也没被撑大，所以图上只有点没有椭圆而图注却写着有椭圆 ——
    "少画了一层"和"画对了"看起来一样。用 `ellipse_points()`，坐标落盘成
    `pca_ellipse.csv`，半径系数显式写在代码里。
    半径用 `sqrt(qchisq(level, 2))`（= 2.45 SD），**不用** ggplot 默认的
    `sqrt(2*qf(level, 2, n-2))` —— 后者在 n=3 时是 6.16 SD，椭圆比数据范围大 6 倍。
    代价是低估了小样本下协方差的不确定性，**图注必须写明它是视觉参考不是检验**。
20. **图例键的填充色要显式给。** 节点是 `shape = 21` + `colour = "white"`，
    而 `fill` 映射在另一个 scale 上；size 图例的键继承不到 `fill`，
    于是画成"白描边 + 无填充"，在白底上**完全隐形**（实测图例只剩标题）。
    用 `override.aes = list(fill = ...)` 补回来。
    `pretty()` 还会给出数据范围外的断点（degree 从 2 起却标 0），要过滤。
21. **标签放不放得下是算出来的，不是看出来的。** 一行标签要 `fontsize + min_gap`
    点的垂直空间，画布能给 `height_in * 72 * panel_frac` 点 —— 走
    `common.R` 的 `fits_labels()` / `label_budget()` / `decide_rownames()`，
    **决定和算式一起进日志**。放不下就整张不显示行名，不缩字号硬塞。
    实测 `top50_heatmap` 原来写的是硬编码的 `length(genes) <= 60`，和画布高度
    毫无关系：50 个基因、5pt 字号、5.75in 高时每行只剩约 3.6px 间隙。
    **行名一旦隐藏，基因身份就只剩表能提供** —— 必须同时落盘
    `top50_heatmap_genes.csv`，而且要是**显示顺序**（行聚类自己算再传给
    pheatmap，保证两边同一棵树）。规则 16 在热图上同样适用。
    **决策本身也要落盘**（`label_decisions.csv`，走 `write_label_decisions()`）：
    日志里有同样的算式，但 CI 日志会滚掉，文件不会。
    **验收项不能写成 `has(csv) || has(pdf)`** —— pdf 一定会产出，那样这条
    等于没写（实测就是这样空转了）。判据要读 `label_decisions.csv` 的真实
    `shown` 值：只有确实藏了行名才强制要求对照表。
22. **WGCNA 只用肿瘤组；LASSO 的终点必须由 config 显式指定。**
    - 带上正常样本的话，第一个模块必然是"肿瘤 vs 正常"轴，而 DEG 已经答过
      那件事了。WGCNA 要回答的是癌组织**内部**的异质性。
    - 终点自动配对在字段名不规整时一定配错，而配错不报错，只会算出错的
      C-index。实测 GSE20685 是 `event_death` 和 `follow_up_duration (years)`，
      名字里没有共同词。
    - **不要用 `survival::concordance()` 的公式接口报 C-index。** 实测它在
      `Surv(time, event) ~ risk` 下返回 1 - Harrell C（训练集 0.121 vs 交叉验证
      0.793，正好互补），而 0.121 看着像个正常数字，不会引起怀疑。
      用 `07_lasso.R` 里的 `harrell_c()`，定义写在注释里。
    - 报 C-index 要报三个：训练集、交叉验证、**外部验证**。
      只报训练集等于没验证。
    - 重复 CV 选出多少个基因的**分布**要报（实测 [16, 3, 3, 22, 3]，
      只换 foldid 就差 7 倍）。只报"最终签名 N 个基因"是把不稳定性藏起来。
23. **`blockwiseModules` 必须在临时挂载 WGCNA 的情况下调用 —— 这是本仓库唯一
    一处 `library()`，且必须 `on.exit` 立刻 detach。**
    `blockwiseModules` 内部用 `do.call(corFnc, ...)` 算 KME，而 `corFnc` 来自
    包内常量 `.corFnc = c("cor", "bicor", "cor")` —— 是个**字符串**，按名字查找。
    不 attach 任何包时它解析到 `stats::cor`，后者没有 `weights.x` / `weights.y` /
    `cosine` 参数。实测报错：
    `unused arguments (weights.x = NULL, weights.y = NULL, cosine = FALSE)`。
    **传 `corFnc = WGCNA::cor` 没用** —— 读 1.74 源码确认 `blockwiseModules`
    形参表里没有 `corFnc`（只有 `corType`），参数掉进 `...`，而 KME 那段用包内
    常量、不看 `...`，报错一字不变。从外面没有参数能改。
    挂载前先记 `"package:WGCNA" %in% search()`，避免把调用方原有状态拆掉。
    脚本其余所有调用仍然写全名。
24. **可选步骤失败不等于"这一步不适用"。** 06/07 是 `required = FALSE`，
    但它们失败时 job 仍然是绿的 —— 实测第一次跑 WGCNA 崩了，CI 全绿。
    所以两个脚本在**每一条退出路径**上都要写状态文件：
    真跑了写 `status = "ok"`，不适用写 `not_applicable` / `not_configured` /
    `too_few_events` 加 `reason`。验收项 `settled()` 检查的就是这份记录，
    文件不存在或 `status` 缺失 = FAIL。

## 25. 图幅按毫米，宽度夹在标准栏宽内

参考规范：K-Dense `scientific-visualization` skill。
姊妹项目 `scrna-pipeline-skill` / `spatial-pipeline-skill` 用**同一组数值**，
三部分文档的图幅因此可比。

**期刊栏宽是按毫米规定的**，英寸是排版软件内部单位。写英寸时"这图多宽"
要靠换算才知道，写毫米时一眼能对上投稿要求。

| 常量 | 值 | 用途 |
|---|---|---|
| `W_SINGLE` | 89 mm | 单栏 |
| `W_ONE_HALF` | 136 mm | 一栏半 |
| `W_DOUBLE` | 183 mm | 双栏（通栏）|

`mm()` 把毫米转成 `save_pdf()` 要的英寸。

**实测修之前有 8~10 英寸（203~254 mm）的图 —— 装不进任何期刊的一页。**
现在全部夹到 183 mm 以内。宽度随类别数增长的图（`ora_w`）尤其要注意。

**`fig_width` 必须等于画布宽度。** 它是副标题的**折行宽度**（规则 17），
两者不一致时折行算错、副标题被静默裁掉，而图上完全看不出来。
所以画布宽和 `fig_width` 现在用**同一个表达式**，不可能再各自漂移。

### 尺寸符号的静态检查

`check_r_syntax.mjs` 会检查两件事：

1. 用到的 `W_*` 符号在 `common.R` 里有定义 —— 挡住 `W_DOUBL` 这类拼写错
2. 用到尺寸符号的脚本确实加载了 `common.R`

**为什么需要：** R 里用到未定义的名字要到运行时才炸，静态检查看不见。
姊妹项目 Python 侧实测漏 import 一个 `W_SINGLE`，`py_compile` 照样报
"语法通过"，白跑一整轮 CI。R 侧同理。

## 26. 每轮运行必须留下可追溯的运行清单（模块零）

参考规范：三大部分整合文档的「模块零」（§0.2–§0.4）。姊妹项目
`scrna-pipeline-skill/AGENTS.md` 规则 16 有完整说明，这里只写 R 侧特有的。

`common.R` 的清单层产出 `results/<GSE>/run_manifest.json`，接口与 Python 侧
**同名同义**（`init_manifest` / `capture_versions` / `record_input` /
`record_params` / `record_decision` / `record_human_review` /
`record_cross_language` / `manifest_summary`），这样三部分的清单可以并排读。

**R 侧有三处必须不同，不是风格问题：**

1. **`NA` 要写成 JSON `null`，不能写成 `"NA"` 字符串。** `jsonlite` 默认把
   `NA` 序列化成字符串 `"NA"`，而 Python 侧写的是 `null` —— 两边不一致时
   "这个工具没装"在读的人看来是"装了一个叫 NA 的工具"。
   所以有独立的 `write_manifest()`，显式 `na = "null"`。
2. **版本要用 `utils::installed.packages()` 全量，不能用
   `sessionInfo()$otherPkgs`。** 本仓库的脚本一律 `pkg::fun()` 写全名、
   不 `attach`，所以 `sessionInfo()$otherPkgs` 是**空的** ——
   用它等于什么都没记（见「禁止」里对 `library()` 的禁令）。
3. **哈希优先 `digest::digest(algo="sha256")`，退回 `tools::md5sum`。**
   `digest` **不在 CI 的 R 包列表里**，所以实际走的是 md5 分支 ——
   这时必须把 `hash_algo` 一起记下来。**算法不同的哈希不可直接比较**，
   不写算法等于给了个无法验证的值。

**人工复核节点**（`geo_availability` / `group_labels` / `outlier_removal` /
`signature_genes` / `virtual_perturbation`）默认 `pending`，**不算失败**。
其中 `signature_genes` 与 `virtual_perturbation` 是 `required = FALSE` ——
前者只在有随访终点时才有意义，后者是保留框架（见 §1.7/§1.8）。

**`init_manifest` 必须清掉上一轮**（同规则 14 的道理）：上轮的清单冒充本轮，
比没有清单更糟。输入哈希在**步骤跑完之后**才登记 —— 可选步骤这轮有没有
产物，跑完才知道。

### 决策链要有调用点，光有函数等于没有

**接口定义了但一次都没调用，清单里的 `decisions` 就是空数组 —— 而空数组
和"这一轮没有任何决策"长得一模一样。** 实测三个仓库的
`record_decision()` 都定义了，但只有 Part 2 真的在调，
R 侧两个仓库一次都没调过。

现在 R 侧的调用点是 `main_analysis.R` 的 `record_geo_decisions()`，
放在**所有步骤跑完之后、验收之前**，记**七条**：`design_mode` / `deg_mode` /
`tf_sources` / `survival_horizons` / `part2_handoff` / `optional_steps` /
`named_tools`。

> **这段文字写错过一次**：加 `named_tools` 那条时只改了代码，
> 这里还写着"六条"并只列了六个。**数数与清单是两处，改一处就会分叉** ——
> 所以验收项不检查"有几条"，检查的是**该有的节点名在不在**
> （见 `references/module0.md` §5）。

**每条 answer 都是从已落盘的状态文件里「取」的，不是重新推理一遍。**
这一点是硬要求：自己再判一遍就会和真正的执行结果分叉，
而**分叉出来的那份看起来同样合理** —— 读者没有任何办法发现。

两个具体的坑（都写进代码注释了）：

- **`deg_mode` 有两套取值**：04/05/07 走 `sel$mode` 给
  `"fdr"`/`"ranked_fallback"`，08 自己判给 `"significant"`/`"ranked_fallback"`。
  只有回退那一个同名同义，所以判据只认 `ranked_fallback`，
  并且四个状态文件依次尝试；取不到写 `unknown`，**不猜**。
- **状态文件有两种形状**：`{"status": ...}`（05/06/07/10）与
  `{"go": {...}, "kegg": {...}}`（04）。无条件取 `$status` 在 04 上拿到
  `NULL` —— 验收那边已经因为同样的形状问题崩过一次
  （`$ operator is invalid for atomic vectors`）。

## 27. §1.5 的 TRRUST / ChEA3：能做的做，做不了的如实写

**两个库都没有 CRAN/Bioconductor 包**（实测查过 CRAN 上 `TRRUST` /
`trrust` / `ChEA3` / `chea3` 四个名字，全无）。所以不能靠 `install.packages`
解决，只能走数据/API。

| | 落地方式 | 状态 |
|---|---|---|
| **TRRUST** | 官方 TSV，base R `download.file`（**不引新依赖**），缓存到 `data/<GSE>/trrust_cache/`，记 sha256 | 已接入，作为 dorothea 的**独立交叉验证** |
| **ChEA3** | Ma'ayan Lab web 服务，调用需 `httr`/`curl`（CI 包列表里没有） | **未接入**，理由写进 `limitations` |

**TRRUST 是交叉验证，不替换主来源。** 主来源仍是 dorothea。

**比对必须报三件事，缺一不可：** 配对重叠（Jaccard）、
**重叠部分的 mor 符号一致率**、各自独有部分的大小。
只报"两个库各有几万条关系"没有信息量 —— 数量大不等于一致。
**mor 冲突是最该看的数**：两个库对同一对 TF-靶基因的调控方向就不一致时，
任何基于 mor 定符号的活性打分都要打折。冲突例子也要落盘（最多 20 条），
只给比例没法核对。

**解析必须校验格式，否则上游改格式会静默出错结果。** 实测 TRRUST v2
人类文件是 **9396 行 / 795 个 TF / 4 列无表头**（TF / 靶基因 / 模式 / PMID），
模式只有 `Activation` / `Repression` / `Unknown` 三种。列数不对或出现
未知取值时**直接放弃并 WARN**，不要硬算。

**`Unknown` 映射成 `NA`，不能映射成 `0`。** `0` 在加权平均里等于
"这条关系不贡献"，而 `NA` 会被显式跳过 —— 前者把"不知道"混进了分母。
实测 9396 条里有 **4325 条是 Unknown**（占 46%），这个区别很大。

**不要写物种分支。** 本仓库从头到尾硬编码人类（门禁 `stop()` 非人源、
`org.Hs.eg.db`、`species = 9606`），写一个从配置读物种的分支是**永远走不到
的死代码**，而且 `check_r_syntax.mjs` 会报"引用了配置里不存在的字段"
（实测报过）。真要用小鼠，先改门禁再加分支。

**算法本地可验证。** R 跑不了，但可以把解析与比对逻辑逐行转写成 Python，
在**真实 TRRUST 文件**上跑 —— 这验证算法，不验证 R 语法，两者都不能省。
实测这一步抓到过三次**测试期望写错**（TRRUST 同一对有多行 PMID 导致
去重后数量变少、`dropna` 丢掉 Unknown 导致可比数变小、`round(...,4)`
让 1e-9 容差失败），三次都不是代码错。

## 28. R 没有隐式字符串拼接 —— 别照 Python 的习惯写

```r
method = ("第一段"
          "第二段")        # <- Error: unexpected string constant
```

Python / C 会把相邻字面量接起来，**R 不会**，这是语法错误。
`c("a", "b")` 和 `paste0("a", "b")` 是逗号分隔的多参数，合法；
区别就在**有没有逗号**。跨行折字符串一律 `paste0()`。

**这个错实测炸过一次 CI（run 35482359507，`09_export_targets.R:258`），
而且 `check_r_syntax.mjs` 没挡住。** 原因是它的 `stripLiterals()` 把每个
字符串换成 `""`，于是 `("a" "b")` 变成 `("" "")`，括号配平完全正常 ——
**做语法检查的函数把语法错误本身擦掉了**。

所以现在多了一个 `checkImplicitConcat()`，**在原始源码上**查（不是
`stripLiterals` 的结果），判据是"两个字符串字面量之间除了空白和注释
什么都没有"。实测它能报出上面的写法，也不会误报 `c(...)` / `paste0(...)`
里的逗号分隔。

**教训不止于此：** 一个检查器"通过了"只说明它检查的那件事通过了。
`check_r_syntax.mjs` 自己的输出就写着"这不等于 R 能跑通"—— 这次又验证了一遍。
新写一类 R 语法结构时，顺手想想现有检查器**是不是根本看不见它**。

## 29. §1.6 的 C-index 只是半个答案：还要报"在哪个时间段有用"和"绝对风险准不准"

`07_lasso.R` 报的 Harrell C-index 是**一个数**。它答不了两个临床上真正的问题：

| 问题 | 工具 | 产物 |
|---|---|---|
| **在哪个时间点**有区分度？ | `timeROC` | `time_roc.csv`（各时间点 AUC + 95% CI） |
| 预测的**绝对风险**准不准？ | `rms` | `calibration.csv` + `rms_validate.csv` |

一个签名完全可能整体 C-index 0.70、而 5 年 AUC 只有 0.55 —— 前者会被读成
"模型有效"，后者才是决策点上的真实表现（5 年生存是乳腺癌的决策点）。

**区分度和校准是两件事。** 模型可以把所有人排序排得很准（AUC 高），
同时把每个人的风险都高估一倍（校准差）—— 后者直接决定"要不要化疗"这类
阈值判断。**只报 AUC 等于只报了一半。**

### 落在新脚本 `10_survival_diagnostics.R`，不是塞进 07

**07 有 700 多行，任何一行动到 `risk_train` 都会让 C-index 变，而 C-index
变了不容易察觉。** 新脚本只**读** `lasso_risk_scores.csv`（07 的产物），
一行都不改 07 —— 这样 07 已经验证过的数值**不可能**因为这一步而变。

新增分析优先走"读上一步的产物表"这个模式，不要为了复用中间变量去改老脚本。

### 四条不能省的约定

1. **时间单位不跨队列共享。** 主队列 `overall survival time_days` 是**天**、
   验证队列 `follow_up_duration (years)` 是**年**。同一个数字 1825 差 365 倍。
   所以时间点是**每个队列各自按事件时间分位数算的**，不是全局常量 ——
   共用会把验证队列的 AUC 算在一个几乎没人随访到的点上。
2. **观测生存用 KM，不用"有没有事件"的二元比例。** 到时间点 t 时删失的人
   不是"没发生事件"，直接算比例会**系统性低估风险** —— 而低估的方向恰好
   让校准曲线看起来"预测偏高"，是会被误读成模型保守的。
3. **校准斜率取 `index.corrected`，不是 `index.orig`。** 在训练集上算斜率
   必然接近 1（模型就是在这批人身上拟合的）；`rms::validate` 的自助法
   估计"这个 1 里有多少是乐观偏差"，减掉之后才是真实斜率。
   报 `index.orig` 是自我循环。
4. **时间点后事件不足 10 个就不报 AUC**（`MIN_EVENTS_AT_HORIZON`）。
   事件少时时间依赖 AUC 的方差极大，3 个事件能给出 0.95 或 0.30。
   实测 GSE42568 的 35 个事件只够报 **2 个**时间点（75 分位之后只剩
   约 8.75 个事件）—— **这是正确行为，不是漏了一个**，被丢弃的时间点
   写进 `notes`。

### 本地能验证的部分：把算法转写成 Python

R 只在云端跑，所以算法错误不在本地挡掉就要烧一整轮 CI。
`_smoke_survival_diag.py` 把 `pick_horizons` / KM-at-horizon / `cut` 分组
逐行转写，20 条断言。**实测 5 条第一版断言全写错了，R 代码是对的：**

- 写死"35 事件给 3 个时间点" —— 实际只够 2 个（门槛过滤生效）
- `km_at` 期望"10 个里 2 个事件 -> KM 0.8" —— 末次事件时风险集只剩 1 人，
  KM 必然降到 0；KM 不是 `1 - 事件数/n`
- `cut` 的组数写成 5 —— R 的 k 个断点给 **k-1** 个区间

**"测试失败"和"代码有错"是两件事。** 上一轮（`_smoke_export_targets.py`）
也是 3 条断言错、0 条代码错。所以失败时先核对期望是怎么来的。

### `timeROC` 用了**未声明的依赖**，必须临时挂载 `survival`

**实测（run 35483432635）：两个队列的 timeROC 全部失败，报
`could not find function "Surv"`。** 查 CRAN 上 timeROC 0.4.1 的元数据：

```
Depends:  R (>= 2.10)                 <- 没有 survival
Imports:  pec (>= 2.4.4), mvtnorm     <- 没有 survival
Suggests: survival, timereg           <- survival 在这里
NAMESPACE: import(pec); import(mvtnorm)
           **没有任何 importFrom(survival, ...)**
```

它内部按名字调用 `Surv`，却只把 `survival` 写成 `Suggests`。
正常用法能跑通是因为用户先 `library(survival)` 了，`Surv` 恰好在搜索路径上。

**`pkg::fun()` 不 attach 任何东西**（连 `Depends` 都不 attach，更别说
`Suggests`），所以本仓库的调用约定下它必然失败。从外面没有参数能改。

这是本仓库**第二处** `library()`（第一处是规则 23 的 `blockwiseModules`），
同样必须：先记 `"package:survival" %in% search()`，`on.exit` 立刻
`detach(..., unload = FALSE)`。其余调用仍写全名。

**推广：** 引入任何新 R 包前，先看它的 `NAMESPACE` 里
**`importFrom` 有没有覆盖它自己用到的函数**。只用 `Suggests` 声明的包
在本仓库的"不 attach"约定下会静默失败 —— 而失败点在包内部，
报错信息不会提到依赖声明。

### 崩溃不能伪装成"这一步不适用" —— 我自己又踩了一次

第一版的 `too_few_events` 分支**无条件**写这个状态，于是上面那次
timeROC 崩溃被记成"事件数不足"，验收照常 PASS。**这正是规则 24
（可选步骤失败不等于这一步不适用）说的那个坑。**

现在分三态：

| 状态 | 含义 | 验收 |
|---|---|---|
| `ok` | 真跑了 | 要求 `n_time_points > 0` 且 CSV 存在 |
| `too_few_events` / `not_configured` | **有理由地没跑**，且 `timeROC` 没报错 | PASS（可见） |
| `failed` / `partial_error` | **崩溃了**，`reason` 里带错误原文 | **FAIL** |

判据是"`roc_notes` 里有没有出现 `失败:`" —— 先把崩溃挑出来，
剩下的才允许落到"没跑"。**顺序反了就会把崩溃洗成合法跳过。**

**推广：** 任何"因为不适用所以没做"的分支，都要先排除"其实是崩了"。
写这类分支时问一句：**如果这一步的代码坏了，它会落到哪个状态？**

### `timeROC` 的 `inference` 是 **list**，不是 matrix

修好 `Surv` 之后（run 35483927979）AUC 算出来了，**但 5 个时间点的
`95%CI` 全是 `[NA, NA]`** —— 而 `[FAIL] §1.6 时间依赖 AUC 带置信区间`
把它抓住了。点估计全对、区间全丢，从数字上看很像"SE 太小"或"事件太少"。

第一版按 matrix 取（`rownames == "SE"` / `inf[1, ]` / `is.numeric(inf)`），
**三个分支全不匹配**。读 timeROC 0.4.1 源码（`R/timeROC_3.R` 末尾）确认：

```r
inference <- list(mat_iid_rep_2 = mat_iid_rep,       # <-> AUC_2
                  mat_iid_rep_1 = mat_iid_rep_star,  # <-> AUC_1
                  vect_sd_1     = vetc_sestar,       # <-> AUC_1 的 SE
                  vect_sd_2     = vetc_se,           # <-> AUC_2 的 SE
                  vect_iid_comp_time = ...)
```

无竞争风险时返回 `ipcwsurvivalROC`，此时 `AUC = AUC_1` → **SE 取
`vect_sd_1`**。这不是猜的：包自己的 `confint.ipcwsurvivalROC()` 第一行就是
`se <- object$inference$vect_sd_1[!is.na(object$AUC)]`。

**读源码比猜结构便宜。** 一个"看起来对但取不到值"的字段访问不报错 ——
它给 NA，而 NA 往下游走会变成"这个时间点算不出来"的样子。

### 点估计区间自己算，不用 `confint()`

`confint()` 为了**同时置信带**要跑 `n.sim = 2000` 次 `rnorm()` —— 那是随机
过程，会引入一个新的需要设种子的来源（规则 11）。点估计区间是闭式的：

```
CI = AUC ± qnorm(0.975) · SE
```

用闭式可让这一步保持确定性，不新增 RNG 源。

### 验收判据不能用"至少一个"绕过结构性缺失

第一版写的是 `sum(is.finite(ci_low) & is.finite(ci_high)) > 0L`。
那次确实 FAIL 了（5 行全 NA），**但只要有一行侥幸有值就会 PASS**。

所有时间点来自**同一个估计量**（同一个 SE 向量），要么全有要么全无 ——
所以"至少一个"这个判据在结构上就是错的。现在要求 `all(is.finite(...))`：
**报了一个时间点的 AUC，就必须给出它的区间。**

**推广：** 判据要匹配数据生成过程的**结构**。当某个量在结构上不可能部分
存在时，"至少一个"就是给了自己一条永远走得到的后门。

### 跑通之后的实测结果（run 35484353422，不得在后续改动中变化）

```
training (GSE42568, 天)   2 个时间点
  t=526   AUC=0.9099  95%CI=[0.8268, 0.9930]  风险集 94  事件 26
  t=1044  AUC=0.8769  95%CI=[0.7949, 0.9589]  风险集 81  事件 17
GSE20685 (年)             3 个时间点
  t=2.1   AUC=0.7577  95%CI=[0.6599, 0.8555]  风险集 306 事件 62
  t=4     AUC=0.7120  95%CI=[0.6300, 0.7941]  风险集 284 事件 40
  t=6.3   AUC=0.6960  95%CI=[0.6264, 0.7655]  风险集 236 事件 20

校准（rms::validate，index.corrected）
  training  Dxy 0.7505（原始 0.7584）  斜率 0.9440（原始 1.0000）
  GSE20685  Dxy 0.3372（原始 0.3423）  斜率 1.0204（原始 1.0000）
完成: 5 个时间点, 8 行校准, 2 个队列的 rms 校正
```

**这张表就是规则 29 开头那句话的证据。** 外部验证的整体 C-index 是 0.672，
但 AUC 从 2.1 年的 0.7577 掉到 6.3 年的 0.6960 —— 区间下界也从 0.6599
收到 0.6264。**"模型有效"这个说法在 2 年成立得比较硬，在 6 年只是勉强。**

而**校准斜率 1.0204 说明验证队列上绝对风险基本无偏**（训练集校正后 0.9440，
略高估）。区分度一般但校准良好 —— 这正是"只报 C-index 等于只报一半"的意思：
0.672 会被读成"模型一般"，而它其实**排序能力一般、但风险数值可以直接用**。

**训练集两个时间点的 AUC 明显高于验证队列，不要拿 0.91 说事** ——
那是同一批人上拟合的模型（EPV 只有 2.2 的那个 16 基因签名，
见规则 7）。**报 §1.6 的数必须以验证队列为准。**

## 代码约定

- R 脚本结构：bootstrap 块 → 辅助函数 → `run_XX(cfg)` → `if (!GEO_ORCHESTRATED())` 自执行块。
  这个模式让脚本既能被 `Rscript` 单独跑，也能被编排器 `source()`。
- 所有路径来自 `cfg$output$results_dir` / `cfg$output$data_dir`，不要硬编码 `results/`。
- 日志用 `log_info` / `log_warn` / `log_error`，不要用裸 `cat()`。
- 步骤失败必须 `stop()`，让编排器捕获并记录，不要 `tryCatch` 后静默继续。
- 可选步骤（富集、PPI）的失败要在自己的状态 JSON 里留下 `reason`。

## 验证

**R 代码只在云端跑 —— 开发机不装 R 是设计，不是障碍。**
本地能跑的都是静态检查；真正的执行、出图、验收在 GitHub Actions。
所以改完 R 代码的循环是：**本地静态检查 → 推 → 看日志 → 改 → 再推**。

```bash
# —— 本地（都不需要 R 运行时）——

# R 语法与括号配平 + 副标题折行 + 配置字段引用一致性
node tools/check_r_syntax.mjs

# GEO 数据集合规性（门禁 + 真实分组取值）
node scripts/find_dataset.mjs check GSE42568

# 样本相关结构（是否分组与全局表达位移混杂）
node tools/check_sample_structure.mjs GSE42568

# 有没有随访终点、有多少个事件、EPV 换算出的签名基因数上限
node tools/check_clinical_endpoints.mjs GSE42568

# 图不是空白的（独立解码 PNG 像素）；参数是**数据集目录**
node tools/check_figures.mjs results/GSE42568

# 配色仍然"一个颜色一个含义"（解析 common.R 的实际色值重算）
node tools/check_palette.mjs

# —— 云端（真正的端到端）——
gh workflow run geo_analysis.yml -f dataset=GSE42568
gh run watch --repo liubarryteb12/geo-normal-pipeline-skill
```

**本地静态检查不等于能跑通。** `check_r_syntax.mjs` 自己的输出就写着这句话。
它挡得住语法错和漏折行的副标题，挡不住"参数传错类型"和"返回的是 list 不是
向量"—— 这两类都实测发生过，只有在云端才暴露。

**每轮只解决日志明确指出的那件事。** 不要凭"看起来可能有问题"改代码：
曾把浮点末位漂移误判成 commit 引起的，白改一轮（见规则 12）。
`references/troubleshooting.md` 有完整的云端循环与排查表。

CI 在 GitHub Actions 上跑 `geo_analysis.yml`，`timeout-minutes: 30` 是硬上限。
**实测（run 35438543496）**：Install R packages 66s（增量）、Run analysis
GSE64790 244s / GSE42568 400s，暖缓存整轮 6m10s / 8m49s。
全冷缓存下装包约 720s，整轮约 20 分钟 —— 所以上限写 30 而不是 20：
被掐死的 job 存不下缓存，下一轮又是冷缓存，会变成"每次都超时"的死循环。

> **push 时两个数据集各跑一个 job**（矩阵），手动触发时只跑指定的那个。
> 产物 artifact 名带数据集（`geo-results-GSE42568`），下载下来不会混。
> 平台注释缓存落在 `data/<GSE>/geo_cache`，**CI 里不持久化**，每轮重下。

> **改 `packages` 列表必须同时把缓存键 `rlib-<os>-bioc-vN` 递增。**
> `actions/cache` 的 key 一旦存在就不再写回，沿用旧 key 会让新装的包每次运行都被丢掉、
> 重新装一遍。当前是 **`bioc-v5`**（v2 加了 `fgsea`；v3 加了 `WGCNA` /
> `glmnet` / `survival` / `matrixStats`；v4 加了 `dorothea`；
> **v5 加了 `timeROC` / `rms`**，见规则 29）。
> 递增后第一次运行会因为 `restore-keys` 前缀命中旧缓存而只增量安装（实测 66s）；
> 之后恢复暖缓存速度。

> **验收不等于验图。** `check_acceptance()` 只看文件在不在，看不出图是不是空白。
> 出图代码的静默失败（设备开了又关、绘图没执行）会产出"存在、大小正常、纯白"的图。
> 所以 CI 里额外有 `check_figures.mjs` 这一道。

## 禁止

- 提交 API key、token 或任何凭据
- 在 `results/` 或 `data/` 里提交运行产物（`.gitignore` 已排除）
- 把上游 k-dense `scientific-agent-skills` 库的内容复制进本仓库
- **`library()` / `require()` / `attach()`** —— 只有**两处**例外，都必须先记
  `"package:X" %in% search()`、`on.exit` 立刻 detach：
  1. 规则 23 的 `blockwiseModules`（`WGCNA::cor` 被包内常量按名字查找）
  2. 规则 29 的 `timeROC::timeROC`（`Surv` 被包内代码按名字调用，
     而 `survival` 只在它的 `Suggests` 里）

  其余一律 `pkg::fun()` 写全名：attach 会遮蔽 `stats::filter` / `stats::lag` /
  `dplyr::filter` 之类的同名函数，而遮蔽**不报错**，只是让某个调用悄悄换了实现。

  **每加一处例外都要问：这是"包内部按名字找东西"吗？** 是，才允许。
  只是"写全名太麻烦"不是理由。
