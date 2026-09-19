# AGENTS.md — 仓库约定

本仓库是一个 GEO 乳腺癌小样本芯片数据挖掘流水线，同时是一个 agent skill。
改代码前先读 [`EXPERIMENTAL_DESIGN.md`](EXPERIMENTAL_DESIGN.md)。

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
gh run watch --repo liubarryteb12/geo-brca-microarray-skill
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
> 重新装一遍。当前是 **`bioc-v3`**（v2 加了 `fgsea`；v3 加了 `WGCNA` /
> `glmnet` / `survival` / `matrixStats`）。
> 递增后第一次运行会因为 `restore-keys` 前缀命中旧缓存而只增量安装（实测 66s）；
> 之后恢复暖缓存速度。

> **验收不等于验图。** `check_acceptance()` 只看文件在不在，看不出图是不是空白。
> 出图代码的静默失败（设备开了又关、绘图没执行）会产出"存在、大小正常、纯白"的图。
> 所以 CI 里额外有 `check_figures.mjs` 这一道。

## 禁止

- 提交 API key、token 或任何凭据
- 在 `results/` 或 `data/` 里提交运行产物（`.gitignore` 已排除）
- 把上游 k-dense `scientific-agent-skills` 库的内容复制进本仓库
- **`library()` / `require()` / `attach()`** —— 唯一例外是规则 23 里
  `blockwiseModules` 那次，且必须 `on.exit` 立刻 detach。
  其余一律 `pkg::fun()` 写全名：attach 会遮蔽 `stats::filter` / `stats::lag` /
  `dplyr::filter` 之类的同名函数，而遮蔽**不报错**，只是让某个调用悄悄换了实现。
