# AGENTS.md — 仓库约定

本仓库是一个 GEO 乳腺癌小样本芯片数据挖掘流水线，同时是一个 agent skill。
改代码前先读 [`EXPERIMENTAL_DESIGN.md`](EXPERIMENTAL_DESIGN.md)。

## 硬性规则

1. **不要绕过 `00_validate_inputs.R`。** 物种 / 数据类型 / 样本量 / 分组四项门禁是
   这个流水线唯一防止"用错数据得出结论"的机制。不要为了让某个数据集跑通而放宽它。
2. **不要把"没有结果"写成"没有富集"。** KEGG 空结果、STRING 回退都必须通过
   `results/enrichment_status.json` / `results/ppi_status.json` 记录原因，
   结论中引用该原因，不得升级为生物学结论。**`deg_mode: ranked_fallback` 时尤其注意**：
   只能说"在最显著的 N 个基因里富集到……"，不能说"显著差异基因富集到……"（§2.10）。
3. **不要为 n<10 的结果编造机制解释。** 见设计文档 §2.2。肿瘤 vs 全组织正常的差异
   主要反映组织成分，不是肿瘤特异事件。
4. **不要凭空开关配对分析。** `paired: true` 必须在 `pairs` 里**显式声明**配对关系，
   且能拿出依据（年龄/患者编号/`Series_overall_design`），不能靠解析样本标题后缀。
   见设计文档 §2.1。
5. **换数据集前先查两件事**：`tools/check_sample_structure.mjs`（分组是否与批次混杂）
   和平台注释列（有没有 `GeneSymbol` 之类的基因注释）。两者任一不合格，数据集就不可用，
   见设计文档 §1.2。
6. **新增分析步骤要同时改三处**：脚本、`main_analysis.R` 的 `STEPS`、
   `check_acceptance()` 的验收项。漏掉后两处会让步骤静默不执行。
7. **富集分析的方法学判据来自 K-Dense `pathway-enrichment` skill，不要凭直觉改**：
   - 有**完整排序表**就用 preranked GSEA，**不要卡阈值跑 ORA**（初版犯过这个错，
     灵敏度差两个数量级）。排序指标用 limma 的 moderated `t`，不用 log2FC。
   - ORA 必须**按上/下调分开跑**，合并会丢掉方向（实测 `PI3K-Akt` 其实全是下调的）。
   - 背景集默认 `detected`（实测基因集），不用全基因组。
   - GO 条目必须**基因重叠去冗余**后再报告，不要罗列同一簇的近义条目。
8. **不要报 post-hoc observed power。** 要报就报固定 n 下的 MDE（敏感性分析）。
   `n < 10` 时必须看 `pvalue_histogram`：峰在 1 或 U 形说明设计有问题，
   那时候连排序表都不能用。判据来自 K-Dense `bulk-rnaseq` / `statistical-power`。
9. **引入任何随机调用都必须紧挨着它 `set.seed(cfg$analysis$seed)`。** 已知四个源：
   `impute.knn`、`fgsea`（`gseGO(seed=)` **不可靠**，必须自己设 RNG）、
   `layout_with_fr`、**`ggrepel::geom_text_repel`**（`seed` 默认是 `NA` 不是 `NULL`）。
   `analysis.seed` 不可删。**验证方式是连跑两轮比对 SHA256**，
   不是看一眼日志说"应该没问题"。
10. **不要靠设种子解决一切。** 多线程 BLAS 的归约顺序会让浮点末位分叉，
    设种子无用，只能把 `OMP_NUM_THREADS` / `OPENBLAS_NUM_THREADS` /
    `MKL_NUM_THREADS` 钉为 1。实测 `deg_table.csv` 曾出现
    `6.00193941779545e-05` vs `...546e-05`。**报"逐字节一致"之前先确认
    是什么机制在保证它** —— 之前几轮的一致有一半是运气。
11. **颜色只能有一个含义，且判据是量化的：色相相差 15° 以内视为同一颜色。**
    改任何色值前先跑 `node tools/check_palette.mjs`，它从 `common.R` 解析实际值重算。
    不要在某个脚本里就地写 `"#C1443C"` 之类的字面量 —— 一律走 `PAL$*`。
    色板来源固定：方向 = **ColorBrewer RdBu** 两端，连续 = **viridis**，
    分类 = Okabe-Ito 变体。**viridis 必须截去暗端** —— `#365C8D` 距 down 蓝仅 3.1°，
    深色点会被读成"下调"。门禁就是为这条设的。
12. **出图代码的错误不得逃逸到方法级的 `tryCatch`。** 实测踩过：画图代码因为
    图缺 `weight` 边属性而报错，被 STRING 分支的 `tryCatch` 当成"STRING 失败"接住，
    **一个画图 bug 静默换掉了分析方法**，而状态 JSON 里看着一切正常。
    绘图要单独兜住，方法本身如实记录。
13. **图上的 hub 是子网络的 hub。** PPI 图按**过滤后**子网络的 degree 排环序，
    `hub_genes.csv` 排的是**全网络**，两者前列基因不同（核心环是增殖模块，
    全网络前列是 GAPDH / CD34 / IGF1）。副标题必须写明环序来自哪个网络 ——
    否则读者会把核心环当成"hub 基因"的答案，而那个文件给的是另一批基因。
14. **布局要落盘。** 从 PNG 反推"第 3 环是不是真的在外圈"是猜。
    `ppi_plot_layout.csv` 记录每个节点的环号、半径、角度、坐标、degree 与模块，
    环结构因此是可核对的数据而不是视觉印象。

## 代码约定

- R 脚本结构：bootstrap 块 → 辅助函数 → `run_XX(cfg)` → `if (!GEO_ORCHESTRATED())` 自执行块。
  这个模式让脚本既能被 `Rscript` 单独跑，也能被编排器 `source()`。
- 所有路径来自 `cfg$output$results_dir` / `cfg$output$data_dir`，不要硬编码 `results/`。
- 日志用 `log_info` / `log_warn` / `log_error`，不要用裸 `cat()`。
- 步骤失败必须 `stop()`，让编排器捕获并记录，不要 `tryCatch` 后静默继续。
- 可选步骤（富集、PPI）的失败要在自己的状态 JSON 里留下 `reason`。

## 验证

```bash
# R 语法与括号配平（不需要 R 运行时）
node tools/check_r_syntax.mjs

# GEO 数据集合规性（门禁 + 真实分组取值）
node scripts/find_dataset.mjs check GSE64790

# 样本相关结构（是否分组与全局表达位移混杂，不需要 R）
node tools/check_sample_structure.mjs GSE64790

# 图不是空白的（独立解码 PNG 像素，不需要 R）
node tools/check_figures.mjs results

# 配色仍然"一个颜色一个含义"（解析 common.R 的实际色值重算，不需要 R）
node tools/check_palette.mjs

# 端到端（需要 R + Bioconductor）
Rscript scripts/main_analysis.R --config assets/config.yml
```

CI 在 GitHub Actions 上跑 `geo_analysis.yml`，`timeout-minutes: 20` 是硬上限。
实测：冷缓存 14m58s，暖缓存 3m33s（R 库由 `actions/cache` 缓存）。

> **改 `packages` 列表必须同时把缓存键 `rlib-<os>-bioc-vN` 递增。**
> `actions/cache` 的 key 一旦存在就不再写回，沿用旧 key 会让新装的包每次运行都被丢掉、
> 重新装一遍。当前是 `bioc-v2`（加入 `fgsea`）。
> 递增后第一次运行会因为 `restore-keys` 前缀命中旧缓存而只增量安装，约 6 分钟；
> 之后恢复暖缓存速度。

> **验收不等于验图。** `check_acceptance()` 只看文件在不在，看不出图是不是空白。
> 出图代码的静默失败（设备开了又关、绘图没执行）会产出"存在、大小正常、纯白"的图。
> 所以 CI 里额外有 `check_figures.mjs` 这一道。

## 禁止

- 提交 API key、token 或任何凭据
- 在 `results/` 或 `data/` 里提交运行产物（`.gitignore` 已排除）
- 把上游 k-dense `scientific-agent-skills` 库的内容复制进本仓库
