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
9. **引入任何随机调用都必须紧挨着它 `set.seed(cfg$analysis$seed)`。** 已知三个源：
   `impute.knn`、`fgsea`（`gseGO(seed=)` **不可靠**，必须自己设 RNG）、
   `layout_with_fr`。`analysis.seed` 不可删。**验证方式是连跑两轮比对 SHA256**，
   不是看一眼日志说"应该没问题"。

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
