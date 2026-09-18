# AGENTS.md — 仓库约定

本仓库是一个 GEO 乳腺癌小样本芯片数据挖掘流水线，同时是一个 agent skill。
改代码前先读 [`EXPERIMENTAL_DESIGN.md`](EXPERIMENTAL_DESIGN.md)。

## 硬性规则

1. **不要绕过 `00_validate_inputs.R`。** 物种 / 数据类型 / 样本量 / 分组四项门禁是
   这个流水线唯一防止"用错数据得出结论"的机制。不要为了让某个数据集跑通而放宽它。
2. **不要把"没有结果"写成"没有富集"。** KEGG 空结果、STRING 回退都必须通过
   `results/enrichment_status.json` / `results/ppi_status.json` 记录原因，
   结论中引用该原因，不得升级为生物学结论。
3. **不要为 n=9 的结果编造机制解释。** 见设计文档 §2.2。肿瘤 vs 全组织正常的差异
   主要反映组织成分，不是肿瘤特异事件。
4. **不要改 `paired: false`**，除非能拿出配对依据。见设计文档 §2.1。
5. **新增分析步骤要同时改三处**：脚本、`main_analysis.R` 的 `STEPS`、
   `check_acceptance()` 的验收项。漏掉后两处会让步骤静默不执行。

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

# GEO 数据集合规性
node scripts/find_dataset.mjs check GSE92252

# 端到端（需要 R + Bioconductor）
Rscript scripts/main_analysis.R --config assets/config.yml
```

CI 在 GitHub Actions 上跑 `geo_analysis.yml`，`timeout-minutes: 20` 是硬上限。

## 禁止

- 提交 API key、token 或任何凭据
- 在 `results/` 或 `data/` 里提交运行产物（`.gitignore` 已排除）
- 把上游 k-dense `scientific-agent-skills` 库的内容复制进本仓库
