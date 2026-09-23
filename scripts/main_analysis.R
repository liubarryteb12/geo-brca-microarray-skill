#!/usr/bin/env Rscript
# ============================================================================
# main_analysis.R — 流水线编排器
# ============================================================================
# spec 的 workflow 步骤。按顺序执行 00 -> 05，每步独立 tryCatch：
#   * 必需步骤失败 -> 记录后中止（后续步骤依赖它）
#   * 可选步骤失败 -> 记录后继续
# 最后按 spec 的 acceptance_criteria 逐项校验产物，并写 results/state.json。
#
# 用法: Rscript scripts/main_analysis.R [--config assets/config.yml]
# 退出码: 0 = 全部必需步骤与验收项通过；1 = 有必需项失败
# ============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
})

# ---- bootstrap -------------------------------------------------------------
local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  here <- if (length(file_arg) > 0L) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
  } else {
    file.path(getwd(), "scripts")
  }
  assign(".geo_scripts_dir", here, envir = globalenv())
  if (exists("load_config", mode = "function")) return(invisible(NULL))
  cand <- c(file.path(here, "lib", "common.R"),
            file.path(here, "scripts", "lib", "common.R"),
            file.path(here, "..", "lib", "common.R"))
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L) stop("找不到 lib/common.R；请从仓库根目录运行")
  source(hit[[1L]])
})

# 让被 source 的步骤脚本不要自动执行 main()
options(geo.orchestrated = TRUE)

for (f in c("00_validate_inputs.R", "01_download_clean.R", "02_qc_pca_correlation.R",
            "03_deg.R", "04_heatmap_enrichment.R", "05_ppi.R",
            "06_wgcna.R", "07_lasso.R", "08_tf_regulation.R",
            "09_export_targets.R", "10_survival_diagnostics.R")) {
  p <- file.path(.geo_scripts_dir, f)
  if (!file.exists(p)) stop(sprintf("缺少步骤脚本: %s", p))
  source(p)
}

# ---- 步骤定义 --------------------------------------------------------------
# required = TRUE 的步骤失败会让整个运行以非零码退出
#
# **06/07 是可选步骤，但"可选"不等于"可以静默不做"** —— 两个脚本都会在
# 不适用时写一份 wgcna_status.json / lasso_status.json 说明原因
# （design_mode 不对、事件数不够、终点没配）。验收项检查的是**那份记录**，
# 不是"有没有出图"。所以 n=6 的 GSE64790 也能通过验收，而日志里说清了为什么。
STEPS <- list(
  list(id = "validate_inputs",    fn = run_00_validate_inputs,     required = TRUE),
  list(id = "download_clean",     fn = run_01_download_clean,      required = TRUE),
  list(id = "qc_pca_correlation", fn = run_02_qc_pca_correlation,  required = TRUE),
  list(id = "limma_deg",          fn = run_03_deg,                 required = TRUE),
  list(id = "deg_heatmap",        fn = run_04a_heatmap,            required = TRUE),
  list(id = "go_kegg_enrich",     fn = run_04b_enrichment,         required = FALSE),
  list(id = "ppi_string",         fn = run_05_ppi,                 required = FALSE),
  list(id = "wgcna",              fn = run_06_wgcna,               required = FALSE),
  list(id = "lasso_cox",          fn = run_07_lasso,               required = FALSE),
  list(id = "tf_regulation",      fn = run_08_tf_regulation,       required = FALSE),
  # §1.6 的时间依赖 AUC（timeROC）与校准（rms）。**只读 lasso_risk_scores.csv**，
  # 不动 07 的任何一行 —— 所以 07 已经验证过的 C-index 不可能因为这一步而变。
  list(id = "survival_diagnostics", fn = run_10_survival_diagnostics, required = FALSE),
  # §1.7/§1.8 的 Part 1 侧：把候选靶基因整理成 CSV 交给 Part 2。
  # 放在最后 —— 它要汇总前面所有步骤的产物。
  list(id = "export_targets",     fn = run_09_export_targets,      required = FALSE)
)

# ---- 模块零：运行清单（§0.3 / §0.4）---------------------------------------
#
# 规范来源：用户整合文档「模块零：语言与运行时规范」。
# 产物：results/<GSE>/run_manifest.json
#
# **为什么不塞进 state.json：** state.json 记的是"这一步跑没跑成"，每步重写；
# manifest 记的是"本轮是在什么条件下跑出来的"，是证据，写入后不该再变。
#
# 文档 §1「本部分人工复核节点」。**默认 pending，不是 confirmed** ——
# 自动化流水线不能替人签字，把未确认的节点记成已确认，等于把复核节点
# 变成摆设。验收里作为**可见但不阻断**的项列出（required 只标"这节点
# 是否适用本数据集"）。
HUMAN_REVIEW_NODES <- list(
  list(id = "geo_availability",     label = "GEO 数据集可用性终判",
       required = TRUE),
  list(id = "group_labels",         label = "分组标签推断结果确认",
       required = TRUE),
  list(id = "outlier_removal",      label = "任何 outlier 样本删除决定",
       required = TRUE),
  list(id = "signature_genes",      label = "预后模型最终基因集确定",
       required = FALSE),
  list(id = "virtual_perturbation", label = "虚拟敲除/过表达靶基因的生物学合理性",
       required = FALSE)
)

# 需要登记哈希的输入（相对 data_dir）。(文件名, 中文说明, 是否必需)
INPUT_FILES <- list(
  list(f = "expr_raw.rds",     d = "原始表达矩阵", required = TRUE),
  list(f = "expr_clean.rds",   d = "清洗后矩阵",   required = TRUE),
  list(f = "group.csv",        d = "分组表",       required = TRUE),
  list(f = "clean_stats.json", d = "清洗统计",     required = TRUE),
  list(f = "feature_mode.json", d = "特征模式",    required = TRUE),
  list(f = "clinical.csv",     d = "临床表",       required = FALSE)
)

# ---- §0.4 Agent 决策链 ------------------------------------------------------
#
# **决策链不是日志。** 日志随 CI 滚动消失，而"为什么这样选"是
# **结论的适用范围** —— 半年后拿到 deg_table.csv 的人必须能直接看到
# 门禁怎么判的、哪些点名工具真跑了、哪些时间点为什么没报。
#
# **只汇总已落盘的状态，不重新推理。** 这里如果自己再判一遍，
# 就会和真正的执行结果分叉 —— 而分叉出来的那份看起来同样合理。
# 所以每条 answer 都从状态 JSON 里**取**，取不到就如实写"取不到"。
record_geo_decisions <- function(cfg) {
  res <- cfg$output$results_dir
  rd <- function(f) {
    p <- file.path(res, f)
    if (!file.exists(p)) return(NULL)
    tryCatch(jsonlite::fromJSON(p, simplifyVector = FALSE), error = function(e) NULL)
  }
  # 00_validate_inputs.R 把门禁结论写在 **data/**（输入侧的记录），不是 results/
  rd_data <- function(f) {
    p <- file.path(cfg$output$data_dir, f)
    if (!file.exists(p)) return(NULL)
    tryCatch(jsonlite::fromJSON(p, simplifyVector = FALSE), error = function(e) NULL)
  }
  g <- function(x, k, default = NULL) if (is.null(x)) default else (x[[k]] %||% default)
  n_of <- function(x) if (is.null(x)) 0L else length(x)

  # 1. 门禁与设计模式 —— 唯一防止"用错数据得出结论"的机制（硬性规则 1）
  gm <- rd_data("geo_metadata.json")
  record_decision(
    cfg, "design_mode",
    "这个数据集该走 small_sample 还是 cohort 门禁？",
    sprintf("design_mode = %s（阈值 adj.P<%g, |log2FC|>%g）",
            cfg$design_mode, cfg$thresholds$adj_p, cfg$thresholds$log2fc),
    evidence = sprintf(
      paste0("门禁结论来自 00_validate_inputs.R（%s）。",
             "两种设计的降级路径与措辞约束不同，**不靠调阈值迁就数据集** —— ",
             "规模不够就跑 small_sample，够就显式写 cohort。"),
      if (is.null(gm)) "geo_metadata.json 取不到"
      else sprintf("n_samples=%s，每组 %s，分组 %s",
                   g(gm, "n_samples", "?"), g(gm, "min_per_group", "?"),
                   paste(unlist(g(gm, "group_counts", list())), collapse = "/")))
  )

  # 2. 差异表达模式 —— 决定下游措辞能不能说"显著差异基因"
  #
  # **没有独立的 deg_mode 文件**，它散在几个步骤的状态里，而且**两套取值**：
  #   04/05/07 走 `sel$mode` -> "fdr" / "ranked_fallback"
  #   08 自己判 -> "significant" / "ranked_fallback"
  # 只有回退那一个是同名同义的，所以判据只认 `ranked_fallback` 这个字符串。
  # 取不到就说取不到 —— **不猜**，下游措辞按最保守的来。
  tf0 <- rd("tf_status.json")
  deg_mode <- "unknown"
  deg_src <- ""
  for (cand in c("enrichment_status.json", "ppi_status.json",
                 "wgcna_status.json", "tf_status.json")) {
    v <- g(rd(cand), "deg_mode")
    if (is.character(v) && length(v) == 1L && nzchar(v)) {
      deg_mode <- v
      deg_src <- cand
      break
    }
  }
  record_decision(
    cfg, "deg_mode",
    "有没有 FDR 显著基因？下游该说『显著差异基因富集到』还是『最显著的 N 个基因里富集到』？",
    sprintf("deg_mode = %s%s", deg_mode,
            if (nzchar(deg_src)) sprintf("（取自 %s）", deg_src) else ""),
    evidence = if (identical(deg_mode, "ranked_fallback")) {
      "FDR 显著基因为 0，走排序表降级。**措辞必须收敛为『在最显著的 N 个基因里富集到……』**（硬性规则 3）—— 说成『显著差异基因富集到』是把降级藏起来。"
    } else if (identical(deg_mode, "unknown")) {
      "四个状态文件里都取不到 deg_mode（相关步骤未跑或状态缺失）—— **不猜**，下游措辞按最保守的来。"
    } else {
      "有 FDR 显著基因（取值 fdr / significant 两套，含义相同），可直接说『差异基因富集到……』。"
    }
  )

  # 3. §1.5 TRRUST / ChEA3 —— 点名工具的落地边界（硬性规则 27）
  tf <- tf0
  tr <- g(tf, "trrust")
  cmp <- g(tf, "trrust_vs_dorothea")
  tr_ok <- is.list(tr) && identical(g(tr, "status", ""), "ok")
  record_decision(
    cfg, "tf_sources",
    "§1.5 点名的 TRRUST / ChEA3 用上了没有？",
    sprintf("TRRUST %s；ChEA3 未接入（需要 httr/curl，CI 包列表里没有）",
            if (is.null(tr)) "状态取不到"
            else if (tr_ok) "已接入（作为 dorothea 的独立交叉验证）"
            else sprintf("未接入（%s）", g(tr, "reason", "理由缺失"))),
    evidence = if (isTRUE(g(cmp, "compared", FALSE))) {
      sprintf(paste0("TRRUST **是交叉验证，不替换主来源**（主来源仍是 dorothea）。",
                     "配对重叠 Jaccard=%s（共有 %s 对），重叠部分 mor 符号一致率=%s —— ",
                     "**mor 冲突是最该看的数**：两个库对同一对 TF-靶基因方向不一致时，",
                     "任何基于 mor 定符号的活性打分都要打折。"),
              g(cmp, "jaccard", "?"), g(cmp, "n_common_pairs", "?"),
              g(cmp, "mor_agreement", "?"))
    } else {
      sprintf(paste0("两个库都没有 CRAN/Bioconductor 包（CRAN 上 TRRUST / trrust / ",
                     "ChEA3 / chea3 四个名字实测全无），只能走数据/API。TRRUST 走官方 TSV + ",
                     "base R download.file（不引新依赖）；比对未成立的原因：%s"),
              g(cmp, "reason", "（状态取不到）"))
    }
  )

  # 4. §1.6 时间点取舍 —— "在哪个时间段有用"（硬性规则 29）
  sv <- rd("survival_diagnostics_status.json")
  if (!is.null(sv)) {
    record_decision(
      cfg, "survival_horizons",
      "§1.6 的时间依赖 AUC 报了哪几个时间点？没报的为什么？",
      sprintf("status=%s，报了 %s 个时间点",
              g(sv, "status", "?"), g(sv, "n_time_points", "0")),
      evidence = sprintf(
        "**时间点后事件不足 %s 个就不报 AUC**（事件少时方差极大，3 个事件能给 0.95 也能给 0.30）。时间单位不跨队列共享：主队列是天、验证队列是年，各自按事件时间分位数取。丢弃的时间点记在 roc_notes：%s",
        g(sv, "min_events_at_horizon", "?"),
        if (n_of(g(sv, "roc_notes")) == 0L) "（无）"
        else paste(utils::head(unlist(g(sv, "roc_notes")), 4), collapse = " / "))
    )
  }

  # 5. 跨部分交接 —— 本仓库唯一一处真正的跨语言转换（§0.2）
  tgt <- rd("part2_targets_status.json")
  if (!is.null(tgt)) {
    record_decision(
      cfg, "part2_handoff",
      "交给 Part 2（Python）的靶基因表里有什么、丢了什么？",
      sprintf("交出 %s 个基因，其中 %s 个带 logFC",
              g(tgt, "n_genes", "?"), g(tgt, "n_with_logfc", "?")),
      evidence = paste0(
        "只走 CSV（§0.2）。**丢失字段逐条登记在清单的 cross_language 里** —— ",
        "Part 2 用不到它们，但记下来才知道将来要用时回 Part 1 的哪张表取。",
        if (identical(as.integer(g(tgt, "n_with_logfc", 0L)), 0L))
          " **本轮没有 logFC**：Part 2 的 signature_alignment 会是空的，那是预期不是 bug。"
        else "")
    )
  }

  # 6. 可选步骤的终态 —— "没做"和"做了没问题"必须长得不一样（硬性规则 24）
  #
  # 状态文件有**两种形状**：`{"status": "..."}`（05/06/07/10）和
  # `{"go": {"status": ...}, "kegg": {...}}`（04）。无条件取 `$status`
  # 会在 04 上拿到 NULL —— 验收那边已经因为同样的形状问题崩过一次。
  one_status <- function(s) {
    if (is.null(s)) return("状态文件不存在")
    if (is.character(s$status) && length(s$status) == 1L) return(s$status)
    parts <- character(0)
    for (k in c("go", "kegg")) {
      node <- s[[k]]
      if (is.list(node) && is.character(node$status)) {
        parts <- c(parts, sprintf("%s=%s", k, node$status))
      }
    }
    if (length(parts) > 0L) paste(parts, collapse = "+") else "status 字段缺失"
  }
  opt <- list(
    list(f = "enrichment_status.json", n = "§1.3 GO/KEGG 富集"),
    list(f = "wgcna_status.json",      n = "§1.4 WGCNA 共表达模块"),
    list(f = "ppi_status.json",        n = "§1.5 PPI / 共表达网络"),
    list(f = "lasso_status.json",      n = "§1.6 LASSO/Cox 预后模型")
  )
  st_txt <- vapply(opt, function(o) sprintf("%s=%s", o$n, one_status(rd(o$f))),
                   character(1))
  record_decision(
    cfg, "optional_steps",
    "可选步骤这一轮是『真跑了』『有理由地没跑』还是『崩了』？",
    paste(st_txt, collapse = "；"),
    evidence = paste(
      "not_applicable / not_configured / too_few_events / package_missing /",
      "not_done / not_run 是**有理由地没跑**（PASS，可见）；failed / partial_error /",
      "schema_error 是**崩了**（FAIL）。两者混在一起时，崩溃会被洗成合法跳过 ——",
      "实测踩过两次（07_grn 漏 import、timeROC 缺 Surv）。")
  )

  # 7. 文档点名工具的缺口登记 —— "装不上"和"不归本仓库管"是两件事
  #
  # `KEY_PACKAGES` 已经记了"装没装"（未装的记 NA，硬性规则 26），
  # 但**没记为什么**。而"CRAN 上根本没有这个包"和"规范把这一节划给了 Part 2"
  # 的后续动作完全不同：前者等上游，后者去另一个仓库找。
  #
  # **`used` 从已落盘的状态文件里取，不在这里重新判** ——
  # 判据是 `tf_status.json` 的 `trrust.status == "ok"`。
  tr_used <- is.list(tr) && identical(g(tr, "status", ""), "ok")
  nt <- record_named_tools(cfg, used = if (tr_used) "TRRUST" else character(0))
  n_used <- sum(vapply(nt, function(x) isTRUE(x$used), logical(1)))
  n_missing_reason <- sum(vapply(
    nt, function(x) !isTRUE(x$used) && !nzchar(x$reason %||% ""), logical(1)))
  record_decision(
    cfg, "named_tools",
    "文档 §1 点名的工具，哪些真的产出了结果？没用的为什么？",
    sprintf("登记 %d 个；实际产出结果的是 %s",
            length(nt),
            if (n_used == 0L) "无" else
              paste(names(nt)[vapply(nt, function(x) isTRUE(x$used), logical(1))],
                    collapse = " / ")),
    evidence = paste0(
      "逐条理由写进清单的 `named_tools`（kind 取值三仓库同义：r_package / ",
      "not_on_cran / web_service / python_part）。",
      "**§1.7 / §1.8 的 scTenifoldKnk / PerturbNet / RegVelo 不是\"装不上\"，",
      "是规范把这两节划给了 Part 2** —— 本仓库只做 CSV 交接；",
      "§1.5 的 TRRUST 不是 R 包但**真的用上了**（官方 TSV，作为 dorothea 的交叉验证），",
      "ChEA3 是 web 服务、无本地包，未接入。",
      sprintf(" 本轮缺理由的登记：%d 个。", n_missing_reason))
  )

  m <- read_manifest(cfg)
  log_info(sprintf("决策链已登记：%d 条（写入 %s）",
                   length(m$decisions %||% list()), MANIFEST_NAME))
}

# ---- 验收项（直接对应 spec 的 acceptance_criteria）-------------------------
check_acceptance <- function(cfg) {
  res <- cfg$output$results_dir
  has <- function(f) file.exists(file.path(res, f))
  # 有些验收项针对 data/（输入侧的记录，不是分析产物）
  has_data <- function(f) file.exists(file.path(cfg$output$data_dir, f))

  read_status <- function(f) {
    p <- file.path(res, f)
    if (!file.exists(p)) return(NULL)
    tryCatch(jsonlite::fromJSON(p, simplifyVector = FALSE), error = function(e) NULL)
  }
  enrich <- read_status("enrichment_status.json")
  ppi <- read_status("ppi_status.json")
  wgcna <- read_status("wgcna_status.json")
  lasso <- read_status("lasso_status.json")
  tf <- read_status("tf_status.json")
  targets <- read_status("part2_targets_status.json")
  # 富集/PPI 允许"为空/回退"，但必须留下原因记录
  # 状态文件有两种形状：
  #   enrichment_status.json -> {"go": {"status": ...}, "kegg": {...}}
  #   ppi_status.json        -> {"status": "...", "reason": "..."}   <- status 在顶层
  # 所以 key 取到的可能是 list（再取 $status），也可能直接就是字符串。
  # 早期版本无条件写 s[[key]]$status，遇到 ppi_status.json 就报
  # "$ operator is invalid for atomic vectors" 并让整轮验收崩掉。
  documented <- function(s, key) {
    if (is.null(s) || is.null(s[[key]])) return(FALSE)
    node <- s[[key]]
    st <- if (is.list(node)) node$status else node
    !is.null(st) && !identical(st, "not_run")
  }
  # WGCNA / LASSO 的"做了"判据是**状态文件里有一个已定论的状态**。
  # 允许的终态：ok（真跑了）、not_applicable / not_configured /
  # too_few_events / package_missing / endpoint_error（有理由地没跑）。
  # 不允许：文件不存在，或者 status 还是 not_run —— 那说明脚本压根没执行到，
  # 而不是"这一步不适用"。
  # not_done 也是**已定论**：明确说了"做不了 + 为什么"，比编个代理指标好。
  settled <- function(s) {
    if (is.null(s) || is.null(s$status)) return(FALSE)
    s$status %in% c("ok", "not_applicable", "not_configured", "too_few_events",
                    "too_few_samples", "package_missing", "endpoint_error",
                    "empty_signature", "not_done")
  }

  # name 用 file.path(res, ...) 而不是字面量 "results/..." —— 产物现在按数据集
  # 分目录（results/GSE42568/），写死的前缀会让验收日志指向一个不存在的路径，
  # 排查时先被误导一次。
  chk <- function(f, required = TRUE, ok = has(f)) {
    list(name = file.path(res, f), ok = ok, required = required)
  }

  # 从 label_decisions.csv 读某张图的**某一类**标签有没有藏。
  # 返回 TRUE / FALSE，读不到就返回 NA（**不知道，不等于没藏**）。
  #
  # **必须按 `label` 过滤，不能只看 `figure`。** 一张图可以有多类标签
  # （`top50_heatmap` 现在同时记"热图基因"行名与"热图样本名（列）"列名）。
  # 原来取该图的**所有**记录再 `any(shown)`，于是"行名藏了、列名显示了"
  # 会被 `any` 判成"没藏" —— 行名对照表那条验收项就不再要求 CSV，
  # 而它正是这条检查存在的理由。**判据要匹配它真正在问的那件事。**
  label_hidden <- function(fig, what = NULL) {
    p <- file.path(res, "label_decisions.csv")
    if (!file.exists(p)) return(NA)
    d <- tryCatch(utils::read.csv(p, stringsAsFactors = FALSE),
                  error = function(e) NULL)
    if (is.null(d) || !all(c("figure", "shown") %in% colnames(d))) return(NA)
    r <- d[d$figure == fig, , drop = FALSE]
    if (!is.null(what) && "label" %in% colnames(d)) {
      r <- r[r$label == what, , drop = FALSE]
    }
    if (nrow(r) == 0L) return(NA)
    !any(r$shown)
  }

  checks <- list(
    chk("01-02-01-unit1-boxplot-before.pdf"),
    chk("01-02-01-unit2-boxplot-after.pdf"),
    chk("01-02-02-unit1-density-before.pdf"),
    chk("01-02-02-unit2-density-after.pdf"),
    chk("01-02-03-unit1-pca-plot.pdf"),
    # 椭圆坐标落盘：从 PNG 反推"椭圆画没画、多大"是猜（实测 stat_ellipse 在 n=3 时
    # 静默产出空数据，图上只有点）。有这张表就能直接核对。
    chk("pca_ellipse.csv", required = FALSE),
    chk("01-02-04-unit1-correlation-heatmap.pdf"),
    chk("correlation_matrix.csv"),
    chk("deg_table.csv"),
    chk("01-03-01-unit1-volcano-plot.pdf"),
    chk("01-04-01-unit1-top50-heatmap.pdf"),
    # DE 之后的 QC 关卡（bulk-rnaseq skill）：p 值分布要"均匀 + 0 附近有峰"。
    # 它区分"功效不足"和"设计有问题"，对 n<10 的设计尤其关键，所以是必需项。
    chk("01-03-02-unit1-pvalue-histogram.pdf"),
    # **这两个名字原来写错了。** 写的是 `GO_dotplot.pdf` / `KEGG_dotplot.pdf`，
    # 而实际产物从来没有这两个名字 —— 所以 `has(...)` 那一半**永远是 false**，
    # 这条检查一直只靠 `|| documented(...)` 通过。现在改成真实图名。
    # （改对只会把 FAIL 变成 PASS，不会新增失败：原来能过的仍然能过。）
    list(name = "GO 富集（结果或空原因）",
         ok = (has("01-04-04-unit1-up-ora-dotplot.pdf") && has("01-04-04-unit2-down-ora-dotplot.pdf")) || documented(enrich, "go"), required = TRUE),
    list(name = "KEGG 富集（结果或空原因）",
         ok = (has("01-04-05-unit1-up-ora-dotplot.pdf") && has("01-04-05-unit2-down-ora-dotplot.pdf")) || documented(enrich, "kegg"), required = TRUE),
    # GSEA 是弱功效数据集的主力方法，但基因集数据库/网络问题可能让它拿不到结果，
    # 所以按可选步骤处理，只要有记录在案的状态即可。
    list(name = "preranked GSEA / GO（结果或原因）",
         ok = has("01-04-02-unit1-gsea-go-dotplot.pdf") || documented(enrich, "gsea_go"),   required = FALSE),
    list(name = "preranked GSEA / KEGG（结果或原因）",
         ok = has("01-04-03-unit1-gsea-kegg-dotplot.pdf") || documented(enrich, "gsea_kegg"), required = FALSE),
    list(name = "PPI 网络（结果或回退原因）",
         ok = has("01-05-01-unit1-ppi-network.png") || documented(ppi, "status"),  required = FALSE),
    # 热图行名放不下时会被隐藏 —— 那就必须有一张表能还原"第 N 行是哪个基因"。
    #
    # **原来这条是空转的**：`ok = has(csv) || has(pdf)`，而 pdf 一定会产出，
    # 所以它永远 PASS，等于没写。判据改成读 label_decisions.csv 的真实决策：
    # 只有**确实藏了行名**才要求对照表；读不到决策文件则退回"有 csv 就算过"
    # （不因为决策文件缺失而误报失败，那条另有验收项管）。
    list(name = "top50_heatmap_genes.csv（行名隐藏时的对照表）",
         # 没藏行名 -> 不需要这张表；藏了 -> 必须有；读不到决策 -> 不误报失败
         # **按 label 过滤**：同一张图还记着列名的决策，不过滤的话
         # "行名藏了、列名显示了"会被 any() 判成"没藏"（见 label_hidden 注释）。
         ok = !isTRUE(label_hidden("01-04-01-unit1-top50-heatmap", "热图基因")) ||
              has("top50_heatmap_genes.csv"),
         required = TRUE),
    # 列名同理：121 个样本在 183 mm 里放不下（可容纳约 42 个），列名会被藏。
    # **藏了就必须能查回"图上第 N 列是哪个样本"** —— 规则 21 在列方向上是同一件事。
    list(name = "top50_heatmap_samples.csv（列名隐藏时的对照表）",
         ok = !isTRUE(label_hidden("01-04-01-unit1-top50-heatmap", "热图样本名（列）")) ||
              has("top50_heatmap_samples.csv"),
         required = TRUE),
    # 决策本身要落盘。日志里有同样的算式，但 CI 日志会滚掉，文件不会。
    list(name = "label_decisions.csv（标签决策，可核对为什么藏了行名）",
         ok = has("label_decisions.csv"), required = FALSE),
    list(name = "PPI 图注文件（图上写了什么，可检索）",
         ok = has("PPI_network_caption.txt") || !has("01-05-01-unit1-ppi-network.png"),
         required = FALSE),
    list(name = "WGCNA（结果或不适用原因）", ok = settled(wgcna), required = FALSE),
    # 批次与分组的混杂评估。**判据是"有记录"，不是"有批次"** ——
    # 找不到批次字段也是一个结论（说明查过），必须落盘；静默跳过才是问题。
    list(name = "batch_assessment.json（批次与分组混杂评估）",
         ok = has_data("batch_assessment.json"), required = FALSE),
    # PH 假设检验。**判据要分三种情况**，否则会在绿 job 上印出误导性的 FAIL：
    #   - 有 cox_zph.csv            -> 真跑了
    #   - ph_assumption 有终态      -> 跑了但失败/不适用，有理由
    #   - LASSO 整步就没适用        -> 无从检验，理由在 step 级别
    # 只有"LASSO 说 ok 却没有 ph_assumption"才算 FAIL —— 那说明这一段没执行到。
    list(name = "PH 假设检验 cox.zph（结果或不适用原因）",
         ok = has("cox_zph.csv") || documented(lasso, "ph_assumption") ||
              (settled(lasso) && !identical(lasso$status, "ok")),
         required = FALSE),
    list(name = "LASSO-Cox（结果或不适用原因）", ok = settled(lasso), required = FALSE),
    # TF 调控。判据同样是"有一个已定论的状态"，不是"有没有出图" ——
    # 没有 dorothea 又没有回退表时写 not_done + reason 也算定论，
    # 但**文件不存在或 status 缺失就是 FAIL**（AGENTS.md 规则 24）。
    list(name = "TF 调控（结果或不可得原因）", ok = settled(tf), required = FALSE),
    # **验收不等于验图。** 上一版只检查 tf_status.json 在不在，而那次
    # `figure_written` 其实是 FALSE（ggsave 收到 dpi=NULL 报错被 tryCatch 接住），
    # 图一张没出，验收却全绿。所以这里读 status 里的真实字段，
    # 并且**核对文件真的在磁盘上**。
    list(name = "TF 图确实产出（status 与磁盘一致）",
         ok = local({
           if (is.null(tf)) return(FALSE)          # 没有状态文件 = FAIL
           if (!identical(tf$status, "ok")) return(TRUE)  # 有理由地没做 = 通过
           figs <- unlist(tf$figures)
           # status=ok 就必须真有图，而且文件要在磁盘上
           isTRUE(tf$figure_written) && length(figs) > 0L &&
             all(vapply(figs, function(f) file.exists(file.path(res, f)),
                        logical(1)))
         }),
         required = FALSE),
    # 有了 TF 结果就必须说明它**不是**什么。这条是 honesty 类：
    # 一个 TF 富集表很容易被读成"这个 TF 在调控这些基因"，而
    # bulk 的 TF mRNA 水平与其蛋白活性经常不相关。
    list(name = "TF 状态写明了方法学限定（not_decoupler / limitations）",
         ok = !is.null(tf) && !is.null(tf$not_decoupler) &&
              !is.null(tf$limitations) && length(tf$limitations) >= 3L,
         required = FALSE),
    # ---- §1.5 TRRUST 交叉验证 -----------------------------------------------
    # **不阻断，但必须可见。** TRRUST 是独立第二个库；它没跑成时
    # "dorothea 跑通了"会被当成"§1.5 做完了"，而实际上少了独立背书。
    list(name = "TRRUST 交叉验证（§1.5）",
         ok = !is.null(tf) && !is.null(tf$trrust) &&
              identical(tf$trrust$status, "ok"),
         required = FALSE),
    # **一致性必须量化，不能只说"两个库高度一致"。**
    # 尤其是 mor 符号一致率 —— 两个库对同一对关系方向不一致时，
    # 任何基于 mor 定符号的活性打分都要打折，而这件事只有算出来才知道。
    list(name = "TRRUST 与 dorothea 的差异已量化（含 mor 符号一致率）",
         ok = local({
           if (is.null(tf) || is.null(tf$trrust_vs_dorothea)) return(FALSE)
           c <- tf$trrust_vs_dorothea
           if (!isTRUE(c$compared)) {
             # 没比成也算"已定论"，但必须给出原因（同规则 24）
             return(!is.null(c$reason) && nzchar(as.character(c$reason)))
           }
           !is.null(c$n_common_pairs) && !is.null(c$jaccard) &&
             !is.null(c$n_mor_comparable) && !is.null(c$mor_agreement)
         }),
         required = FALSE),
    # ---- §1.7/§1.8 交接给 Part 2 的候选靶基因 --------------------------------
    # **Part 1 对 §1.7/§1.8 的全部职责就是这张表。** 规范把虚拟扰动的
    # 计算放在 Part 2（Python + 单细胞），候选靶基因由 Part 1 产出；
    # §0.2 规定跨部分只走 CSV。所以这张表在不在，决定了 Part 2 是
    # "用 Part 1 的证据做扰动"还是"用它自己的调控子凑一个候选集"。
    list(name = "Part 2 候选靶基因交接表（§1.7/§1.8）",
         ok = settled(targets) &&
              (is.null(targets) || !identical(targets$status, "ok") ||
               file.exists(file.path(res, targets$output_file %||% "part2_targets.csv"))),
         required = FALSE),
    # **没有 logFC 的交接表是半成品。** Part 2 的 signature_alignment
    # （预测扰动方向 vs 疾病签名方向）全靠这一列；缺了它那一列只能是空的，
    # 而"空"和"算出来是 0"完全不同。所以单独报这一条。
    list(name = "交接表带上了 logFC（Part 2 算 signature_alignment 要用）",
         ok = local({
           if (is.null(targets) || !identical(targets$status, "ok")) return(TRUE)
           n <- targets$n_genes %||% 0
           k <- targets$n_with_logfc %||% 0
           n > 0 && k > 0
         }),
         required = FALSE),
    # 点名工具的缺口登记：**每个没用上的都必须有理由**
    #
    # 判据是"每个不可用的都有原因"，**不是"全都不可用"** ——
    # 将来某个工具能装了，这条应该自动变 PASS 而不是 FAIL
    # （姊妹项目 Python 侧同一条，见 scrna AGENTS 规则 17）。
    list(name = "§1 点名工具的缺口已登记（清单 named_tools，逐条理由）",
         ok = local({
           m <- read_manifest(cfg)
           nt <- m$named_tools
           if (!is.list(nt) || length(nt) == 0L) return(FALSE)
           # 用了的（TRRUST）和没用但有理由的都算"说清了"
           all(vapply(nt, function(x) {
             isTRUE(x$used) || nzchar(x$reason %||% "")
           }, logical(1)))
         }),
         required = FALSE),
    # 跨部分交接必须在清单里留痕（§0.2：转换前后、丢了什么字段）
    list(name = "跨部分交接已登记（清单 cross_language，含丢失字段）",
         ok = local({
           m <- read_manifest(cfg)
           cl <- m$cross_language
           if (is.null(cl) || length(cl) == 0L) return(FALSE)
           e <- cl[[length(cl)]]
           !is.null(e$src) && !is.null(e$dst) && !is.null(e$format) &&
             length(e$lost_fields %||% list()) > 0L
         }),
         required = FALSE),
    # ---- §1.6 时间依赖 AUC 与校准（10_survival_diagnostics.R）----------------
    # **C-index 是一个数，它不告诉你模型在哪个时间段有用。** 一个签名完全
    # 可能整体 C-index 0.70、而 5 年 AUC 只有 0.55。这一步就是去报那个数。
    #
    # 判据是"**要么报了，要么留下了不报的原因**" —— 事件不够时不报是
    # 正确的（时间依赖 AUC 在事件少时方差极大），但那必须写下来。
    list(name = "§1.6 时间依赖 AUC 已报或已说明不报（timeROC）",
         ok = local({
           s <- read_status("survival_diagnostics_status.json")
           if (is.null(s)) return(FALSE)
           st <- s$status
           # **`failed` / `partial_error` 是崩溃，不是"这一步不适用"。**
           # 实测踩过：timeROC 因 `Surv` 找不到而全部失败，脚本却写成
           # `too_few_events` + "事件数不足"，验收照常 PASS —— 一次崩溃
           # 伪装成了合法的跳过。所以这里显式把这两个状态判成 FAIL。
           if (identical(st, "failed") || identical(st, "partial_error")) {
             return(FALSE)
           }
           if (identical(st, "ok")) {
             return(!is.null(s$n_time_points) && s$n_time_points > 0L &&
                      file.exists(file.path(res, "time_roc.csv")))
           }
           # too_few_events / not_configured / schema_error 都算已定论，
           # 但必须给出 reason —— 否则读者分不清"没事件"和"脚本没跑到"
           !is.null(st) && !identical(st, "not_run") &&
             !is.null(s$reason) && nzchar(as.character(s$reason))
         }),
         required = FALSE),
    # **区分度和校准是两件事。** 模型可以把人排序排得很准（AUC 高），
    # 同时把每个人的绝对风险高估一倍（校准差）—— 后者直接决定
    # "要不要化疗"这类阈值判断。只报 AUC 等于只报了一半。
    list(name = "§1.6 校准已做或已说明不做（rms / Cox 基线 + KM）",
         ok = local({
           s <- read_status("survival_diagnostics_status.json")
           if (is.null(s)) return(FALSE)
           if (!identical(s$status, "ok")) return(TRUE)   # 上一条已判过
           n_cal <- s$n_calibration_rows %||% 0
           n_rms <- s$n_rms_rows %||% 0
           (n_cal > 0 || !is.null(s$calibration_notes)) &&
             (n_rms > 0 || !is.null(s$rms_notes))
         }),
         required = FALSE),
    # **报 AUC 必须报它的区间。** 点估计单独看会被当成结论：
    # 实测 AUC 0.68 在 35 个事件下的 95% CI 能跨到 0.5 附近。
    list(name = "§1.6 时间依赖 AUC 带置信区间",
         ok = local({
           s <- read_status("survival_diagnostics_status.json")
           if (is.null(s) || !identical(s$status, "ok")) return(TRUE)
           p <- file.path(res, "time_roc.csv")
           if (!file.exists(p)) return(FALSE)
           d <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
           if (nrow(d) == 0L) return(FALSE)
           # **判据是"每一行都有区间"，不是"至少一行有"。**
           # 第一版写的是 `sum(is.finite(...)) > 0L`，而实测（run 35483927979）
           # 全部 5 行都是 NA —— 那次确实 FAIL 了，但只要有一行侥幸有值就会
           # PASS。所有时间点来自同一个估计量（同一个 SE 向量），
           # 要么全有要么全无，所以"至少一个"这个判据在结构上就是错的。
           # 报了一个时间点的 AUC 就必须给出它的区间 —— 点估计单独看没有意义。
           all(is.finite(d$ci_low) & is.finite(d$ci_high))
         }),
         required = FALSE),
    # 这一步的方法学限定必须写明（同 §1.6 的"内部诊断不是外部验证"）
    list(name = "§1.6 诊断的方法学限定已写明（>=5 条）",
         ok = local({
           s <- read_status("survival_diagnostics_status.json")
           if (is.null(s)) return(FALSE)
           length(s$limitations %||% list()) >= 5L
         }),
         required = FALSE),
    # ---- 差距清单补图验收（T-01 Phase C2）----
    # **"整步不适用"要判 PASS，不是 FAIL**（AGENTS 规则 24 的三态原则）。
    # 实测 GSE64790 是 small_sample 数据集、不做生存分析 →
    # survival_diagnostics_status.json 写着 not_applicable，
    # 而这几项原判据 `if (is.null(s)) return(FALSE)` 把它判成红 ——
    # 那是"有理由地没跑"，不是崩溃。判据先看状态：
    #   not_applicable / not_configured / too_few_events → PASS（可见）
    #   failed / 文件缺失                                  → FAIL
    # 写成一个共享的判定函数，四项复用（避免四处各写一遍再漂移）。
    list(name = "图质量·risk plot 三联已出（G4 组）",
         ok = local({
           s <- read_status("survival_diagnostics_status.json")
           if (is.null(s)) return(FALSE)
           if (identical(s$status, "not_applicable") ||
               identical(s$status, "not_configured") ||
               identical(s$status, "too_few_events")) return(TRUE)
           figs <- s$figures %||% character(0)
           sum(c("01-10-05-unit1-risk-scores.png", "01-10-05-unit2-survival-time.png")
               %in% figs) == 2L
         }),
         required = FALSE),
    list(name = "图质量·forest plot 已出且带 HR 表",
         ok = local({
           s <- read_status("survival_diagnostics_status.json")
           if (is.null(s)) return(FALSE)
           if (identical(s$status, "not_applicable") ||
               identical(s$status, "not_configured") ||
               identical(s$status, "too_few_events")) return(TRUE)
           figs <- s$figures %||% character(0)
           "01-10-06-unit1-forest-plot.png" %in% figs &&
             file.exists(file.path(res, "cox_univariate_hr.csv"))
         }),
         required = FALSE),
    list(name = "图质量·nomogram/DCA 已出或已说明原因",
         ok = local({
           s <- read_status("survival_diagnostics_status.json")
           if (is.null(s)) return(FALSE)
           if (identical(s$status, "not_applicable") ||
               identical(s$status, "not_configured") ||
               identical(s$status, "too_few_events")) return(TRUE)
           figs <- s$figures %||% character(0)
           "01-10-04-unit1-dca.png" %in% figs ||
             (!is.null(s$nomogram_dca_note) && nzchar(s$nomogram_dca_note))
         }),
         required = FALSE),
    list(name = "图质量·KM 图带 HR(95%CI) 注释",
         ok = local({
           s <- read_status("survival_diagnostics_status.json")
           if (is.null(s)) return(FALSE)
           if (identical(s$status, "not_applicable") ||
               identical(s$status, "not_configured") ||
               identical(s$status, "too_few_events")) return(TRUE)
           # HR 在 KM subtitle 里 —— 判据是 cox_univariate_hr.csv 的 risk 行存在（forest 步骤落盘）
           p <- file.path(res, "cox_univariate_hr.csv")
           if (!file.exists(p)) return(FALSE)
           d <- utils::read.csv(p)
           any(d$term == "risk score")
         }),
         required = FALSE)
  )

  # deg_table.csv 必须含 spec 要求的四列
  deg_path <- file.path(res, "deg_table.csv")
  if (file.exists(deg_path)) {
    cols <- colnames(utils::read.csv(deg_path, nrows = 1L, check.names = FALSE))
    need <- c("gene", "logFC", "P.Value", "adj.P.Val")
    miss <- setdiff(need, cols)
    checks[[length(checks) + 1L]] <- list(
      name = sprintf("deg_table.csv 含 %s", paste(need, collapse = "/")),
      ok = length(miss) == 0L, required = TRUE
    )
  }
  checks
}

# ---- 主流程 ----------------------------------------------------------------
main <- function() {
  t0 <- Sys.time()
  cfg <- load_config()
  ensure_dirs(cfg)

  log_info("############################################################")
  log_info(sprintf(" GEO 数据挖掘流水线启动 | %s | %s",
                   cfg$dataset_id, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  log_info(sprintf(" 对比: %s vs %s | 阈值: adj.P<%g, |log2FC|>%g",
                   cfg$contrast[1], cfg$contrast[2],
                   cfg$thresholds$adj_p, cfg$thresholds$log2fc))
  log_info("############################################################")

  # ---- 模块零：建立本轮运行清单（§0.3 / §0.4）-----------------------------
  # **必须在任何步骤之前建，且先清掉上一轮** —— 清单描述的是本轮。
  init_manifest(cfg)
  capture_versions(cfg)
  record_params(cfg, list(
    seed = cfg$analysis$seed,
    design_mode = cfg$design_mode,
    contrast = cfg$contrast,
    thresholds = cfg$thresholds,
    dataset_id = cfg$dataset_id
  ))
  for (nd in HUMAN_REVIEW_NODES) {
    record_human_review(cfg, nd$id, required = nd$required, status = "pending",
                        note = sprintf("%s —— 需人工确认，本轮自动化未确认", nd$label))
  }
  log_info(sprintf("运行清单：%s", manifest_path(cfg)))

  aborted <- FALSE
  for (step in STEPS) {
    # 一旦某个必需步骤失败，后续所有步骤都失去输入，无论必需与否都跳过，
    # 否则会看到一串"文件不存在"的次级报错，掩盖真正的原因
    if (aborted) {
      log_warn(sprintf("跳过步骤 %s（前序必需步骤已失败）", step$id))
      record_step(cfg, step$id, "skipped", required = step$required,
                  message = "upstream required step failed")
      next
    }
    log_info(sprintf("---- 开始 %s ----", step$id))
    t <- Sys.time()
    err <- tryCatch({
      step$fn(cfg)
      NULL
    }, error = function(e) conditionMessage(e))
    secs <- as.numeric(difftime(Sys.time(), t, units = "secs"))

    if (is.null(err)) {
      log_info(sprintf("---- %s 完成 (%.1fs) ----", step$id, secs))
      record_step(cfg, step$id, "ok", seconds = secs, required = step$required)
    } else {
      log_error(sprintf("---- %s 失败 (%.1fs): %s ----", step$id, secs, err))
      record_step(cfg, step$id, "failed", seconds = secs, message = err, required = step$required)
      if (isTRUE(step$required)) aborted <- TRUE
    }
  }

  # ---- 模块零：登记输入哈希（§0.4）---------------------------------------
  # 放在所有步骤之后 —— 可选步骤的产物这轮有没有，跑完才知道。
  for (inf in INPUT_FILES) {
    record_input(cfg, file.path(cfg$output$data_dir, inf$f),
                 label = sprintf("%s (%s)", inf$d, inf$f),
                 required = inf$required)
  }
  msum <- manifest_summary(cfg)
  log_info(sprintf("输入登记 %d 项%s", msum$n_inputs,
                   if (length(msum$inputs_missing) > 0L) {
                     sprintf("，缺失 %d 项：%s", length(msum$inputs_missing),
                             paste(msum$inputs_missing, collapse = ", "))
                   } else "，全部就位"))

  # ---- §0.4 Agent 决策链 -------------------------------------------------
  #
  # **决策链不是日志。** 日志随 CI 滚动消失，而"为什么这样选"是
  # **结论的适用范围** —— 半年后拿到 deg_table.csv 的人必须能直接看到
  # 门禁怎么判的、哪些点名工具真跑了、哪些时间点为什么没报，
  # 而不是去翻几十个 JSON 自己拼。
  #
  # 三条判据都来自**已落盘的状态文件**，不是重新推理一遍：
  # 这里只做汇总，不做二次判断 —— 二次判断会和真正的执行结果分叉。
  record_geo_decisions(cfg)

  # ---- 标签决策在这里兜底落盘 ---------------------------------------------
  # **`write_label_decisions()` 原来只在 `04b` 的最后一行调用**，而标签决策
  # 是在**很多步**里产生的（02 的箱线图/PCA/相关性热图、04a 的 top50 行列名、
  # 04b 的富集点图）。04b 是可选步骤 —— 它一旦失败，前面几步的决策就全丢了，
  # 而验收那边靠 `label_decisions.csv` 判断"行名/列名藏了没有、要不要对照表"，
  # 文件没了就会静默退回"读不到决策 -> 不误报失败"，检查等于失效。
  # 在这里再写一次：`.label_decisions` 是跨步骤累积的，所以这一份必然包含全部。
  # （04b 里那次保留 —— 单独跑 04 的脚本时仍然需要它。）
  write_label_decisions(file.path(cfg$output$results_dir, "label_decisions.csv"))

  # ---- 验收 --------------------------------------------------------------
  log_info("=== 验收检查 ===")
  checks <- check_acceptance(cfg)

  # ---- 模块零：运行清单验收（§0.3 / §0.4）--------------------------------
  # 清单缺项不是"分析错了"，而是"这轮跑出来的东西没法追溯"。
  # **人工复核未确认不算失败** —— 默认就是 pending，那是设计如此；
  # 把它判成 FAIL 会让每个 job 都红，反而没人看。但必须可见。
  checks[[length(checks) + 1L]] <- list(
    name = "运行清单存在（run_manifest.json）",
    ok = isTRUE(msum$present), required = TRUE)
  if (isTRUE(msum$present)) {
    checks[[length(checks) + 1L]] <- list(
      name = sprintf("版本记录非空（installed.packages 全量，%d 个）",
                     msum$n_versions),
      ok = msum$n_versions >= 20L, required = TRUE)
    checks[[length(checks) + 1L]] <- list(
      name = sprintf("输入哈希已登记且必需项无缺失（%d 项，可选缺失 %d）",
                     msum$n_inputs, length(msum$inputs_missing) -
                       length(msum$inputs_missing_required)),
      ok = msum$n_inputs >= length(INPUT_FILES) &&
        length(msum$inputs_missing_required) == 0L, required = TRUE)
    checks[[length(checks) + 1L]] <- list(
      name = sprintf("人工复核节点待确认（%d 个，不阻断 job）",
                     length(msum$human_review_pending)),
      ok = TRUE, required = FALSE)
  }

  for (chk in checks) {
    log_info(sprintf("  [%s] %s%s", if (isTRUE(chk$ok)) "PASS" else "FAIL",
                     chk$name, if (isTRUE(chk$required)) "" else " (optional)"))
  }
  failed_required <- Filter(function(c) !isTRUE(c$ok) && isTRUE(c$required), checks)
  failed_optional <- Filter(function(c) !isTRUE(c$ok) && !isTRUE(c$required), checks)

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  state_path <- file.path(cfg$output$results_dir, "state.json")
  state <- if (file.exists(state_path)) {
    jsonlite::fromJSON(state_path, simplifyVector = FALSE)
  } else {
    list(steps = list())
  }
  state$acceptance <- lapply(checks, function(c) list(name = c$name, ok = c$ok, required = c$required))
  state$acceptance_failed_required <- vapply(failed_required, function(c) c$name, character(1))
  state$acceptance_failed_optional <- vapply(failed_optional, function(c) c$name, character(1))
  state$total_seconds <- round(elapsed, 1)
  state$dataset_id <- cfg$dataset_id
  # **让产物自带出处。** 同一个文件名会在不同 commit 上产出不同内容，
  # 光看名字分不出手上这份是哪个版本 —— 实测就是这么把两轮 CI 的图当成同一张的。
  # 工作区的图按 <名字>__<短SHA>.png 命名；这里把出处写进 JSON，
  # 这样从 artifact zip 里解出来的结果也能自证版本，不需要比对哈希。
  state$build <- list(
    commit = Sys.getenv("GITHUB_SHA", unset = NA_character_),
    commit_short = substr(Sys.getenv("GITHUB_SHA", unset = ""), 1, 7),
    run_id = Sys.getenv("GITHUB_RUN_ID", unset = NA_character_),
    run_attempt = Sys.getenv("GITHUB_RUN_ATTEMPT", unset = NA_character_),
    ref = Sys.getenv("GITHUB_REF_NAME", unset = NA_character_),
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  )
  state$finished_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  state$manifest <- msum
  write_json(state_path, state)

  log_info("############################################################")
  log_info(sprintf(" 流水线结束 | 用时 %.1fs | 必需项失败 %d | 可选项失败 %d",
                   elapsed, length(failed_required), length(failed_optional)))
  log_info("############################################################")

  if (length(failed_required) > 0L) {
    log_error(sprintf("必需验收项未通过: %s",
                      paste(vapply(failed_required, function(c) c$name, character(1)), collapse = "; ")))
    quit(save = "no", status = 1L)
  }
  quit(save = "no", status = 0L)
}

main()
