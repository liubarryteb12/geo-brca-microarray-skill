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
            "03_deg.R", "04_heatmap_enrichment.R", "05_ppi.R")) {
  p <- file.path(.geo_scripts_dir, f)
  if (!file.exists(p)) stop(sprintf("缺少步骤脚本: %s", p))
  source(p)
}

# ---- 步骤定义 --------------------------------------------------------------
# required = TRUE 的步骤失败会让整个运行以非零码退出
STEPS <- list(
  list(id = "validate_inputs",    fn = run_00_validate_inputs,     required = TRUE),
  list(id = "download_clean",     fn = run_01_download_clean,      required = TRUE),
  list(id = "qc_pca_correlation", fn = run_02_qc_pca_correlation,  required = TRUE),
  list(id = "limma_deg",          fn = run_03_deg,                 required = TRUE),
  list(id = "deg_heatmap",        fn = run_04a_heatmap,            required = TRUE),
  list(id = "go_kegg_enrich",     fn = run_04b_enrichment,         required = FALSE),
  list(id = "ppi_string",         fn = run_05_ppi,                 required = FALSE)
)

# ---- 验收项（直接对应 spec 的 acceptance_criteria）-------------------------
check_acceptance <- function(cfg) {
  res <- cfg$output$results_dir
  has <- function(f) file.exists(file.path(res, f))

  read_status <- function(f) {
    p <- file.path(res, f)
    if (!file.exists(p)) return(NULL)
    tryCatch(jsonlite::fromJSON(p, simplifyVector = FALSE), error = function(e) NULL)
  }
  enrich <- read_status("enrichment_status.json")
  ppi <- read_status("ppi_status.json")
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

  checks <- list(
    list(name = "results/boxplot_before_after.pdf", ok = has("boxplot_before_after.pdf"), required = TRUE),
    list(name = "results/density_plot.pdf",         ok = has("density_plot.pdf"),         required = TRUE),
    list(name = "results/pca_plot.pdf",             ok = has("pca_plot.pdf"),             required = TRUE),
    list(name = "results/correlation_heatmap.pdf",  ok = has("correlation_heatmap.pdf"),  required = TRUE),
    list(name = "results/correlation_matrix.csv",   ok = has("correlation_matrix.csv"),   required = TRUE),
    list(name = "results/deg_table.csv",            ok = has("deg_table.csv"),            required = TRUE),
    list(name = "results/volcano_plot.pdf",         ok = has("volcano_plot.pdf"),         required = TRUE),
    list(name = "results/top50_heatmap.pdf",        ok = has("top50_heatmap.pdf"),        required = TRUE),
    list(name = "GO 富集（结果或空原因）",
         ok = has("GO_dotplot.pdf") || documented(enrich, "go"),   required = TRUE),
    list(name = "KEGG 富集（结果或空原因）",
         ok = has("KEGG_dotplot.pdf") || documented(enrich, "kegg"), required = TRUE),
    list(name = "PPI 网络（结果或回退原因）",
         ok = has("PPI_network.png") || documented(ppi, "status"),  required = FALSE)
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

  # ---- 验收 --------------------------------------------------------------
  log_info("=== 验收检查 ===")
  checks <- check_acceptance(cfg)
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
  state$finished_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
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
