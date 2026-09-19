# ============================================================================
# lib/common.R — 共享工具：配置加载、日志、JSON 状态、GEO SOFT 抓取
# ============================================================================
# 被 scripts/00..05 与 main_analysis.R 共同 source()。
# 不依赖任何 Bioconductor 包，只用 base R + jsonlite + yaml。
# ============================================================================

# 编排器标志：main_analysis.R 会设为 TRUE，脚本据此决定是否自动执行 main()
GEO_ORCHESTRATED <- function() isTRUE(getOption("geo.orchestrated", FALSE))

# ---- 日志 ------------------------------------------------------------------

.geo_log <- function(level, ...) {
  msg <- paste0(..., collapse = "")
  cat(sprintf("[%s] %-5s %s\n", format(Sys.time(), "%H:%M:%S"), level, msg))
  flush(stdout())
}

log_info <- function(...) .geo_log("INFO", ...)
log_warn <- function(...) .geo_log("WARN", ...)
log_error <- function(...) .geo_log("ERROR", ...)

# ---- 配置 ------------------------------------------------------------------

#' 解析 `--config <path>` 命令行参数
#' @return 配置文件路径（默认 assets/config.yml）
parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  cfg <- "assets/config.yml"
  i <- 1L
  while (i <= length(argv)) {
    if (argv[i] %in% c("--config", "-c")) {
      if (i == length(argv)) stop("--config 需要一个路径参数")
      cfg <- argv[i + 1L]
      i <- i + 2L
    } else if (argv[i] %in% c("--help", "-h")) {
      cat("用法: Rscript <script>.R [--config assets/config.yml]\n")
      quit(save = "no", status = 0L)
    } else {
      # 允许直接传位置参数作为配置路径
      cfg <- argv[i]
      i <- i + 1L
    }
  }
  cfg
}

#' 读取并校验配置文件
#' @param path 配置 YAML 路径
#' @return 配置 list（已做必填项与取值合法性校验）
load_config <- function(path = parse_args()) {
  if (!file.exists(path)) stop(sprintf("配置文件不存在: %s", path))
  cfg <- yaml::read_yaml(path)

  required <- c("dataset_id", "group_field", "group_values", "contrast")
  missing <- setdiff(required, names(cfg))
  if (length(missing) > 0L) {
    stop(sprintf("配置缺少必填字段: %s", paste(missing, collapse = ", ")))
  }
  if (length(cfg$contrast) != 2L) {
    stop("contrast 必须是长度为 2 的向量: [分子, 分母]")
  }
  unknown <- setdiff(unlist(cfg$contrast), names(cfg$group_values))
  if (length(unknown) > 0L) {
    stop(sprintf("contrast 引用了 group_values 中不存在的组: %s", paste(unknown, collapse = ", ")))
  }
  if (!is.null(cfg$group_values) && length(cfg$group_values) < 2L) {
    stop("group_values 至少需要两个组才能做对比")
  }

  # 默认值，避免下游到处写 is.null 判断
  cfg$platform_id        <- cfg$platform_id %||% ""
  cfg$paired             <- isTRUE(cfg$paired)
  cfg$thresholds         <- utils::modifyList(
    list(adj_p = 0.05, log2fc = 1.0, string_score = 400, outlier_cor = 0.8),
    cfg$thresholds %||% list()
  )
  cfg$analysis           <- utils::modifyList(
    list(top_heatmap_genes = 50, impute_k = 10, pca_top_genes = 2000, hub_gene_count = 10),
    cfg$analysis %||% list()
  )
  cfg$enrichment         <- utils::modifyList(
    list(universe = "genome", p_adjust = "BH", pvalue_cutoff = 0.05,
         qvalue_cutoff = 0.2, top_terms = 15, ont = "BP", kegg_organism = "hsa"),
    cfg$enrichment %||% list()
  )
  cfg$output             <- utils::modifyList(
    list(results_dir = "results", data_dir = "data"),
    cfg$output %||% list()
  )
  cfg
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' 建立输出目录
ensure_dirs <- function(cfg) {
  for (d in c(cfg$output$results_dir, cfg$output$data_dir)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(TRUE)
}

# ---- JSON / 状态 -----------------------------------------------------------

#' 写 JSON（原子写：先写临时文件再改名，避免半截文件）
write_json <- function(path, obj) {
  tmp <- paste0(path, ".tmp")
  writeLines(jsonlite::toJSON(obj, auto_unbox = TRUE, pretty = TRUE, null = "null"), tmp)
  if (file.exists(path)) unlink(path)
  file.rename(tmp, path)
  invisible(path)
}

#' 记录单个步骤的执行结果到 results/state.json
record_step <- function(cfg, id, status, seconds = NA_real_, message = "", required = TRUE) {
  path <- file.path(cfg$output$results_dir, "state.json")
  state <- if (file.exists(path)) {
    tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE), error = function(e) list(steps = list()))
  } else {
    list(steps = list())
  }
  state$steps <- Filter(function(s) !identical(s$id, id), state$steps %||% list())
  state$steps[[length(state$steps) + 1L]] <- list(
    id = id, status = status, required = required,
    seconds = round(seconds, 1), message = message
  )
  failed <- Filter(function(s) identical(s$status, "failed") && isTRUE(s$required), state$steps)
  state$required_failed <- vapply(failed, function(s) s$id, character(1))
  state$optional_failed <- vapply(
    Filter(function(s) identical(s$status, "failed") && !isTRUE(s$required), state$steps),
    function(s) s$id, character(1)
  )
  state$updated <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  write_json(path, state)
  invisible(state)
}

# ---- 绘图 ------------------------------------------------------------------

#' 把绘图表达式同时写进 PDF 和 PNG，保证设备一定关闭
#'
#' **为什么要同时出 PNG：** 产物是打包成 artifact zip 下载的，PDF 在 zip 里
#' 不能直接预览 —— 拿到 artifact 的人得先解压再找 PDF 阅读器。PNG 可以直接看。
#' PDF 保留是因为它是矢量图，放大不失真；PNG 是为了能一眼看到。
#'
#' **实现要点：** 绘图代码要跑两遍（一次 PDF、一次 PNG），所以必须用
#' `substitute()` 抓住**未求值**的表达式再 `eval()` 两次。
#' 不能写成 `force(expr); force(expr)` —— R 的 promise 有记忆，
#' 第二次 force 直接返回缓存值，**PNG 设备会开了又关、什么都不画**，得到一张空白图。
#'
#' 每次绘图都用一个独立的 `render()` 开关设备，`on.exit` 才精确对应这一次开设备；
#' 若在 `save_pdf` 主体里 `on.exit(add = TRUE)` 两次，退出时会多关一次设备，
#' 可能把调用方的设备一起关掉。
#'
#' PNG 走 `ragg`（若装了）或 `png(type="cairo")`，两者抗锯齿都更好；
#' 都没有时退回默认设备，仍然出图，不因为画质问题中断流程。
#'
#' @param path  PDF 输出路径；同名 `.png` 会写在旁边
#' @param expr  绘图表达式，会在调用者的环境里求值两次
save_pdf <- function(path, expr, width = 8, height = 6, dpi = 150) {
  code <- substitute(expr)
  env  <- parent.frame()

  render <- function(open_dev) {
    open_dev()
    on.exit(grDevices::dev.off(), add = TRUE)
    eval(code, envir = env)
  }

  render(function() grDevices::pdf(path, width = width, height = height))

  png_path <- sub("\\.pdf$", ".png", path)
  tryCatch(
    render(function() {
      if (requireNamespace("ragg", quietly = TRUE)) {
        ragg::agg_png(png_path, width = width, height = height, units = "in", res = dpi)
      } else if (capabilities("cairo")) {
        grDevices::png(png_path, width = width, height = height, units = "in",
                       res = dpi, type = "cairo")
      } else {
        grDevices::png(png_path, width = width * dpi, height = height * dpi, res = dpi)
      }
    }),
    error = function(e) {
      # 出图失败不能中断分析：PDF 已经拿到了，PNG 只是方便预览
      log_warn(sprintf("PNG 输出失败（PDF 已生成）: %s", conditionMessage(e)))
    }
  )
  invisible(path)
}

# ---- GEO SOFT 抓取（纯 base R，不需要 Bioconductor）------------------------

#' 从 GEO SOFT 文本接口抓取元数据
#'
#' 用 `form=text&view=brief` 拿到的是纯文本，比下载完整 SOFT family 小几个数量级，
#' 而且不需要 GEOquery —— 这样前置校验可以在装任何 Bioconductor 包之前跑完。
#'
#' @param acc GSE 或 GSM 编号
#' @param targ "self" 取 series 级，"gsm" 取全部样本级
#' @return 字符向量，每行一条 SOFT 记录
fetch_geo_soft <- function(acc, targ = "self", view = "brief", retries = 3L) {
  url <- sprintf(
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=%s&form=text&view=%s",
    acc, targ, view
  )
  last_err <- NULL
  for (attempt in seq_len(retries)) {
    con <- NULL
    res <- tryCatch({
      con <- url(url, open = "rt")
      lines <- readLines(con, warn = FALSE)
      if (length(lines) == 0L) stop("返回为空")
      lines
    }, error = function(e) {
      last_err <<- conditionMessage(e)
      NULL
    }, finally = {
      if (!is.null(con)) try(close(con), silent = TRUE)
    })
    if (!is.null(res)) return(res)
    log_warn(sprintf("GEO 抓取失败 (第 %d/%d 次): %s", attempt, retries, last_err))
    Sys.sleep(2 * attempt)
  }
  stop(sprintf("无法从 GEO 获取 %s 的元数据: %s", acc, last_err))
}

#' 从 SOFT 行中提取某个键的全部取值
#' @param lines SOFT 行
#' @param key 不含前导 "!" 的键名，例如 "Sample_characteristics_ch1"
soft_values <- function(lines, key) {
  pattern <- paste0("^!", key, " = ")
  hits <- grep(pattern, lines, value = TRUE)
  sub(pattern, "", hits)
}

#' 提取 SOFT 中第一个匹配键的值
soft_value <- function(lines, key, default = NA_character_) {
  v <- soft_values(lines, key)
  if (length(v) == 0L) default else v[[1L]]
}

# ---- 基因/统计小工具 -------------------------------------------------------

#' 按行 Z-score（用于热图）
row_zscore <- function(m) {
  t(scale(t(as.matrix(m))))
}

#' 为下游（富集、PPI）挑选差异基因
#'
#' **为什么需要降级路径：** spec 要求样本数 < 10，而在这个量级上，
#' 对全基因组（约 1.6 万个基因）做 BH 校正几乎不可能有任何基因通过 ——
#' GSE64790（n=6，3 对配对）实测：1,456 个基因 raw P < 0.05 且 |log2FC| > 1，
#' 但最小的 adj.P 也有 0.394。这是样本量本身的限制，不是分析错误。
#'
#' 所以：FDR 显著基因够用时用 FDR；不够时退回「raw P 排序前 N 个（仍要求
#' |log2FC| > 阈值）」。**降级必须被标注**，mode 会写进
#' enrichment_status.json / ppi_status.json，结论里不得把它当成显著差异基因。
select_degs <- function(deg, cfg, min_genes = 5L) {
  padj <- cfg$thresholds$adj_p
  lfc  <- cfg$thresholds$log2fc
  # 样本数只用于把降级原因写清楚；从 03 步的摘要里取，取不到就算了
  sp <- file.path(cfg$output$results_dir, "deg_summary.json")
  n_s <- if (file.exists(sp)) {
    tryCatch(jsonlite::fromJSON(sp)$n_samples, error = function(e) NA_integer_)
  } else NA_integer_
  if (is.null(n_s) || length(n_s) == 0L) n_s <- NA_integer_

  sig <- deg[!is.na(deg$adj.P.Val) & deg$adj.P.Val < padj & abs(deg$logFC) > lfc, , drop = FALSE]
  if (nrow(sig) >= min_genes) {
    return(list(
      genes = unique(sig$gene), mode = "fdr", n = nrow(sig), table = sig,
      reason = sprintf("adj.P < %g 且 |log2FC| > %g，共 %d 个", padj, lfc, nrow(sig))
    ))
  }

  top_n <- cfg$analysis$ranked_fallback_genes
  if (is.null(top_n)) top_n <- 500L
  cand <- deg[!is.na(deg$P.Value) & abs(deg$logFC) > lfc, , drop = FALSE]
  cand <- cand[order(cand$P.Value), , drop = FALSE]
  cand <- utils::head(cand, top_n)
  list(
    genes = unique(cand$gene), mode = "ranked_fallback", n = nrow(cand), table = cand,
    reason = sprintf(paste0(
      "FDR 显著基因仅 %d 个（需 >= %d）。样本数 %s 下对 %d 个基因做 BH 校正过严，",
      "最小的 adj.P 为 %.3f。退回按 raw P 排序、|log2FC| > %g 的前 %d 个基因。",
      "**这是假设生成，不是显著差异基因清单。**"),
      nrow(sig), min_genes, if (is.na(n_s)) "很少" else as.character(n_s),
      nrow(deg), if (nrow(deg)) min(deg$adj.P.Val, na.rm = TRUE) else NA_real_,
      lfc, nrow(cand))
  )
}

#' 安全地调用一个可选包；缺失时返回 NULL 而不是报错
require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    log_warn(sprintf("可选包 %s 未安装，跳过相关分析", pkg))
    return(FALSE)
  }
  TRUE
}
