# ============================================================================
# 01_download_clean.R — 下载表达矩阵并清洗
# ============================================================================
# spec 的 fetch_geo + clean_expression 步骤。
#
# 流程：GEOquery 下载 -> 探针映射到基因 symbol -> 去全零 -> KNN 填补 ->
#       quantile 标准化 -> 与 group.csv 对齐
#
# 输出：data/expr_raw.rds（标准化前，供 QC 对比）
#       data/expr_clean.rds / data/expr_clean.csv（标准化后）
#       data/clean_stats.json
# ============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(limma)
})

local({
  if (exists("load_config", mode = "function")) return(invisible(NULL))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  here <- if (length(file_arg) > 0L) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
  } else {
    getwd()
  }
  cand <- c(file.path(here, "lib", "common.R"),
            file.path(here, "scripts", "lib", "common.R"),
            file.path(here, "..", "lib", "common.R"))
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L) stop("找不到 lib/common.R；请从仓库根目录运行")
  source(hit[[1L]])
})

# 不同平台注释里基因 symbol 的列名差异很大，按优先级依次尝试
SYMBOL_COLUMNS <- c("GENE_SYMBOL", "Gene Symbol", "GeneSymbol", "GENE", "Symbol",
                    "gene_symbol", "symbol", "GENE_NAME", "ILMN_Gene")

#' 从平台注释中挑出基因 symbol 列
pick_symbol_column <- function(fdata) {
  for (col in SYMBOL_COLUMNS) {
    if (col %in% colnames(fdata)) {
      v <- as.character(fdata[[col]])
      if (sum(nzchar(v) & !is.na(v)) > 100L) {
        log_info(sprintf("使用平台注释列 '%s' 作为基因 symbol 来源", col))
        return(col)
      }
    }
  }
  NULL
}

#' 同一 symbol 的多探针取方差最大者
collapse_to_symbol <- function(expr, symbols) {
  keep <- !is.na(symbols) & nzchar(symbols) & symbols != "---"
  expr <- expr[keep, , drop = FALSE]
  symbols <- symbols[keep]
  if (nrow(expr) == 0L) stop("探针映射后没有任何基因，请检查平台注释")

  vars <- matrixStats_rownanvar(expr)
  ord <- order(symbols, -vars)
  expr <- expr[ord, , drop = FALSE]
  symbols <- symbols[ord]
  dup <- duplicated(symbols)
  log_info(sprintf("多探针折叠: %d 个探针 -> %d 个基因（取方差最大探针）",
                   length(symbols), length(unique(symbols))))
  expr <- expr[!dup, , drop = FALSE]
  rownames(expr) <- symbols[!dup]
  expr
}

#' base R 的按行方差（忽略 NA），避免额外依赖 matrixStats
matrixStats_rownanvar <- function(m) {
  apply(m, 1L, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2L) 0 else stats::var(x)
  })
}

run_01_download_clean <- function(cfg) {
  log_info("=== 步骤 01：下载与清洗 ===")
  ensure_dirs(cfg)

  group <- utils::read.csv(file.path(cfg$output$data_dir, "group.csv"),
                           stringsAsFactors = FALSE)
  if (nrow(group) == 0L) stop("group.csv 为空，请先运行 00_validate_inputs.R")

  cache_dir <- file.path(cfg$output$data_dir, "geo_cache")
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- 1. 下载 ------------------------------------------------------------
  log_info(sprintf("从 GEO 下载 %s（含平台注释）...", cfg$dataset_id))
  eset <- tryCatch(
    GEOquery::getGEO(cfg$dataset_id, GSEMatrix = TRUE, getGPL = TRUE,
                     destdir = cache_dir, AnnotGPL = FALSE),
    error = function(e) stop(sprintf("GEO 下载失败: %s", conditionMessage(e)))
  )
  if (is.list(eset)) {
    if (length(eset) > 1L) {
      log_warn(sprintf("该 series 含 %d 个平台，取第一个", length(eset)))
    }
    eset <- eset[[1L]]
  }
  log_info(sprintf("下载完成: %d 个探针 x %d 个样本", nrow(eset), ncol(eset)))

  # ---- 2. 取表达矩阵与样本 ID ---------------------------------------------
  expr <- Biobase::exprs(eset)
  pd <- Biobase::pData(eset)
  gsm <- if ("geo_accession" %in% colnames(pd)) as.character(pd$geo_accession) else rownames(pd)
  if (length(gsm) != ncol(expr)) stop("样本数与表达矩阵列数不一致")
  colnames(expr) <- gsm

  # 双色 Agilent 芯片返回的是 log2 ratio；单色芯片若明显未取对数则补取
  if (stats::median(expr, na.rm = TRUE) > 50) {
    log_warn("表达值中位数 > 50，判定为未取对数，执行 log2(x + 1)")
    expr[expr < 0] <- NA
    expr <- log2(expr + 1)
  }

  # ---- 3. 探针 -> 基因 symbol ---------------------------------------------
  fdata <- Biobase::fData(eset)
  sym_col <- pick_symbol_column(fdata)
  if (is.null(sym_col)) {
    stop(sprintf(
      "平台 %s 的注释中没有可用的基因 symbol 列（已尝试: %s）。\n  可用列: %s\n  GO/KEGG 富集需要 symbol，无法继续。",
      Biobase::annotation(eset), paste(SYMBOL_COLUMNS, collapse = ", "),
      paste(colnames(fdata), collapse = ", ")))
  }
  expr <- collapse_to_symbol(expr, as.character(fdata[[sym_col]]))

  # ---- 4. 去全零 / 全 NA 基因 ---------------------------------------------
  before <- nrow(expr)
  keep <- rowSums(!is.na(expr)) >= 2L & apply(expr, 1L, function(x) {
    v <- x[!is.na(x)]
    length(v) > 0L && stats::sd(v) > 0
  })
  expr <- expr[keep, , drop = FALSE]
  log_info(sprintf("过滤无变异基因: %d -> %d", before, nrow(expr)))
  if (nrow(expr) < 100L) stop("过滤后剩余基因过少，数据可能有问题")

  # ---- 5. KNN 填补 --------------------------------------------------------
  n_missing <- sum(is.na(expr))
  if (n_missing > 0L) {
    if (!requireNamespace("impute", quietly = TRUE)) {
      stop(sprintf("存在 %d 个缺失值但 impute 包不可用，无法按 spec 做 KNN 填补", n_missing))
    }
    log_info(sprintf("KNN 填补 %d 个缺失值 (%.2f%%)，k=%d",
                     n_missing, 100 * n_missing / length(expr), cfg$analysis$impute_k))
    expr <- impute::impute.knn(expr, k = cfg$analysis$impute_k)$data
  } else {
    log_info("无缺失值，跳过 KNN 填补")
  }

  # 标准化前矩阵留档，供 QC 画 before/after
  expr_raw <- expr

  # ---- 6. quantile 标准化 -------------------------------------------------
  expr <- limma::normalizeBetweenArrays(expr, method = "quantile")
  log_info("已执行 quantile 标准化")

  # ---- 7. 与分组对齐（spec: expr columns match meta rownames）------------
  missing_samples <- setdiff(group$gsm, colnames(expr))
  if (length(missing_samples) > 0L) {
    stop(sprintf("表达矩阵缺少分组文件中的样本: %s", paste(missing_samples, collapse = ", ")))
  }
  extra <- setdiff(colnames(expr), group$gsm)
  if (length(extra) > 0L) {
    log_warn(sprintf("表达矩阵含分组文件之外的样本，已丢弃: %s", paste(extra, collapse = ", ")))
  }
  expr <- expr[, group$gsm, drop = FALSE]
  expr_raw <- expr_raw[, group$gsm, drop = FALSE]
  stopifnot(identical(colnames(expr), group$gsm))
  log_info(sprintf("样本对齐完成: %d 个样本按 group.csv 顺序排列", ncol(expr)))

  # ---- 8. 落盘 ------------------------------------------------------------
  saveRDS(expr_raw, file.path(cfg$output$data_dir, "expr_raw.rds"))
  saveRDS(expr, file.path(cfg$output$data_dir, "expr_clean.rds"))
  utils::write.csv(data.frame(gene = rownames(expr), expr, check.names = FALSE),
                   file.path(cfg$output$data_dir, "expr_clean.csv"), row.names = FALSE)
  write_json(file.path(cfg$output$data_dir, "clean_stats.json"), list(
    platform = Biobase::annotation(eset),
    symbol_column = sym_col,
    genes_final = nrow(expr),
    samples = ncol(expr),
    missing_imputed = n_missing,
    normalization = "quantile (limma::normalizeBetweenArrays)",
    probe_collapse = "max variance per symbol"
  ))

  log_info(sprintf("已写出 expr_raw.rds / expr_clean.rds / expr_clean.csv（%d 基因 x %d 样本）",
                   nrow(expr), ncol(expr)))
  invisible(expr)
}

if (!GEO_ORCHESTRATED()) {
  cfg <- load_config()
  run_01_download_clean(cfg)
}
