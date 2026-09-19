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
                    "gene_symbol", "symbol", "ILMN_Gene")
# 没有 symbol 列时的替代途径
ACCESSION_COLUMNS <- c("GB_ACC", "GB_ACCESSION", "ACCESSION", "GenBank", "GB_LIST", "REFSEQ")
GENENAME_COLUMNS  <- c("DESCRIPTION", "Gene Title", "GENE_NAME", "gene_assignment_name")

# 判定「symbol 模式是否可用」的两个条件，满足**任一**即可：
#   1) 探针覆盖率 >= MIN_SYMBOL_COVERAGE
#   2) 唯一 symbol 数 >= MIN_SYMBOL_GENES
#
# 为什么需要条件 2：覆盖率是「带注释的探针 / 全部探针」，对 lncRNA 芯片、
# 外显子芯片这类平台天然偏低 —— 大部分探针本来就对应非编码转录本，没有 symbol
# 是正常的。GSE64790（Agilent lncRNA V4.0）覆盖率只有 35%，但仍有 21,812 个
# 带 symbol 的探针，做 GO/KEGG/STRING 绰绰有余。只看覆盖率会把它误判成探针模式。
MIN_SYMBOL_COVERAGE <- 0.5
MIN_SYMBOL_GENES    <- 5000L

#' 从平台注释中挑出基因 symbol 列
pick_symbol_column <- function(fdata) {
  for (col in SYMBOL_COLUMNS) {
    if (col %in% colnames(fdata)) {
      v <- as.character(fdata[[col]])
      if (sum(nzchar(v) & !is.na(v)) > 100L) return(col)
    }
  }
  NULL
}

#' 在 fdata 中按优先级找第一个可用列，返回其字符向量
first_column <- function(fdata, candidates) {
  for (col in candidates) {
    if (col %in% colnames(fdata)) {
      v <- as.character(fdata[[col]])
      if (sum(nzchar(v) & !is.na(v) & v != "---") > 100L) {
        return(list(column = col, values = v))
      }
    }
  }
  NULL
}

#' 用 org.Hs.eg.db 把一批 key 映射到 SYMBOL
#' @return 命名向量 key -> symbol（只含成功映射的）
lookup_symbols <- function(keys, keytype) {
  keys <- unique(keys[!is.na(keys) & nzchar(keys)])
  if (length(keys) == 0L) return(stats::setNames(character(0), character(0)))
  hit <- suppressWarnings(tryCatch(
    AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = keys,
                          keytype = keytype, columns = "SYMBOL"),
    error = function(e) {
      log_warn(sprintf("%s 映射失败: %s", keytype, conditionMessage(e)))
      NULL
    }
  ))
  if (is.null(hit) || nrow(hit) == 0L) return(stats::setNames(character(0), character(0)))
  hit <- hit[!is.na(hit$SYMBOL), , drop = FALSE]
  stats::setNames(as.character(hit$SYMBOL), as.character(hit[[keytype]]))
}

#' 基因名归一化：去括号内容、去所有非字母数字、转小写
#'
#' GPL16025 的 DESCRIPTION 是 2007 年前后的旧基因名，与今天的 GENENAME 经常只差
#' 标点或一个括号补充。例如：
#'   "SH3-domain binding protein 2"                      vs "SH3 domain binding protein 2"
#'   "Rap guanine nucleotide exchange factor (GEF) 2"    vs "Rap guanine nucleotide exchange factor 2"
#'   "solute carrier family 15 (oligopeptide transporter), member 1"
#'                                                       vs "solute carrier family 15 member 1"
#' 归一化后这三组都能对上，精确匹配则全部落空。
normalize_gene_name <- function(x) {
  x <- tolower(x)
  x <- gsub("\\([^)]*\\)", " ", x)
  gsub("[^a-z0-9]+", "", x)
}

.gene_name_cache <- new.env(parent = emptyenv())

#' 构建 GENENAME -> SYMBOL 映射（精确 + 归一化），进程内只算一次
genename_map <- function() {
  if (!is.null(.gene_name_cache$map)) return(.gene_name_cache$map)
  m <- tryCatch({
    db <- org.Hs.eg.db::org.Hs.eg.db
    gn_keys <- AnnotationDbi::keys(db, keytype = "GENENAME")
    hit <- suppressWarnings(AnnotationDbi::select(db, keys = gn_keys, keytype = "GENENAME",
                                                  columns = "SYMBOL"))
    hit <- hit[!is.na(hit$SYMBOL) & !is.na(hit$GENENAME), , drop = FALSE]
    exact <- stats::setNames(as.character(hit$SYMBOL), as.character(hit$GENENAME))
    exact <- exact[!duplicated(names(exact))]

    nk <- normalize_gene_name(names(exact))
    keep <- nzchar(nk)
    norm <- stats::setNames(unname(exact)[keep], nk[keep])
    norm <- norm[!duplicated(names(norm))]

    log_info(sprintf("GENENAME 映射表: %d 条精确 + %d 条归一化", length(exact), length(norm)))
    list(exact = exact, norm = norm)
  }, error = function(e) {
    log_warn(sprintf("构建 GENENAME 映射表失败: %s", conditionMessage(e)))
    list(exact = stats::setNames(character(0), character(0)),
         norm  = stats::setNames(character(0), character(0)))
  })
  .gene_name_cache$map <- m
  m
}

#' 三级探针 -> 基因 symbol 映射
#'
#' GPL16025（NimbleGen）这类平台只有 ID / GB_ACC / DESCRIPTION 三列，没有 symbol，
#' 也没有 GEO curated 注释（`GPL16025.annot.gz` 返回 404）。所以除了直接读 symbol 列，
#' 还要能走 GenBank accession（ACCNUM）、RefSeq（REFSEQ）和基因全名（GENENAME）。
#' 全部途径都拿不到足够覆盖率时返回 probe 模式 —— 基因层面分析照常做，但下游必须
#' 跳过 GO/KEGG 并说明原因，而不是拿探针 ID 冒充基因去富集。
map_features_to_symbols <- function(ids, fdata) {
  n <- length(ids)
  coverage <- function(v) sum(!is.na(v) & nzchar(v)) / n
  results <- list()
  add <- function(symbols, method) {
    results[[length(results) + 1L]] <<- list(symbols = symbols, method = method,
                                             coverage = coverage(symbols))
  }
  blank <- function() rep(NA_character_, n)

  # ---- 途径 1：平台注释自带 symbol 列 ------------------------------------
  col <- pick_symbol_column(fdata)
  if (!is.null(col)) {
    s <- as.character(fdata[[col]])
    s[!nzchar(s) | s == "---"] <- NA
    add(s, sprintf("platform annotation column '%s'", col))
  }

  has_db <- requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
            requireNamespace("AnnotationDbi", quietly = TRUE)
  if (has_db && (is.null(col) || coverage(results[[1L]]$symbols) < MIN_SYMBOL_COVERAGE)) {

    # ---- 途径 2/3：accession -> ACCNUM / REFSEQ --------------------------
    acc <- first_column(fdata, ACCESSION_COLUMNS)
    if (!is.null(acc)) {
      # org.Hs.eg.db 的 ACCNUM 不带版本号后缀
      keys <- sub("\\.[0-9]+$", "", trimws(acc$values))
      keys[!nzchar(keys) | keys == "---"] <- NA
      add(unname(lookup_symbols(keys, "ACCNUM")[keys]),
          sprintf("accession '%s' -> ACCNUM", acc$column))

      is_refseq <- !is.na(keys) & grepl("^(NM_|NR_|XM_|XR_)", keys)
      if (any(is_refseq)) {
        s <- blank()
        s[is_refseq] <- unname(lookup_symbols(keys[is_refseq], "REFSEQ")[keys[is_refseq]])
        add(s, sprintf("RefSeq subset of '%s' -> REFSEQ", acc$column))
      }
    }

    # ---- 途径 4/5：基因全名 -> GENENAME（精确 / 归一化）------------------
    gn <- first_column(fdata, GENENAME_COLUMNS)
    if (!is.null(gn)) {
      keys <- trimws(gn$values)
      keys[!nzchar(keys) | keys == "---"] <- NA
      gm <- genename_map()
      add(unname(gm$exact[keys]),
          sprintf("gene name '%s' -> GENENAME (exact)", gn$column))
      add(unname(gm$norm[normalize_gene_name(keys)]),
          sprintf("gene name '%s' -> GENENAME (normalized)", gn$column))
    }
  }

  if (length(results) == 0L) {
    return(list(mode = "probe", symbols = NULL, method = "none", coverage = 0,
                reason = "平台注释中没有 symbol 列，也没有可用的 accession / 基因全名列"))
  }

  for (r in results) {
    syms <- r$symbols[!is.na(r$symbols) & nzchar(r$symbols)]
    log_info(sprintf("  映射途径 %-52s 覆盖率 %5.1f%%  唯一 symbol %d",
                     r$method, 100 * r$coverage, length(unique(syms))))
  }

  best <- results[[which.max(vapply(results, function(r) r$coverage, numeric(1)))]]
  best_syms <- best$symbols[!is.na(best$symbols) & nzchar(best$symbols)]
  n_genes <- length(unique(best_syms))

  if (best$coverage < MIN_SYMBOL_COVERAGE && n_genes < MIN_SYMBOL_GENES) {
    return(list(mode = "probe", symbols = NULL, method = best$method, coverage = best$coverage,
                mapped_genes = n_genes,
                reason = sprintf(paste0("最佳映射途径 '%s' 覆盖率 %.1f%%（阈值 %.0f%%）",
                                        "且唯一 symbol 仅 %d 个（阈值 %d），不足以做基因层面分析"),
                                 best$method, 100 * best$coverage, 100 * MIN_SYMBOL_COVERAGE,
                                 n_genes, MIN_SYMBOL_GENES)))
  }
  list(mode = "symbol", symbols = best$symbols, method = best$method,
       coverage = best$coverage, mapped_genes = n_genes, reason = NA_character_)
}

#' 抓取平台注释表
#'
#' **不要用 GEOquery 的 getGPL=TRUE。** 对 GPL16025，那条路会下载 182 MB 的
#' `GPL16025_family.soft.gz`（里面含该平台上千个 GSM 的完整记录），解析又慢又吃内存。
#' 而 GEO 的 CGI `view=full` 返回**同样完整的 45,033 行**注释表，只有 2.6 MB。
#'
#' 结果缓存为 RDS，重复运行不再联网。
fetch_platform_annotation <- function(gpl_id, cache_dir) {
  cache <- file.path(cache_dir, sprintf("%s_annotation.rds", gpl_id))
  if (file.exists(cache)) {
    log_info(sprintf("使用缓存的平台注释: %s", cache))
    return(readRDS(cache))
  }

  log_info(sprintf("抓取 %s 注释表（CGI view=full，约 2.6 MB）...", gpl_id))
  lines <- fetch_geo_soft(gpl_id, targ = "self", view = "full")
  begin <- grep("^!platform_table_begin", lines)
  end   <- grep("^!platform_table_end", lines)
  if (length(begin) == 0L || length(end) == 0L) {
    stop(sprintf("%s 的 SOFT 中没有平台注释表", gpl_id))
  }
  header <- strsplit(lines[begin[1L] + 1L], "\t", fixed = TRUE)[[1L]]
  body <- lines[(begin[1L] + 2L):(end[1L] - 1L)]
  fields <- strsplit(body, "\t", fixed = TRUE)

  # 少数行字段数不足，补 NA 保证能拼成矩形
  ncol_expected <- length(header)
  fields <- lapply(fields, function(f) {
    length(f) <- ncol_expected
    f
  })
  df <- as.data.frame(do.call(rbind, fields), stringsAsFactors = FALSE)
  colnames(df) <- header
  df[] <- lapply(df, function(x) { x[is.na(x)] <- ""; trimws(x) })

  log_info(sprintf("平台注释: %d 行 x %d 列 (%s)", nrow(df), ncol(df),
                   paste(header, collapse = ", ")))
  saveRDS(df, cache)
  df
}

#' 同一 symbol 的多探针取方差最大者
collapse_to_symbol <- function(expr, symbols) {
  keep <- !is.na(symbols) & nzchar(symbols) & symbols != "---"
  expr <- expr[keep, , drop = FALSE]
  symbols <- symbols[keep]
  if (nrow(expr) == 0L) stop("探针映射后没有任何基因，请检查平台注释")

  vars <- row_variance(expr)
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
row_variance <- function(m) {
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
  # getGPL=FALSE：平台注释单独用轻量接口取，见 fetch_platform_annotation()
  log_info(sprintf("从 GEO 下载 %s 表达矩阵...", cfg$dataset_id))
  eset <- tryCatch(
    GEOquery::getGEO(cfg$dataset_id, GSEMatrix = TRUE, getGPL = FALSE,
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

  # 单色芯片返回的已是 log2 强度；若明显未取对数（中位数 > 50）则补取
  if (stats::median(expr, na.rm = TRUE) > 50) {
    log_warn("表达值中位数 > 50，判定为未取对数，执行 log2(x + 1)")
    expr[expr < 0] <- NA
    expr <- log2(expr + 1)
  }

  # ---- 3. 平台注释 -> 探针映射 -------------------------------------------
  gpl_id <- Biobase::annotation(eset)
  if (is.null(gpl_id) || !nzchar(gpl_id) || identical(gpl_id, "NA")) {
    gpl_id <- unique(as.character(Biobase::pData(eset)$platform_id))[1L]
  }
  log_info(sprintf("平台: %s", gpl_id))

  fdata <- tryCatch(
    fetch_platform_annotation(gpl_id, cache_dir),
    error = function(e) {
      log_warn(sprintf("轻量注释抓取失败，回退到 GEOquery getGPL=TRUE（会下载大文件）: %s",
                       conditionMessage(e)))
      eset <<- GEOquery::getGEO(cfg$dataset_id, GSEMatrix = TRUE, getGPL = TRUE,
                                destdir = cache_dir, AnnotGPL = FALSE)
      if (is.list(eset)) eset <- eset[[1L]]
      Biobase::fData(eset)
    }
  )

  # 把注释对齐到表达矩阵的探针顺序
  if ("ID" %in% colnames(fdata)) {
    idx <- match(rownames(expr), as.character(fdata$ID))
    hit <- sum(!is.na(idx))
    log_info(sprintf("注释对齐: %d/%d 个探针在平台注释中找到", hit, nrow(expr)))
    if (hit < 0.5 * nrow(expr)) {
      log_warn("超过一半探针在平台注释中找不到，映射结果可能不可靠")
    }
    fdata <- fdata[idx, , drop = FALSE]
    rownames(fdata) <- rownames(expr)
  } else {
    log_warn(sprintf("平台注释没有 ID 列（实际列: %s），按行号对齐",
                     paste(colnames(fdata), collapse = ", ")))
  }

  mapping <- map_features_to_symbols(rownames(expr), fdata)
  log_info(sprintf("特征映射途径: %s（覆盖率 %.1f%%）", mapping$method, 100 * mapping$coverage))

  if (identical(mapping$mode, "symbol")) {
    expr <- collapse_to_symbol(expr, mapping$symbols)
    feature_ids <- rownames(expr)
  } else {
    # 拿不到 symbol 就退回探针层面：QC/PCA/相关性/DEG/热图仍然成立，
    # 但下游 04b 必须跳过 GO/KEGG 并说明原因，不能拿探针 ID 冒充基因去富集。
    log_warn(sprintf("无法映射到基因 symbol（%s）", mapping$reason))
    log_warn("退回探针层面分析；GO/KEGG 富集将被跳过并在 enrichment_status.json 中说明原因")
    feature_ids <- rownames(expr)
  }

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
    platform_annotation_columns = colnames(fdata),
    feature_mode = mapping$mode,
    feature_id_type = if (identical(mapping$mode, "symbol")) "gene symbol" else "platform probe ID",
    mapping_method = mapping$method,
    mapping_coverage = round(mapping$coverage, 4),
    mapping_unique_symbols = mapping$mapped_genes,
    mapping_reason = mapping$reason,
    genes_final = nrow(expr),
    samples = ncol(expr),
    missing_imputed = n_missing,
    normalization = "quantile (limma::normalizeBetweenArrays)",
    probe_collapse = if (identical(mapping$mode, "symbol")) "max variance per symbol" else NA
  ))

  # 供 04b 判断能否做富集
  write_json(file.path(cfg$output$data_dir, "feature_mode.json"), list(
    mode = mapping$mode,
    method = mapping$method,
    coverage = round(mapping$coverage, 4),
    reason = mapping$reason
  ))

  log_info(sprintf("已写出 expr_raw.rds / expr_clean.rds / expr_clean.csv（%d 个特征 x %d 样本）",
                   nrow(expr), ncol(expr)))
  invisible(expr)
}

if (!GEO_ORCHESTRATED()) {
  cfg <- load_config()
  run_01_download_clean(cfg)
}
