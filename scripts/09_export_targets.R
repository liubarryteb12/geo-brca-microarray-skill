# ============================================================================
# 09_export_targets.R — 导出 Part 2 用的候选靶基因（§1.7 / §1.8 交接）
#
# ## 这一步存在的唯一理由
#
# 规范把 §1.7 / §1.8（虚拟敲除 / 过表达）标为**保留框架**、主语言 Python，
# 而候选靶基因由 Part 1 产出。§0.2 规定跨部分交接**只走 CSV**（不走 RDS）——
# 所以 Part 1 的最后一个动作是把候选基因整理成一张 CSV 交出去。
#
# **不交接的后果不是"少一步"，是"下游自己挑基因"。**
# Part 2 的 `08_virtual_perturbation.py` 在没有这张表时会回退到它自己的
# 调控子，那时虚拟扰动的靶基因就不再来自 Part 1 的预后/网络证据 ——
# 而下游状态文件里会写着 `internal_top_regulons`，那是**可以核对**的，
# 比默默换了候选集好。
#
# ## 候选来源与优先级
#
#   1. lasso_signature   lasso_coefficients.csv   —— 预后签名（§1.6）
#   2. wgcna_hub         wgcna_modules.csv        —— 模块成员（§1.3）
#   3. ppi_hub           hub_genes.csv            —— PPI 枢纽（§1.5）
#   4. tf_regulon        tf_regulon_enrichment.csv—— 显著调控子（§1.5）
#
# **同名的基因保留优先级最高的那个来源，但 `n_sources` 记下它被几个来源
# 提到过。** 一个基因被两条独立证据同时点名，和只被一条点名，
# 是很不同的证据强度 —— 去重时把这个信息丢掉，等于把最有价值的信号抹平。
#
# ## 产物
#
#   results/<GSE>/part2_targets.csv        交给 Part 2 的候选靶基因
#   results/<GSE>/part2_targets_status.json 各来源计数、缺失来源、方法学限定
#
# **`logFC` 一并带上**：Part 2 的虚拟扰动要算 `signature_alignment`
# （预测变化方向 vs 疾病签名方向），没有 logFC 那一列只能是空的。
# ============================================================================

suppressPackageStartupMessages({
  library(stats)
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

# 交接文件名。Part 2 的 config 里 perturbation.targets_csv 指向它。
TARGETS_FILE <- "part2_targets.csv"

# 每个来源最多取多少个基因。**必须有上限** —— WGCNA 的模块可能上千个基因，
# 全交出去等于没筛，下游的虚拟扰动会跑成"对全转录组做扰动"。
MAX_PER_SOURCE <- 200L


#' 在数据框里按候选名字找一个列，找不到返回 NULL
#'
#' **不写死列名。** 上游脚本的列名是 `gene` / `tf` / `logFC`，
#' 大小写还不统一；写死一个名字时上游改名不会报错，
#' 只会让这一步静默交出一张空表 —— 而空表看起来像"没有候选基因"。
pick_col <- function(df, candidates) {
  for (nm in candidates) {
    if (nm %in% colnames(df)) return(nm)
    hit <- colnames(df)[tolower(colnames(df)) == tolower(nm)]
    if (length(hit) > 0L) return(hit[[1L]])
  }
  NULL
}


#' 读一个来源，返回 data.frame(gene, source, source_rank) 或 NULL
read_source <- function(path, label, gene_cols, rank_col = NULL,
                        decreasing = TRUE) {
  if (!file.exists(path)) {
    log_warn(sprintf("  来源 %s 缺失（%s）—— 跳过", label, basename(path)))
    return(NULL)
  }
  df <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE,
                                 check.names = FALSE),
                 error = function(e) {
                   log_warn(sprintf("  来源 %s 读取失败: %s", label, conditionMessage(e)))
                   NULL
                 })
  if (is.null(df) || nrow(df) == 0L) {
    log_warn(sprintf("  来源 %s 是空表 —— 跳过", label))
    return(NULL)
  }
  gcol <- pick_col(df, gene_cols)
  if (is.null(gcol)) {
    log_warn(sprintf("  来源 %s 里找不到基因列（试过 %s，实际列：%s）—— 跳过",
                     label, paste(gene_cols, collapse = "/"),
                     paste(colnames(df), collapse = ", ")))
    return(NULL)
  }
  out <- data.frame(gene = as.character(df[[gcol]]), stringsAsFactors = FALSE)
  out <- out[!is.na(out$gene) & nzchar(out$gene), , drop = FALSE]
  if (nrow(out) == 0L) return(NULL)

  # 排序：有排序依据就按它排，否则保持原顺序（上游已经排过了）
  if (!is.null(rank_col) && rank_col %in% colnames(df)) {
    key <- suppressWarnings(as.numeric(df[[rank_col]]))
    ord <- order(if (decreasing) -abs(key) else key, na.last = TRUE)
    out <- out[ord, , drop = FALSE]
  }
  out <- out[!duplicated(out$gene), , drop = FALSE]
  n_total <- nrow(out)
  out <- utils::head(out, MAX_PER_SOURCE)
  out$source <- label
  out$source_rank <- seq_len(nrow(out))
  log_info(sprintf("  来源 %s: %d 个基因（取前 %d）", label, n_total, nrow(out)))
  out
}


run_09_export_targets <- function(cfg) {
  res <- cfg$output$results_dir
  status_path <- file.path(res, "part2_targets_status.json")

  sources <- list(
    list(file = "lasso_coefficients.csv", label = "lasso_signature",
         gene_cols = c("gene", "symbol", "tf"), rank_col = "coef"),
    list(file = "wgcna_modules.csv", label = "wgcna_hub",
         gene_cols = c("gene", "symbol"), rank_col = NULL),
    list(file = "hub_genes.csv", label = "ppi_hub",
         gene_cols = c("gene", "symbol"), rank_col = "degree"),
    list(file = "tf_regulon_enrichment.csv", label = "tf_regulon",
         gene_cols = c("tf", "gene", "symbol"), rank_col = "p_value",
         decreasing = FALSE)
  )

  log_info("收集候选靶基因（交给 Part 2 做虚拟扰动，§1.7/§1.8）")
  parts <- list()
  missing_sources <- character(0)
  for (s in sources) {
    p <- read_source(file.path(res, s$file), s$label, s$gene_cols,
                     s$rank_col, s$decreasing %||% TRUE)
    if (is.null(p)) {
      missing_sources <- c(missing_sources, s$label)
    } else {
      parts[[length(parts) + 1L]] <- p
    }
  }

  if (length(parts) == 0L) {
    status <- list(
      dataset_id = cfg$dataset_id,
      status = "no_sources",
      reason = paste0("四个候选来源一个都没有产出：",
                      paste(missing_sources, collapse = ", "),
                      " —— 上游步骤可能都没跑或都失败了"),
      sources_missing = as.list(missing_sources),
      limitations = c(
        "没有交接表时，Part 2 的虚拟扰动会回退到它自己的调控子，",
        "那时靶基因不再来自 Part 1 的预后/网络证据。"
      )
    )
    write_json(status_path, status)
    log_warn(status$reason)
    return(invisible(status))
  }

  all <- do.call(rbind, parts)

  # **一个基因被几个来源点名 —— 去重前先数。**
  # 这是这张表里信息量最大的一列：两条独立证据同时点名，
  # 和只被一条点名，是很不同的证据强度。
  n_src <- tapply(all$gene, all$gene, function(x) length(unique(all$source)))
  # 优先级 = 来源在 sources 里的次序（lasso 最高）
  prio <- vapply(all$source, function(lb) {
    which(vapply(sources, function(s) identical(s$label, lb), logical(1)))[1L]
  }, integer(1))
  all$..prio <- prio
  all <- all[order(all$..prio, all$source_rank), , drop = FALSE]
  keep <- !duplicated(all$gene)
  out <- all[keep, , drop = FALSE]
  out$n_sources <- as.integer(n_src[out$gene])
  out$..prio <- NULL
  out <- out[order(-out$n_sources, out$source_rank), , drop = FALSE]
  rownames(out) <- NULL

  # ---- 附上 logFC（方向）--------------------------------------------------
  deg_path <- file.path(res, "deg_table.csv")
  n_with_logfc <- 0L
  if (file.exists(deg_path)) {
    deg <- tryCatch(utils::read.csv(deg_path, stringsAsFactors = FALSE,
                                    check.names = FALSE),
                    error = function(e) NULL)
    if (!is.null(deg)) {
      dg <- pick_col(deg, c("gene", "symbol"))
      dfc <- pick_col(deg, c("logFC", "log2FC", "log_fc"))
      if (!is.null(dg) && !is.null(dfc)) {
        m <- suppressWarnings(as.numeric(deg[[dfc]]))
        names(m) <- as.character(deg[[dg]])
        out$logfc <- unname(m[out$gene])
        n_with_logfc <- sum(!is.na(out$logfc))
        log_info(sprintf("  附上 logFC：%d/%d 个基因匹配到（列 %s）",
                         n_with_logfc, nrow(out), dfc))
      } else {
        log_warn(sprintf("  deg_table.csv 里找不到基因列或 logFC 列（实际列：%s）",
                         paste(colnames(deg), collapse = ", ")))
      }
    }
  } else {
    log_warn("  deg_table.csv 不存在")
  }
  if (n_with_logfc == 0L) {
    out$logfc <- NA_real_
    log_warn("  **交接表没有 logFC**：Part 2 的 signature_alignment 会是空的")
  }

  # 列顺序固定，方便下游按名读
  out <- out[, c("gene", "logfc", "source", "source_rank", "n_sources")]
  utils::write.csv(out, file.path(res, TARGETS_FILE), row.names = FALSE)
  log_info(sprintf("交接表: %s（%d 个基因，来自 %d 个来源）",
                   TARGETS_FILE, nrow(out), length(parts)))

  # ---- §0.2 跨语言交接登记 ------------------------------------------------
  # **这是本仓库唯一一处真正的跨部分交接**，所以清单里必须有记录：
  # 交出去多少、什么格式、丢了什么字段。
  record_cross_language(
    cfg,
    src = "geo-normal-pipeline-skill (R, Part 1)",
    dst = "scrna-pipeline-skill (Python, Part 2)",
    format = "csv",
    tool = "utils::write.csv / pandas.read_csv",
    before = list(n_genes_all_sources = nrow(all),
                  n_sources = length(parts)),
    after = list(n_genes = nrow(out), n_with_logfc = n_with_logfc),
    lost = c("module_label (wgcna 模块归属未交接)",
             "degree (PPI 度数未交接)",
             "p_adj_bh (调控子显著性未交接)",
             "coef (LASSO 系数未交接)"),
    note = paste0("只交基因名 + logFC + 来源标签。丢失的字段在 Part 2 用不到，",
                  "但**记录下来了** —— 将来要用时知道回 Part 1 的哪张表取。")
  )

  status <- list(
    dataset_id = cfg$dataset_id,
    status = "ok",
    output_file = TARGETS_FILE,
    n_genes = nrow(out),
    n_with_logfc = n_with_logfc,
    n_sources_used = length(parts),
    sources_used = as.list(unique(all$source)),
    sources_missing = as.list(missing_sources),
    n_multi_source_genes = sum(out$n_sources > 1L),
    per_source_counts = as.list(table(all$source)),
    top_genes = utils::head(out, 15),
    # **R 没有隐式字符串拼接**（不像 Python/C）。写成
    #   method = ("第一段" "第二段")
    # 是语法错误：`unexpected string constant`。
    # 必须显式 paste0()。这个错在 CI 上炸过一次（run 35482359507）。
    method = paste0("汇总 Part 1 的四类候选来源（LASSO 签名 / WGCNA 模块 / ",
                    "PPI 枢纽 / 显著调控子），按优先级去重，附上 logFC"),
    limitations = c(
      "**去重保留了优先级最高的来源，但 n_sources 记下了它被几个来源点名** ——",
      "      被多条证据同时点名的基因证据更强，这一列不要忽略",
      "**每个来源最多取 200 个基因**（MAX_PER_SOURCE）：WGCNA 模块可能上千个，",
      "      全交出去等于没筛，下游会变成对全转录组做扰动",
      "logFC 来自 Part 1 的肿瘤 vs 正常比较，**不是预后方向的效应量** ——",
      "      用它算 signature_alignment 时要知道这一点",
      if (n_with_logfc == 0L)
        "**本轮没有 logFC**：Part 2 的 signature_alignment 会是空的"
      else
        sprintf("logFC 匹配到 %d/%d 个基因", n_with_logfc, nrow(out)),
      "候选来源本身的可靠性上界：LASSO 受 EPV 限制（见 §1.6 的过拟合警告）、",
      "      WGCNA 是共表达模块不是因果、PPI 枢纽依赖 STRING 覆盖率"
    ),
    outputs = c(TARGETS_FILE)
  )
  write_json(status_path, status)
  log_info("候选靶基因交接表导出完成")
  invisible(status)
}


if (!GEO_ORCHESTRATED()) {
  cfg <- parse_args()
  run_09_export_targets(cfg)
}
