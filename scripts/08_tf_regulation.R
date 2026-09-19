# ============================================================================
# 08_tf_regulation.R — 转录因子调控分析
# ============================================================================
# 用户文档把 TF 调控列为 bulk 部分的分析项。本脚本补齐它。
#
# ## 两个分析，回答两个不同的问题
#
#   1. **调控子富集**（regulon enrichment，Fisher 精确检验）
#      问：某个 TF 的靶基因在差异基因里是否**过度出现**？
#      这是主分析 —— 它不依赖样本量，n<10 时仍然成立。
#
#   2. **调控子活性分数**（regulon activity score，按样本）
#      问：某个 TF 的靶基因整体在两组之间是否**系统性地位移**？
#      分数 = 靶基因 z-score 的加权均值（mode of regulation 定符号）。
#      **n<10 时这个分数噪声很大**，所以它只作为佐证，
#      并在产物里写明样本量警告。
#
# ## 为什么不用 decoupleR 的 ULM
#
# decoupleR 是更标准的做法（单样本线性模型 + 置换检验），但它多一个
# Bioconductor 依赖，而它算出的"活性"与本脚本的加权 z-score
# 在**排序上高度一致**。这里选依赖更少的路，把置换检验留给自己做。
#
# **这不等于 decoupleR。** 产物里写明了差异。
#
# ## 数据来源与回退
#
#   主路径：dorothea::dorothea_hs（confidence A/B/C 的调控子）
#   回退：assets/tf_regulons_minimal.tsv（少量手工整理的 TF-靶基因对）
#         —— **回退时 coverage 标为 minimal，结论不能当全基因组结论**
#   都没有：写 status = "not_done" + reason，不产出任何"结果"
#
# ## 产物
#
#   results/<GSE>/tf_regulon_enrichment.csv   每个 TF 的 Fisher 检验
#   results/<GSE>/tf_activity_by_sample.csv   每个样本 × TF 的活性分数
#   results/<GSE>/tf_activity_group_test.csv  组间比较
#   results/<GSE>/tf_regulon_enrichment.png        调控子富集条形图
#   results/<GSE>/tf_activity_group_difference.png 调控子活性组间效应量
#   results/<GSE>/tf_status.json                   状态与方法学限定
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

# 用哪些置信度的调控子。
# A = 人工审编（ChEA/TRRUST/文献），B = 中等，C = 预测。
# **不用 D/E** —— 那两档基本是基序预测，假阳性率高，
# 而富集分析对假阳性特别敏感（每个 TF 的靶基因集合被稀释）。
TF_CONFIDENCE <- c("A", "B", "C")

# 一个 TF 至少要这么多靶基因落在检测到的基因集里，才做检验。
# 太少的话 Fisher 的自由度不足，p 值没有意义。
MIN_TARGETS_IN_DATA <- 5L

# 一个 TF 的活性分数至少要这么多靶基因才计算（噪声控制）
MIN_TARGETS_FOR_ACTIVITY <- 5L


#' 读入 TF-靶基因调控子表
#'
#' 返回 list(regulons = data.frame(tf, target, mor, confidence), source, coverage)
#' mor = mode of regulation（+1 激活 / -1 抑制）
load_tf_regulons <- function(cfg) {
  # ---- 主路径：dorothea ----
  if (requireNamespace("dorothea", quietly = TRUE)) {
    db <- dorothea::dorothea_hs
    # 列名在不同版本间可能是 tf/target/mor/confidence
    need <- c("tf", "target", "mor", "confidence")
    if (all(need %in% colnames(db))) {
      keep <- db$confidence %in% TF_CONFIDENCE
      db <- db[keep, need, drop = FALSE]
      log_info(sprintf("dorothea: %d 个 TF, %d 条 TF-靶基因关系（置信度 %s）",
                       length(unique(db$tf)), nrow(db),
                       paste(TF_CONFIDENCE, collapse = "/")))
      return(list(regulons = as.data.frame(db),
                  source = "dorothea::dorothea_hs",
                  coverage = "genome_wide",
                  confidence_levels = TF_CONFIDENCE))
    }
    log_warn(sprintf("dorothea 的列名不是预期的 %s，实际 %s —— 转回退路径",
                     paste(need, collapse = "/"),
                     paste(colnames(db), collapse = "/")))
  } else {
    log_warn("dorothea 未安装 —— 转回退路径")
  }

  # ---- 回退：随仓库附带的最小表 ----
  # 不假设 cwd 是仓库根 —— 从 --file= 反推脚本位置，再往上找 assets/。
  # 找不到就按 cwd 试一次（编排器从根目录运行时走这条）。
  cand <- character(0)
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    here <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
    cand <- c(cand, file.path(here, "..", "assets", "tf_regulons_minimal.tsv"),
              file.path(here, "assets", "tf_regulons_minimal.tsv"))
  }
  cand <- c(cand, file.path("assets", "tf_regulons_minimal.tsv"))
  cand <- unique(cand)
  fb <- cand[file.exists(cand)]
  if (length(fb) > 0L) {
    fb <- fb[[1L]]
    # **`comment.char = "#"` 不能省。** `read.delim` 默认 `comment.char = ""`，
    # 也就是**不把 `#` 当注释** —— 表头那 20 行说明文字会被当成数据行读进来，
    # 然后 `tf` 列里出现一堆 `# ===...`，靶基因匹配全部落空。
    # 而这个失败不报错，只是富集结果为空。
    db <- utils::read.delim(fb, stringsAsFactors = FALSE, comment.char = "#")
    if (all(c("tf", "target", "mor") %in% colnames(db))) {
      db$confidence <- "minimal"
      log_warn(sprintf(
        "**回退到最小 TF 表**：%d 个 TF, %d 条关系。覆盖度远低于 dorothea，",
        length(unique(db$tf)), nrow(db)))
      log_warn("结论只能作为**示例**，不能当全基因组的 TF 调控结论。")
      return(list(regulons = db,
                  source = "assets/tf_regulons_minimal.tsv (回退)",
                  coverage = "minimal",
                  confidence_levels = "minimal"))
    }
  }

  # ---- 都没有 ----
  NULL
}


#' 调控子富集：Fisher 精确检验
#'
#' 2x2 表：
#'              在 DEG 里   不在 DEG 里
#'   是靶基因        a            b
#'   不是靶基因      c            d
#'
#' 背景集用**检测到的基因**（detected），不用全基因组 —— 与
#' 04_heatmap_enrichment.R 的 ORA 保持一致（AGENTS.md 规则 9）。
regulon_enrichment <- function(regulons, deg_genes, detected_genes) {
  det <- unique(detected_genes)
  deg <- intersect(unique(deg_genes), det)
  n_det <- length(det)
  n_deg <- length(deg)
  if (n_deg == 0L) return(NULL)

  tfs <- sort(unique(regulons$tf))
  rows <- vector("list", length(tfs))
  k <- 0L
  for (tf in tfs) {
    tg <- intersect(unique(regulons$target[regulons$tf == tf]), det)
    if (length(tg) < MIN_TARGETS_IN_DATA) next
    a <- length(intersect(tg, deg))
    b <- length(tg) - a
    cc <- n_deg - a
    d <- n_det - n_deg - b
    if (d < 0L) next
    ft <- stats::fisher.test(matrix(c(a, b, cc, d), nrow = 2L),
                             alternative = "greater")
    k <- k + 1L
    rows[[k]] <- data.frame(
      tf = tf,
      n_targets_in_data = length(tg),
      n_targets_in_deg = a,
      expected_in_deg = round(length(tg) * n_deg / n_det, 3),
      odds_ratio = round(unname(ft$estimate), 4),
      p_value = ft$p.value,
      stringsAsFactors = FALSE)
  }
  if (k == 0L) return(NULL)
  out <- do.call(rbind, rows[seq_len(k)])
  out$p_adj_bh <- stats::p.adjust(out$p_value, method = "BH")
  out <- out[order(out$p_value), , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' 调控子活性分数（按样本）
#'
#' 分数 = 该 TF 靶基因 z-score 的加权均值，符号由 mor 决定。
#' 这是经典的 regulon score，**不是 decoupleR 的 ULM**。
regulon_activity <- function(expr, regulons, detected_genes) {
  # 按基因（行）z-score：让不同基因可比
  z <- t(scale(t(expr)))
  z[!is.finite(z)] <- 0
  det <- intersect(rownames(z), detected_genes)
  z <- z[det, , drop = FALSE]

  tfs <- sort(unique(regulons$tf))
  scores <- list()
  for (tf in tfs) {
    sub <- regulons[regulons$tf == tf, , drop = FALSE]
    sub <- sub[sub$target %in% rownames(z), , drop = FALSE]
    if (nrow(sub) < MIN_TARGETS_FOR_ACTIVITY) next
    m <- z[sub$target, , drop = FALSE]
    # mor 定符号；缺失 mor 时按 +1（激活）处理并记录
    mor <- ifelse(is.na(sub$mor), 1, sign(sub$mor))
    scores[[tf]] <- as.numeric(crossprod(mor, m) / sum(abs(mor)))
  }
  if (length(scores) == 0L) return(NULL)
  out <- as.data.frame(do.call(cbind, scores), stringsAsFactors = FALSE)
  colnames(out) <- names(scores)
  rownames(out) <- colnames(z)
  out
}


#' 组间比较：每个 TF 的活性分数在两组间是否有差异
#'
#' **n<10 时这里只报效应量，不把 p 值当结论。**
#' 用 Wilcoxon（不假设正态）与组间均值差（效应量）两个都报。
group_test <- function(act, group) {
  grp <- as.character(group)
  lv <- sort(unique(grp))
  if (length(lv) != 2L) return(NULL)
  rows <- vector("list", ncol(act))
  for (j in seq_len(ncol(act))) {
    v <- act[[j]]
    g1 <- v[grp == lv[1]]
    g2 <- v[grp == lv[2]]
    wt <- suppressWarnings(stats::wilcox.test(g1, g2, exact = FALSE))
    rows[[j]] <- data.frame(
      tf = colnames(act)[j],
      mean_group1 = round(mean(g1), 4),
      mean_group2 = round(mean(g2), 4),
      # 效应量：Cohen's d（小样本下也只是参考）
      cohens_d = round((mean(g2) - mean(g1)) /
                         sqrt((stats::var(g1) + stats::var(g2)) / 2), 4),
      delta_mean = round(mean(g2) - mean(g1), 4),
      p_value = wt$p.value,
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out$p_adj_bh <- stats::p.adjust(out$p_value, method = "BH")
  out$group1 <- lv[1]
  out$group2 <- lv[2]
  out[order(out$p_value), , drop = FALSE]
}


run_08_tf_regulation <- function(cfg) {
  set.seed(cfg$analysis$seed)
  res <- cfg$output$results_dir
  dat <- cfg$output$data_dir
  status_path <- file.path(res, "tf_status.json")

  # ---- 输入 ---------------------------------------------------------------
  deg_path <- file.path(res, "deg_table.csv")
  expr_path <- file.path(dat, "expr_clean.rds")
  if (!file.exists(deg_path) || !file.exists(expr_path)) {
    write_json(status_path, list(
      status = "not_configured",
      reason = paste("缺少输入：",
                     if (!file.exists(deg_path)) "deg_table.csv" else "",
                     if (!file.exists(expr_path)) "expr_clean.rds" else ""),
      note = "TF 调控分析需要 limma 的差异结果与清洗后的表达矩阵"))
    log_warn("TF 调控：缺输入，记 not_configured")
    return(invisible(NULL))
  }

  deg <- utils::read.csv(deg_path, stringsAsFactors = FALSE)
  expr <- readRDS(expr_path)
  group <- utils::read.csv(file.path(dat, "group.csv"),
                           stringsAsFactors = FALSE)

  # ---- 调控子 -------------------------------------------------------------
  tfdb <- load_tf_regulons(cfg)
  if (is.null(tfdb)) {
    write_json(status_path, list(
      status = "not_done",
      reason = paste(
        "既没有 dorothea 包，也没有 assets/tf_regulons_minimal.tsv。",
        "**没有 TF-靶基因关系就无法做调控分析** —— 不用代理指标冒充。",
        "装 dorothea：BiocManager::install('dorothea')"),
      method = NULL))
    log_warn("TF 调控：没有调控子来源，记 not_done（不用代理指标冒充）")
    return(invisible(NULL))
  }

  detected <- rownames(expr)

  # ---- 分析 1：调控子富集 -------------------------------------------------
  # 用**显著 DEG**；没有显著基因时用最显著的 N 个（与 deg_mode 的措辞约束一致）
  deg_sig <- deg$gene[!is.na(deg$adj.P.Val) & deg$adj.P.Val < 0.05]
  deg_mode <- if (length(deg_sig) > 0L) "significant" else "ranked_fallback"
  if (deg_mode == "ranked_fallback") {
    n_fb <- min(200L, nrow(deg))
    deg_sig <- deg$gene[seq_len(n_fb)]
    log_warn(sprintf(
      "没有 FDR<0.05 的基因 —— 回退到最显著的 %d 个（deg_mode=ranked_fallback）。",
      n_fb))
    log_warn("**措辞约束**：只能说『在最显著的 N 个基因里富集到……』，")
    log_warn("不能说『显著差异基因富集到……』（AGENTS.md 规则 3 / 设计文档 §2.10）")
  }

  enr <- regulon_enrichment(tfdb$regulons, deg_sig, detected)
  if (!is.null(enr)) {
    utils::write.csv(enr, file.path(res, "tf_regulon_enrichment.csv"),
                     row.names = FALSE)
    n_sig <- sum(enr$p_adj_bh < 0.05, na.rm = TRUE)
    log_info(sprintf("调控子富集: %d 个 TF 可检验，BH<0.05 的 %d 个",
                     nrow(enr), n_sig))
    if (nrow(enr) > 0L) {
      log_info(sprintf("  最显著: %s", paste(
        sprintf("%s(p=%.2g)", utils::head(enr$tf, 5),
                utils::head(enr$p_value, 5)), collapse = ", ")))
    }
  } else {
    log_warn("调控子富集：没有 TF 满足最小靶基因数要求")
  }

  # ---- 分析 2：调控子活性 -------------------------------------------------
  act <- regulon_activity(expr, tfdb$regulons, detected)
  act_test <- NULL
  if (!is.null(act)) {
    utils::write.csv(data.frame(sample = rownames(act), act,
                                check.names = FALSE),
                     file.path(res, "tf_activity_by_sample.csv"),
                     row.names = FALSE)
    log_info(sprintf("调控子活性: %d 个样本 x %d 个 TF", nrow(act), ncol(act)))

    # 分组向量与表达矩阵的列对齐（**顺序错了不报错，只会算出错的组间差**）
    g <- group$group[match(rownames(act), group$sample)]
    if (any(is.na(g))) {
      log_warn(sprintf("有 %d 个样本在 group.csv 里找不到分组，已剔除",
                       sum(is.na(g))))
      keep <- !is.na(g)
      act <- act[keep, , drop = FALSE]
      g <- g[keep]
    }
    act_test <- group_test(act, g)
    if (!is.null(act_test)) {
      utils::write.csv(act_test, file.path(res, "tf_activity_group_test.csv"),
                       row.names = FALSE)
      log_info(sprintf("组间比较: %d 个 TF，BH<0.05 的 %d 个",
                       nrow(act_test), sum(act_test$p_adj_bh < 0.05, na.rm = TRUE)))
    }
  } else {
    log_warn("调控子活性：没有 TF 满足最小靶基因数要求")
  }

  # ---- 出图 ---------------------------------------------------------------
  # **绘图单独兜住，不让画图错误影响方法本身的记录**（AGENTS.md 规则 14）。
  # 两张图分开写，不用 patchwork —— 少一个依赖，两张图本来也回答两个问题。
  plot_ok <- tryCatch({
    n_show <- min(20L, if (!is.null(enr)) nrow(enr) else 0L)
    if (n_show >= 1L) {
      top <- utils::head(enr, n_show)
      top$tf <- factor(top$tf, levels = rev(top$tf))
      p1 <- ggplot2::ggplot(top, ggplot2::aes(x = .data$tf,
                                              y = -log10(.data$p_value))) +
        ggplot2::geom_col(fill = PAL$up, width = 0.7) +
        ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                            colour = PAL$muted, linewidth = 0.4) +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = sprintf("TF regulon enrichment among DEGs - %s",
                          cfg$dataset_id),
          subtitle = wrap_subtitle(sprintf(
            "Fisher exact test, background = %d detected genes; %s. Dashed line = p 0.05 (uncorrected)",
            length(detected),
            if (identical(deg_mode, "significant"))
              sprintf("%d FDR<0.05 genes", length(deg_sig))
            else sprintf("ranked_fallback: top %d genes (no FDR<0.05)",
                         length(deg_sig)))),
          x = NULL, y = expression(-log[10](italic(p)))) +
        theme_paper()
      ggplot2::ggsave(file.path(res, "tf_regulon_enrichment.png"), p1,
                      width = 7.0, height = 6.0,
                      dpi = cfg$analysis$figure_dpi, bg = "white")

      if (!is.null(act_test) && nrow(act_test) >= 2L) {
        tt <- utils::head(act_test[order(-abs(act_test$cohens_d)), ], 20L)
        tt$tf <- factor(tt$tf, levels = rev(tt$tf))
        tt$direction <- ifelse(tt$delta_mean > 0,
                               sprintf("higher in %s", tt$group2[1]),
                               sprintf("higher in %s", tt$group1[1]))
        # 方向色板：红 = 在 group2 高，蓝 = 在 group1 高。
        # **同一个红在火山图里也是"上调"** —— 一个颜色一个含义。
        lv2 <- sprintf("higher in %s", tt$group2[1])
        lv1 <- sprintf("higher in %s", tt$group1[1])
        p2 <- ggplot2::ggplot(tt, ggplot2::aes(x = .data$tf,
                                               y = .data$cohens_d,
                                               fill = .data$direction)) +
          ggplot2::geom_col(width = 0.7) +
          ggplot2::coord_flip() +
          ggplot2::scale_fill_manual(values = stats::setNames(
            c(PAL$up, PAL$down), c(lv2, lv1))) +
          ggplot2::labs(
            title = sprintf("TF regulon activity, group difference - %s",
                            cfg$dataset_id),
            subtitle = wrap_subtitle(sprintf(
              "Cohen's d (%s vs %s), min n=%d per group. n is small - read effect size, not p",
              act_test$group2[1], act_test$group1[1], min(table(g)))),
            x = NULL, y = "Cohen's d", fill = NULL) +
          theme_paper()
        ggplot2::ggsave(file.path(res, "tf_activity_group_difference.png"), p2,
                        width = 7.0, height = 6.0,
                        dpi = cfg$analysis$figure_dpi, bg = "white")
      }
      TRUE
    } else {
      FALSE
    }
  }, error = function(e) {
    log_warn(sprintf("TF 图绘制失败（方法本身不受影响）: %s", conditionMessage(e)))
    FALSE
  })

  # ---- 状态 ---------------------------------------------------------------
  n_sig_enr <- if (!is.null(enr)) sum(enr$p_adj_bh < 0.05, na.rm = TRUE) else 0L
  n_sig_act <- if (!is.null(act_test)) {
    sum(act_test$p_adj_bh < 0.05, na.rm = TRUE)
  } else 0L
  n_per_group <- if (exists("g") && length(g) > 0L) {
    as.list(table(g))
  } else {
    NULL
  }

  write_json(status_path, list(
    status = "ok",
    regulon_source = tfdb$source,
    coverage = tfdb$coverage,
    confidence_levels = tfdb$confidence_levels,
    n_tfs_total = length(unique(tfdb$regulons$tf)),
    n_interactions = nrow(tfdb$regulons),
    deg_mode = deg_mode,
    n_genes_in_deg_set = length(deg_sig),
    n_genes_detected = length(detected),
    n_tfs_testable = if (!is.null(enr)) nrow(enr) else 0L,
    n_tfs_enriched_bh05 = n_sig_enr,
    n_tfs_activity_tested = if (!is.null(act_test)) nrow(act_test) else 0L,
    n_tfs_activity_bh05 = n_sig_act,
    n_per_group = n_per_group,
    figure_written = plot_ok,
    method = paste(
      "调控子富集 = Fisher 精确检验（背景集 = 检测到的基因）；",
      "调控子活性 = 靶基因 z-score 的 mor 加权均值"),
    not_decoupler = paste(
      "**这不是 decoupleR 的 ULM。** decoupleR 用单样本线性模型 + 置换检验；",
      "本脚本用加权 z-score。两者在**排序上**通常一致，但 p 值与绝对值不同。",
      "要 decoupleR 的结果请另装该包。"),
    limitations = c(
      sprintf("**每组 n 很小**（%s）。调控子活性分数的组间比较只应看效应量，",
              if (!is.null(n_per_group))
                paste(names(n_per_group), unlist(n_per_group),
                      sep = "=", collapse = ", ") else "未知"),
      "      p 值不可作为结论 —— 小样本下 Wilcoxon 的功效极低。",
      sprintf("调控子覆盖度 = %s。", tfdb$coverage),
      if (identical(tfdb$coverage, "minimal"))
        "     **回退表只覆盖少量 TF，结论只能作为示例。**"
      else "     但仍不含低置信度（D/E）的预测调控子。",
      "TF 调控是**关联**，不是因果。表达相关不等于该 TF 在调控这些靶基因。",
      "bulk 数据里 TF 的 mRNA 水平与其**蛋白活性**经常不相关",
      "（翻译后修饰、核转位都不反映在 mRNA 上）—— 这是本分析的根本限制。",
      "靶基因集合之间大量重叠（一个基因受多个 TF 调控），",
      "所以各 TF 的 p 值**不独立**，BH 校正偏保守。"),
    outputs = c("tf_regulon_enrichment.csv", "tf_activity_by_sample.csv",
                "tf_activity_group_test.csv", "tf_regulon_enrichment.png",
                "tf_activity_group_difference.png")))

  log_info("TF 调控分析完成")
  invisible(NULL)
}


if (!GEO_ORCHESTRATED()) {
  cfg <- parse_args()
  run_08_tf_regulation(cfg)
}
