# ============================================================================
# 03_deg.R — limma 差异表达分析
# ============================================================================
# spec 的 limma_deg 步骤。
#
# 设计矩阵 ~ 0 + group，对比 = contrast[1] - contrast[2]。
# 输出：results/deg_table.csv（含 gene/logFC/AveExpr/t/P.Value/adj.P.Val/B）
#       results/volcano_plot.pdf
#       results/deg_summary.json
# ============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
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

run_03_deg <- function(cfg) {
  log_info("=== 步骤 03：limma 差异表达 ===")
  ensure_dirs(cfg)

  res <- cfg$output$results_dir
  expr <- readRDS(file.path(cfg$output$data_dir, "expr_clean.rds"))
  group <- utils::read.csv(file.path(cfg$output$data_dir, "group.csv"), stringsAsFactors = FALSE)

  numerator   <- cfg$contrast[1L]
  denominator <- cfg$contrast[2L]

  groups <- factor(group$group, levels = unique(group$group))
  if (!all(c(numerator, denominator) %in% levels(groups))) {
    stop(sprintf("分组中找不到 contrast 指定的组: %s / %s", numerator, denominator))
  }
  # 让分母成为参考水平，便于解读
  groups <- stats::relevel(groups, ref = denominator)

  # ---- 0. 设计矩阵：配对 vs 非配对 ----------------------------------------
  #
  # 配对设计把患者作为阻断因子。GSE64790 这类"肿瘤 + 同一患者正常组织"的研究里，
  # 患者间差异往往比肿瘤/正常差异还大，不阻断就会把信号埋进残差。
  #
  # 但配对会消耗自由度（3 对只剩 2 df），且对不完整配对很敏感，所以这里
  # **失败就退回非配对**并记录原因，绝不让它把整条流水线拖垮。
  paired_used <- FALSE
  pair_reason <- NA_character_
  design <- NULL

  if (isTRUE(cfg$paired)) {
    # 只把「能不能配对」交给 tryCatch；设计矩阵与标志位走正常赋值。
    # 早先版本在 tryCatch 表达式里写 paired_used <<- TRUE —— 那个表达式在
    # run_03_deg 的帧里求值，<<- 会跳过当前帧写到外层环境，导致本地
    # paired_used 一直是 FALSE：设计矩阵用了配对，摘要却报 paired=false。
    pat <- tryCatch({
      if (is.null(group$patient)) stop("group.csv 里没有 patient 列")
      if (anyNA(group$patient) || !all(nzchar(group$patient))) stop("存在未配对的样本")
      p <- factor(group$patient)
      if (nlevels(p) < 2L) stop(sprintf("只有 %d 个配对水平，无法阻断", nlevels(p)))
      p
    }, error = function(e) {
      # 这里是嵌套函数，<<- 才会正确落到 run_03_deg 的帧
      pair_reason <<- conditionMessage(e)
      log_warn(sprintf("配对设计不可用，退回非配对: %s", pair_reason))
      NULL
    })

    if (!is.null(pat)) {
      d <- stats::model.matrix(~ 0 + groups + pat)
      if (qr(d)$rank < ncol(d)) {
        pair_reason <- "设计矩阵秩不足（配对不完整？）"
        log_warn(sprintf("配对设计不可用，退回非配对: %s", pair_reason))
      } else {
        design <- d
        paired_used <- TRUE
        log_info(sprintf("配对设计: %d 对，阻断因子 patient (%d 水平)，残差 df = %d",
                         nlevels(pat), nlevels(pat), nrow(d) - qr(d)$rank))
      }
    }
  }

  if (is.null(design)) {
    design <- stats::model.matrix(~ 0 + groups)
  }
  colnames(design) <- sub("^groups", "", colnames(design))
  log_info(sprintf("设计矩阵: %s", paste(colnames(design), collapse = ", ")))
  log_info(sprintf("分组样本数: %s", paste(sprintf("%s=%d", names(table(groups)),
                                                   as.integer(table(groups))), collapse = ", ")))

  # ---- 1. 线性模型 + 经验贝叶斯 -------------------------------------------
  fit <- limma::lmFit(expr, design)
  cont <- limma::makeContrasts(contrasts = sprintf("%s - %s", numerator, denominator),
                               levels = design)
  fit2 <- limma::contrasts.fit(fit, cont)
  # robust 需要 statmod；缺失时退回标准 eBayes，不中断流程
  fit2 <- tryCatch(
    limma::eBayes(fit2, trend = TRUE, robust = TRUE),
    error = function(e) {
      log_warn(sprintf("robust eBayes 不可用，退回标准 eBayes: %s", conditionMessage(e)))
      limma::eBayes(fit2, trend = TRUE)
    }
  )

  # ---- 2. 全基因结果表 ----------------------------------------------------
  tt <- limma::topTable(fit2, number = Inf, sort.by = "P", adjust.method = "BH")
  tt$gene <- rownames(tt)
  tt <- tt[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

  required_cols <- c("gene", "logFC", "P.Value", "adj.P.Val")
  stopifnot(all(required_cols %in% colnames(tt)))
  if (anyNA(tt$adj.P.Val)) {
    log_warn(sprintf("%d 个基因的 adj.P.Val 为 NA（方差为 0 等），已置为 1", sum(is.na(tt$adj.P.Val))))
    tt$adj.P.Val[is.na(tt$adj.P.Val)] <- 1
    tt$P.Value[is.na(tt$P.Value)] <- 1
  }
  utils::write.csv(tt, file.path(res, "deg_table.csv"), row.names = FALSE)

  # ---- 3. 显著基因 --------------------------------------------------------
  padj_cut <- cfg$thresholds$adj_p
  lfc_cut  <- cfg$thresholds$log2fc
  sig <- tt[tt$adj.P.Val < padj_cut & abs(tt$logFC) > lfc_cut, , drop = FALSE]
  sig_up   <- sum(sig$logFC > 0)
  sig_down <- sum(sig$logFC < 0)

  log_info(sprintf("显著 DEG (adj.P < %g 且 |log2FC| > %g): %d 个（上调 %d / 下调 %d）",
                   padj_cut, lfc_cut, nrow(sig), sig_up, sig_down))

  if (nrow(sig) == 0L) {
    log_warn("没有基因通过显著性阈值。下游热图会降级为 top 20 by P，富集会输出空表。")
  }
  if (nrow(sig) < 10L) {
    log_warn(sprintf("显著基因仅 %d 个，GO/KEGG 富集很可能无结果（小样本设计的正常表现）",
                     nrow(sig)))
  }

  # ---- 4. 火山图 ----------------------------------------------------------
  tt$status <- "ns"
  tt$status[tt$adj.P.Val < padj_cut & tt$logFC >  lfc_cut] <- "up"
  tt$status[tt$adj.P.Val < padj_cut & tt$logFC < -lfc_cut] <- "down"
  tt$status <- factor(tt$status, levels = c("up", "down", "ns"))
  tt$neg_log10_p <- -log10(pmax(tt$adj.P.Val, .Machine$double.xmin))

  label_df <- rbind(
    head(tt[tt$status == "up", , drop = FALSE][order(-tt$logFC[tt$status == "up"]), ], 10),
    head(tt[tt$status == "down", , drop = FALSE][order(tt$logFC[tt$status == "down"]), ], 10)
  )

  p_volcano <- ggplot(tt, aes(x = logFC, y = neg_log10_p, colour = status)) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_colour_manual(values = c(up = "#C1443C", down = "#2E5FA3", ns = "grey75"),
                        labels = c(up = sprintf("up (%d)", sig_up),
                                   down = sprintf("down (%d)", sig_down),
                                   ns = "not significant")) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", linewidth = 0.3) +
    labs(title = sprintf("Volcano: %s vs %s (%s)", numerator, denominator, cfg$dataset_id),
         subtitle = sprintf("adj.P < %g, |log2FC| > %g", padj_cut, lfc_cut),
         x = sprintf("log2 fold change (%s / %s)", numerator, denominator),
         y = "-log10 adjusted P value", colour = NULL) +
    theme_bw(base_size = 10)

  if (nrow(label_df) > 0L) {
    p_volcano <- p_volcano +
      ggrepel_labels(label_df)
  }
  save_pdf(file.path(res, "volcano_plot.pdf"), print(p_volcano), width = 8, height = 6.5)
  log_info("已生成 volcano_plot.pdf")

  # ---- 5. p 值分布诊断（DE 之后的 QC 关卡）--------------------------------
  #
  # 来自 K-Dense `bulk-rnaseq` skill 的设计/QC 清单：
  #   "A well-behaved test gives a roughly uniform histogram with a peak near 0
  #    (the true positives). A peak near 1, or a U-shape, signals a problem:
  #    misspecified design, unmodeled batch, or filtering issues.
  #    Fix the design rather than trusting the gene list."
  #
  # 这张图对本设计尤其重要：n=6 时几乎不可能有基因通过 FDR，
  # 光看"0 个显著基因"分不清是**功效不足**还是**设计有问题**。
  # p 值直方图能把这两者区分开 —— 0 附近有峰说明信号真实、只是检不出；
  # 峰在 1 或 U 形说明模型设定错了，那时候连排序表都不能用。
  pv <- tt$P.Value[!is.na(tt$P.Value)]
  p_hist <- ggplot2::ggplot(data.frame(p = pv), ggplot2::aes(x = p)) +
    ggplot2::geom_histogram(bins = 40, boundary = 0, fill = "grey55",
                            colour = "white", linewidth = 0.2) +
    ggplot2::geom_hline(yintercept = length(pv) / 40, linetype = "dashed",
                        colour = "firebrick", linewidth = 0.4) +
    ggplot2::labs(
      title = sprintf("P value distribution: %s vs %s (%s)", numerator, denominator,
                      cfg$dataset_id),
      subtitle = sprintf("%d genes tested; dashed line = uniform expectation. %s",
                         length(pv),
                         if (length(pv) && min(pv) < 1e-3)
                           sprintf("%.0f genes at P < 0.001 (real signal below the null)",
                                   sum(pv < 1e-3))
                         else "no gene below P < 0.001"),
      x = "raw P value", y = "gene count") +
    ggplot2::theme_bw(base_size = 10)
  save_pdf(file.path(res, "pvalue_histogram.pdf"), print(p_hist), width = 7, height = 5)

  # 诊断结论：把"功效不足"和"设计有问题"分开
  n_below_001 <- sum(pv < 0.001)
  # 在正确的原假设下，P < 0.05 的基因数应约为 5%。明显超出 = 有真实信号。
  excess <- sum(pv < 0.05) / max(1, length(pv)) / 0.05
  pval_verdict <- if (excess >= 2) {
    "signal_present_but_underpowered"   # 有信号，只是过不了多重检验
  } else if (excess <= 0.5) {
    "little_or_no_signal"
  } else {
    "inconclusive"
  }
  log_info(sprintf("P 值分布: P<0.05 的基因占 %.1f%%（原假设期望 5%%，超出 %.1f 倍）→ %s",
                   100 * sum(pv < 0.05) / max(1, length(pv)), excess, pval_verdict))
  if (pval_verdict == "signal_present_but_underpowered") {
    log_warn(sprintf(paste0("有真实信号但功效不足：%d 个基因 P < 0.001，",
                            "却没有任何基因通过 FDR。这是样本量的限制，不是设计错误。"),
                     n_below_001))
  }

  # ---- 6. 摘要 ------------------------------------------------------------
  utils::write.csv(sig, file.path(res, "deg_significant.csv"), row.names = FALSE)
  write_json(file.path(res, "deg_summary.json"), list(
    dataset_id = cfg$dataset_id,
    contrast = sprintf("%s - %s", numerator, denominator),
    numerator = numerator, denominator = denominator,
    n_samples = ncol(expr), n_genes_tested = nrow(tt),
    adj_p_cutoff = padj_cut, log2fc_cutoff = lfc_cut,
    n_significant = nrow(sig), n_up = sig_up, n_down = sig_down,
    paired = paired_used,
    paired_requested = isTRUE(cfg$paired),
    paired_fallback_reason = pair_reason,
    residual_df = ncol(expr) - qr(design)$rank,
    n_p_below_0.001 = n_below_001,
    frac_p_below_0.05 = round(sum(pv < 0.05) / max(1, length(pv)), 4),
    pvalue_diagnosis = pval_verdict,
    top10_by_p = head(tt$gene, 10)
  ))

  log_info("已生成 deg_table.csv / deg_significant.csv / deg_summary.json / pvalue_histogram.pdf")
  invisible(tt)
}

#' 火山图基因标注（ggrepel 存在时用它，否则退回 geom_text）
ggrepel_labels <- function(label_df) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(
      data = label_df, aes(label = gene), size = 2.3, max.overlaps = 20,
      segment.size = 0.2, show.legend = FALSE
    )
  } else {
    ggplot2::geom_text(
      data = label_df, aes(label = gene), size = 2.3,
      vjust = -0.6, show.legend = FALSE, check_overlap = TRUE
    )
  }
}

if (!GEO_ORCHESTRATED()) {
  cfg <- load_config()
  run_03_deg(cfg)
}
