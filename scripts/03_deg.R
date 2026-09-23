# ============================================================================
# 03_deg.R — limma 差异表达分析
# ============================================================================
# spec 的 limma_deg 步骤。
#
# 设计矩阵 ~ 0 + group，对比 = contrast[1] - contrast[2]。
# 输出：results/deg_table.csv（含 gene/logFC/AveExpr/t/P.Value/adj.P.Val/B）
#       results/01-03-01-unit1-volcano-plot.pdf
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
  #
  # **颜色只承载方向，显著性用 alpha + 大小承载。**
  #
  # 判据 `adj.P < adj_p` 且 `|log2FC| > log2fc` 是 AND 关系，n=6 时全基因组
  # BH 校正下**没有任何基因**能过 adj.P 这一关（GSE64790 最小 adj.P = 0.394）。
  # 早先的做法是把这批"名义显著"（raw P < 0.05 且倍数过关，但未过 FDR）
  # 单独染成橙色 —— 结果图上**看不到任何上调/下调**，因为红蓝两档是空的。
  #
  # 现在改成两个正交通道：
  #   * `direction`（颜色）：up = 红、down = 蓝、不显著 = 灰
  #   * `tier`（alpha + 大小）：FDR 显著 = 实心大点，名义显著 = 半透明小点
  # 这样上调红、下调蓝一眼可见，而"未过 FDR"这个事实仍然写在图例和副标题里，
  # 不会因为染成方向色就被当成显著。**两个通道编码两个变量，不是一色两义。**
  #
  # 纵轴统一用 **raw P**（这也是火山图更常见的约定）。
  #
  # **不能把 raw P 和 adj.P 混在一条纵轴上。** 曾经想让名义显著的点用 raw P、
  # 其余用 adj.P，但那样名义显著的点会画在 `-log10(adj_p) = 1.3` 这条线**之上**，
  # 看起来像通过了 FDR —— 正好是这张图要避免的误导。
  # 而全部用 adj.P 也不行：GSE64790 的 adj.P 全在 0.394~1 之间，
  # 纵轴范围只有 0~0.4，整个图压成一张饼。
  #
  # 所以纵轴一律 raw P，横线标 `P = 0.05`（含义与轴一致）。
  # FDR 阈值**无法**用一条横线表示 —— BH 校正是逐基因的，没有常数截断。
  # 这一点写在副标题里。
  tt$direction <- ifelse(tt$logFC > 0, "up", "down")
  tt$direction[tt$logFC == 0] <- "ns"
  tt$tier <- "ns"
  tt$tier[!is.na(tt$P.Value) & tt$P.Value < 0.05 & abs(tt$logFC) > lfc_cut] <- "nominal"
  fdr_hit <- !is.na(tt$adj.P.Val) & tt$adj.P.Val < padj_cut & abs(tt$logFC) > lfc_cut
  tt$tier[fdr_hit] <- "fdr"
  n_nominal <- sum(tt$tier == "nominal")
  n_fdr <- sum(tt$tier == "fdr")
  tt$direction[tt$tier == "ns"] <- "ns"
  log_info(sprintf("火山图分层: FDR 显著 %d（up %d / down %d）；名义显著（raw P<0.05，未过 FDR）%d（up %d / down %d）；其余 %d",
                   n_fdr, sig_up, sig_down, n_nominal,
                   sum(tt$tier == "nominal" & tt$direction == "up"),
                   sum(tt$tier == "nominal" & tt$direction == "down"),
                   sum(tt$tier == "ns")))

  tt$plot_y <- -log10(pmax(tt$P.Value, .Machine$double.xmin))

  label_df <- rbind(
    head(tt[tt$tier == "fdr" & tt$direction == "up", , drop = FALSE], 10),
    head(tt[tt$tier == "fdr" & tt$direction == "down", , drop = FALSE], 10)
  )
  # FDR 显著为空时，标注名义显著里最极端的那些，否则图上没有一个基因名
  if (nrow(label_df) == 0L && n_nominal > 0L) {
    cand <- tt[tt$tier == "nominal", , drop = FALSE]
    cand <- cand[order(cand$P.Value), , drop = FALSE]
    label_df <- rbind(head(cand[cand$direction == "up", , drop = FALSE], 8),
                      head(cand[cand$direction == "down", , drop = FALSE], 8))
  }

  # 图层顺序：不显著的先画（当背景），名义显著其次，FDR 显著最后（压在最上面）
  tt$tier <- factor(tt$tier, levels = c("ns", "nominal", "fdr"))
  tt$direction <- factor(tt$direction, levels = c("ns", "down", "up"))

  # **FDR 那句话必须随数据变。**
  # 原来 "No gene passes FDR" 是格式串里的**写死文字**，不随数据变：GSE64790
  # 恰好为真（0 个显著基因），但 GSE42568 有 3795 个基因通过 FDR —— 图上却
  # 仍然声称"没有基因通过 FDR"，而且与同一句里印出的 "min adj.P = 0.000"
  # 自相矛盾。**一张图上的事实陈述不能是常量。**
  # 顺带把 min adj.P 的格式从 %.3f 改成 %.3g：前者把 1e-12 印成 "0.000"，
  # 看上去像"完全没有信号"，而 %.3g 给 "1e-12"。
  min_adjp <- if (nrow(tt)) min(tt$adj.P.Val, na.rm = TRUE) else NA_real_
  fdr_note <- if (!is.finite(min_adjp)) {
    sprintf("FDR status undetermined (min adj.P = NA) for adj.P < %g. ", padj_cut)
  } else if (n_fdr > 0L) {
    sprintf(paste0("%d gene%s pass FDR (adj.P < %g; min adj.P = %.3g) and are ",
                   "drawn opaque; the semi-transparent points are nominal-P only. "),
            n_fdr, if (n_fdr == 1L) "" else "s", padj_cut, min_adjp)
  } else {
    sprintf(paste0("No gene passes FDR (adj.P < %g; min adj.P = %.3g), so the ",
                   "semi-transparent points are exploratory: nominal P only. "),
            padj_cut, min_adjp)
  }

  p_volcano <- ggplot(tt, aes(x = logFC, y = plot_y,
                              colour = direction, alpha = tier, size = tier)) +
    geom_point() +
    scale_colour_manual(
      values = c(up = PAL$up, down = PAL$down, ns = PAL$ns),
      labels = c(up = sprintf("%s up (%d)", numerator, sum(tt$direction == "up")),
                 down = sprintf("%s up (%d)", denominator, sum(tt$direction == "down")),
                 ns = "not significant"),
      name = "direction") +
    scale_alpha_manual(
      values = c(fdr = 0.95, nominal = 0.40, ns = 0.25),
      # **图例标签要短**（用户反馈"图例挤压主图"）：原标签
      # "nominal P < 0.05, NOT FDR-significant (1873)" 在右侧占了近 1/3 画布宽，
      # 主图被挤扁。计数保留（信息量在），措辞压到最短。
      labels = c(fdr = sprintf("FDR < %g (n=%d)", padj_cut, n_fdr),
                 nominal = sprintf("nominal only (n=%d)", n_nominal),
                 ns = sprintf("P >= 0.05 (n=%d)", sum(tt$tier == "ns"))),
      name = "significance") +
    scale_size_manual(values = c(fdr = 1.9, nominal = 1.1, ns = 0.6), guide = "none") +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed",
               linewidth = 0.3, colour = PAL$ink) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               linewidth = 0.3, colour = PAL$ink) +
    labs(title = sprintf("Volcano: %s vs %s (%s)", numerator, denominator, cfg$dataset_id),
         subtitle = wrap_subtitle(paste0(
           sprintf("|log2FC| > %g (vertical); horizontal line = nominal P 0.05. ", lfc_cut),
           sprintf("Red = higher in %s, blue = higher in %s. ", numerator, denominator),
           fdr_note,
           "BH is per-gene, so the FDR cutoff is not a horizontal line."),
           fig_width = W_DOUBLE),
         x = sprintf("log2 fold change (%s / %s)", numerator, denominator),
         y = "-log10 raw P value",
         colour = NULL, alpha = NULL) +
    theme_paper(10) +
    # **图例一律框外右侧、纵向**（用户约定 v2）：legend.direction 管单个图例内部
    # （键竖排）、legend.box 管多个图例之间（图例块竖摞）—— 两个都要设。
    # 字号再压一档，避免标签宽度把主图挤扁。
    ggplot2::theme(legend.direction = "vertical", legend.box = "vertical",
                   legend.text = ggplot2::element_text(size = 7),
                   legend.key.size = ggplot2::unit(0.55, "lines")) +
    guides(colour = guide_legend(order = 1, override.aes = list(alpha = 1, size = 2.2)),
           alpha = guide_legend(order = 2, override.aes = list(colour = PAL$ink, size = 2.2)))

  if (nrow(label_df) > 0L) {
    p_volcano <- p_volcano + ggrepel_labels(label_df, seed = cfg$analysis$seed)
  }
  save_pdf(file.path(res, "01-03-01-unit1-volcano-plot.pdf"), print(p_volcano), width = W_DOUBLE, height = mm(165))
  log_info("已生成 01-03-01-unit1-volcano-plot.pdf")

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
  # 直方图用中性灰、参考线用墨色：这张图里颜色**不承载方向语义**，
  # 用方向色会凭空暗示"蓝色 = 下调"。灰色恰好与 ns 的含义一致（零假设下的主体）。
  p_hist <- ggplot2::ggplot(data.frame(p = pv), ggplot2::aes(x = p)) +
    ggplot2::geom_histogram(bins = 40, boundary = 0, fill = PAL$ns,
                            colour = "white", linewidth = 0.2) +
    ggplot2::geom_hline(yintercept = length(pv) / 40, linetype = "dashed",
                        colour = PAL$ink, linewidth = 0.6) +
    ggplot2::labs(
      title = sprintf("P value distribution: %s vs %s (%s)", numerator, denominator,
                      cfg$dataset_id),
      subtitle = wrap_subtitle(sprintf(
        "%d genes tested; dashed line = uniform expectation. %s",
        length(pv),
        if (length(pv) && min(pv) < 1e-3)
          sprintf("%.0f genes at P < 0.001 (real signal below the null)", sum(pv < 1e-3))
        else "no gene below P < 0.001"),
        fig_width = mm(165)),
      x = "raw P value", y = "gene count") +
    theme_paper(10)
  save_pdf(file.path(res, "01-03-02-unit1-pvalue-histogram.pdf"), print(p_hist), width = mm(165), height = mm(127))

  # 诊断结论：把"功效不足"和"设计有问题"分开
  n_below_001 <- sum(pv < 0.001)
  # 在正确的原假设下，P < 0.05 的基因数应约为 5%。明显超出 = 有真实信号。
  excess <- sum(pv < 0.05) / max(1, length(pv)) / 0.05
  # **"有信号"和"功效不足"是两件事，必须分开判。**
  #
  # 原来写的是 `excess >= 2 -> "signal_present_but_underpowered"` —— 只要检出信号
  # 就贴上"功效不足"的标签，**根本没看有没有基因通过 FDR**。
  # 实测 GSE42568（n=121，3795 个基因过 FDR）也被标成 underpowered：
  # 这个字符串会直接进 deg_summary.json，任何引用它的措辞都会跟着错。
  #
  # 正确的判据是三分：有信号且**一个都没过** FDR 才叫功效不足；
  # 有信号且过了 FDR 就是功效充足。
  n_sig <- nrow(sig)
  pval_verdict <- if (excess >= 2 && n_sig == 0L) {
    "signal_present_but_underpowered"
  } else if (excess >= 2) {
    "signal_present"
  } else if (excess <= 0.5) {
    "little_or_no_signal"
  } else {
    "inconclusive"
  }
  log_info(sprintf("P 值分布: P<0.05 的基因占 %.1f%%（原假设期望 5%%，超出 %.1f 倍），FDR 显著 %d 个 → %s",
                   100 * sum(pv < 0.05) / max(1, length(pv)), excess, n_sig, pval_verdict))
  if (identical(pval_verdict, "signal_present_but_underpowered")) {
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
    n_nominal_only = n_nominal,
    paired = paired_used,
    paired_requested = isTRUE(cfg$paired),
    paired_fallback_reason = pair_reason,
    residual_df = ncol(expr) - qr(design)$rank,
    n_p_lt_0p001 = n_below_001,
    frac_p_lt_0p05 = round(sum(pv < 0.05) / max(1, length(pv)), 4),
    pvalue_diagnosis = pval_verdict,
    top10_by_p = head(tt$gene, 10)
  ))

  log_info("已生成 deg_table.csv / deg_significant.csv / deg_summary.json / 01-03-02-unit1-pvalue-histogram.pdf")
  invisible(tt)
}

#' 火山图基因标注（ggrepel 存在时用它，否则退回 geom_text）
#'
#' **必须传 seed。** `geom_text_repel` 用环境 RNG 做标签排布的模拟退火，
#' 不指定 seed 时位置取决于"此刻"的随机数状态 —— 而 `save_pdf` 会把同一个
#' 绘图表达式求值两次（一次 PDF、一次 PNG），第二次接着第一次消耗过的状态跑，
#' 于是两次的标签位置本来就不同；跨运行更会随上游随机数消耗量漂移。
#' 实测两轮 CI 的 12 张图里 11 张逐字节一致，只有火山图不一致，根因就是这里。
ggrepel_labels <- function(label_df, seed = NULL) {
  # ggrepel 的 seed 默认值是 NA（不是 NULL），NULL 会被当成缺参
  if (is.null(seed)) seed <- NA
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(
      data = label_df, aes(label = gene), size = 2.3, max.overlaps = 20,
      segment.size = 0.2, show.legend = FALSE, seed = seed,
      # 火山图把 alpha 映射到了 tier，标注不能继承 —— 否则名义显著的基因名
      # 会跟着变成半透明，正好是最需要看清的那些
      alpha = 1, fontface = "bold"
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
