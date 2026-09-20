# ============================================================================
# 02_qc_pca_correlation.R — 质控、PCA、样本间相关性
# ============================================================================
# spec 的 qc / pca / sample_correlation 三个步骤。
#
# 输出：results/boxplot_before_after.pdf
#       results/density_plot.pdf
#       results/pca_plot.pdf
#       results/correlation_heatmap.pdf
#       results/correlation_matrix.csv
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
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

#' 把矩阵转成 ggplot 需要的长表
matrix_to_long <- function(m, value_name = "expression") {
  df <- data.frame(
    sample = rep(colnames(m), each = nrow(m)),
    value  = as.vector(m),
    stringsAsFactors = FALSE
  )
  names(df)[names(df) == "value"] <- value_name
  df
}

#' 分组配色已挪到 common.R 的 `group_palette()` ——
#' 04 的热图注释条也要用它，定义在本文件会让 04 单独跑时找不到。
#' 这里曾经是本仓库最明显的一处配色不一致：旧实现把**第一个出现的分组**
#' 涂成蓝色，于是 PCA 上 tumor 是蓝的，而火山图、热图上 tumor（上调）是红的。
#' 同一含义两种颜色，读者要重新学一遍。现在颜色由 `cfg$contrast` 决定。

run_02_qc_pca_correlation <- function(cfg) {
  log_info("=== 步骤 02：QC / PCA / 样本相关性 ===")
  ensure_dirs(cfg)

  res <- cfg$output$results_dir
  expr_raw <- readRDS(file.path(cfg$output$data_dir, "expr_raw.rds"))
  expr <- readRDS(file.path(cfg$output$data_dir, "expr_clean.rds"))
  group <- utils::read.csv(file.path(cfg$output$data_dir, "group.csv"), stringsAsFactors = FALSE)
  groups <- factor(group$group, levels = unique(group$group))
  gp <- group_palette(group$group, cfg)

  # ---- 1. 标准化前后箱线图 ------------------------------------------------
  long <- rbind(
    cbind(matrix_to_long(expr_raw), stage = "before"),
    cbind(matrix_to_long(expr), stage = "after")
  )
  long$stage <- factor(long$stage, levels = c("before", "after"))

  # **x 轴样本名放不放得下，算出来，不靠默认。**
  # 标签旋转 90°，所以每个标签在**横向**占用的正是它的行高
  # `fontsize + min_gap` —— 与行名同一个一维模型，直接复用 decide_rownames()。
  # 关键区别：分面是 `ncol = 2`，**每个面板只有半幅宽**，可用长度必须按
  # 单个面板算，按整幅算会把容量高估一倍。
  # 实测 121 个样本、6pt、半幅 91mm、面板占宽 0.85：可容纳约 25 个 -> 放不下。
  # 原来无条件画 121 个 GSM 号，每个只剩约 0.75mm，糊成一条黑带。
  box_fs <- 6
  show_x <- decide_rownames(ncol(expr), W_DOUBLE / 2, box_fs, "箱线图样本名",
                            panel_frac = 0.85, min_gap = 2.5,
                            figure = "boxplot_before_after")

  p_box <- ggplot(long, aes(x = sample, y = expression, fill = stage)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.25) +
    facet_wrap(~stage, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = c(before = PAL$ns, after = PAL$down)) +
    labs(title = "Expression distribution before / after quantile normalization",
         subtitle = wrap_subtitle(sprintf(
           "%s - %d genes x %d samples. %s", cfg$dataset_id, nrow(expr), ncol(expr),
           if (isTRUE(show_x)) "x 轴标出样本编号。"
           else paste0("样本编号未标出（", ncol(expr),
                       " 个放不下）；逐个样本的身份见 data/", cfg$dataset_id,
                       "/group.csv。")),
           fig_width = W_DOUBLE),
         x = NULL, y = "log2 expression") +
    theme_paper(9) +
    theme(axis.text.x = if (isTRUE(show_x))
            element_text(angle = 90, hjust = 1, vjust = 0.5, size = box_fs)
          else element_blank(),
          # **刻度线跟着标签一起藏。** 121 个刻度线每个只占 0.75mm，密到连成
          # 一条黑带，看上去像图坏了。刻度线存在的意义是给标签定位 ——
          # 标签不画了，它就只剩副作用。
          axis.ticks.x = if (isTRUE(show_x)) element_line() else element_blank(),
          legend.position = "none")
  save_pdf(file.path(res, "boxplot_before_after.pdf"), print(p_box), width = W_DOUBLE, height = mm(140))
  log_info("已生成 boxplot_before_after.pdf")

  # ---- 2. 密度曲线 --------------------------------------------------------
  # 按**分组**上色（不是按样本），这样密度图和 PCA 用的是同一套条件色，
  # 一眼能看出某一组的分布是否整体偏移
  long$group <- groups[match(long$sample, group$gsm)]
  # 线型同样编码分组（颜色不是唯一载体）
  p_density <- ggplot(long, aes(x = expression, colour = group, linetype = group,
                                group = sample)) +
    geom_density(linewidth = 0.4, alpha = 0.85) +
    facet_wrap(~stage, ncol = 1, scales = "free_y") +
    scale_colour_condition(levels(groups), name = NULL) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")[seq_along(levels(groups))],
                          name = NULL) +
    labs(title = "Expression density before / after normalization",
         subtitle = wrap_subtitle(sprintf("coloured and styled by group; one curve per sample (%s)",
                                          cfg$dataset_id), fig_width = W_DOUBLE),
         x = "log2 expression", y = "density") +
    theme_paper(9)
  save_pdf(file.path(res, "density_plot.pdf"), print(p_density), width = W_DOUBLE, height = mm(152))
  log_info("已生成 density_plot.pdf")

  # ---- 3. PCA -------------------------------------------------------------
  top_n <- cfg$analysis$pca_top_genes
  if (is.null(top_n) || top_n <= 0 || top_n > nrow(expr)) {
    pca_input <- expr
  } else {
    v <- apply(expr, 1L, stats::var)
    pca_input <- expr[order(-v)[seq_len(top_n)], , drop = FALSE]
  }
  pca <- stats::prcomp(t(pca_input), center = TRUE, scale. = FALSE)
  var_explained <- (pca$sdev^2) / sum(pca$sdev^2) * 100

  pca_df <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1L],
    PC2 = pca$x[, 2L],
    group = groups
  )
  # 颜色**不是唯一的语义载体**：形状也编码分组。
  # 色盲读者、以及黑白打印时，仍然能分出两组。
  #
  # ---- 组内置信椭圆 -------------------------------------------------------
  #
  # **不用 `stat_ellipse()`。** 实测它在每组 3 个样本时产出空数据 ——
  # 图上只有点、没有椭圆，ggplot 不报错，坐标范围也没被撑大，
  # 所以从图上完全看不出"这一层没画"。自己算坐标还能落盘核对。
  #
  # 半径按 `sqrt(qchisq(0.95, 2))` = 2.45 倍标准差（`type = "norm"` 的口径），
  # 把协方差当作**已知**。ggplot 默认的 `"t"` 用
  # `sqrt(2 * qf(0.95, 2, n-2))`，每组 3 个样本时是 6.16 倍，椭圆比数据范围还大 6 倍、
  # 把点压成中心一小团，所以不用它。代价是**低估**了小样本下协方差本身的不确定性 ——
  # 这一点必须在副标题里写明，不能让它冒充严格的 95% 区间。
  ell_ok <- all(table(groups) >= 3L)
  ell_radius <- sqrt(stats::qchisq(0.95, 2))
  ell_df <- do.call(rbind, lapply(levels(groups), function(g) {
    idx <- which(as.character(groups) == g)
    e <- if (ell_ok) ellipse_points(pca_df$PC1[idx], pca_df$PC2[idx], level = 0.95) else NULL
    if (is.null(e)) return(NULL)
    e$group <- g
    e
  }))
  if (!is.null(ell_df)) {
    utils::write.csv(ell_df[, c("group", "x", "y", "radius_sd")],
                     file.path(res, "pca_ellipse.csv"), row.names = FALSE)
    ell_extent <- vapply(levels(groups), function(g) {
      d <- ell_df[ell_df$group == g, , drop = FALSE]
      max(diff(range(d$x)), diff(range(d$y)))
    }, numeric(1))
    data_extent <- max(diff(range(pca_df$PC1)), diff(range(pca_df$PC2)))
    log_info(sprintf("PCA 椭圆：半径 %.2f SD；椭圆跨度 %s，数据跨度 %.1f（比值 %.2f）",
                     ell_radius, paste(sprintf("%.1f", ell_extent), collapse = " / "),
                     data_extent, max(ell_extent) / data_extent))
  } else {
    log_warn("PCA 椭圆跳过：每组样本不足 3 个，或协方差奇异")
  }

  # **散点标签放不放得下，也算出来，不能无条件画。**
  # 标签贴在点上方（`vjust = -1`）沿 x 排，所以两个标签碰撞的判据是
  # "x 方向间距 < 标签宽度" —— 还是一维模型，可用长度是**面板宽度**。
  # 与行名/列名的唯一区别：每个标签占的是**它自己的字宽**而不是行高：
  #     max(nchar) × 0.5 × fontsize + min_gap
  # （0.5 × 字号 是 wrap_subtitle() 一直在用的估字宽，这里沿用同一个估值。）
  # 所以把 `pt_per_label` 传给 decide_rownames() 的 `fontsize` 形参 ——
  # 它在这个模型里就是"每个标签占的长度"。下面另打一行日志把真实字号说清楚，
  # 免得日志里那个数被当成字号读。
  # 实测同一段代码两种结果：GSE42568 有 121 个 GSM 号 -> 可容纳约 27 个，
  # **不标**；GSE64790 只有 6 个样本 -> 照常标注。
  pt_fs <- 2.4
  pt_per_label <- max(nchar(as.character(pca_df$sample))) * 0.5 * pt_fs
  show_pt <- decide_rownames(nrow(pca_df), mm(165), pt_per_label, "PCA 样本名",
                             panel_frac = 0.85, min_gap = 2.5,
                             figure = "pca_plot")
  log_info(sprintf("PCA 散点标签：%d 个样本、字号 %gpt、最长标签 %d 字符 -> %s",
                   nrow(pca_df), pt_fs, max(nchar(as.character(pca_df$sample))),
                   if (isTRUE(show_pt)) "标注" else "不标注（放不下，会糊成一片）"))

  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = group, shape = group,
                              label = sample))
  if (!is.null(ell_df)) {
    # 椭圆用自己算的坐标画。fill 映射到 group，图例靠 show.legend = FALSE 压掉。
    p_pca <- p_pca + geom_polygon(
      data = ell_df, ggplot2::aes(x = x, y = y, fill = group, group = group),
      colour = NA, alpha = 0.12, show.legend = FALSE, inherit.aes = FALSE)
  }
  p_pca <- p_pca +
    geom_point(size = 3.6, stroke = 0.9)
  if (isTRUE(show_pt)) {
    p_pca <- p_pca +
      geom_text(vjust = -1, size = pt_fs, show.legend = FALSE, colour = PAL$ink)
  }
  p_pca <- p_pca +
    scale_colour_condition(levels(groups), name = NULL) +
    scale_fill_condition(levels(groups), name = NULL) +
    scale_shape_manual(values = c(16, 17, 15, 18, 8)[seq_along(levels(groups))],
                       name = NULL) +
    labs(title = sprintf("PCA of %s", cfg$dataset_id),
         subtitle = wrap_subtitle(sprintf(
           paste0("top %d variable genes; group encoded by colour AND shape. ",
                  "Shaded ellipse = within-group 95%% normal ellipse (radius %.2f SD). ",
                  "The SMALLEST group has %d samples (%s), so its covariance rests on ",
                  "%d df and the ellipse understates the true spread - it is a visual ",
                  "aid, not a test. %s"),
           nrow(pca_input), ell_radius, min(table(groups)),
           paste(sprintf("%s n=%d", names(table(groups)), as.integer(table(groups))),
                 collapse = ", "),
           min(table(groups)) - 1L,
           # **标不出来就要说出来。** 图上少一层而图注不提，读者会以为
           # "这张图没有样本名"是设计如此，而不是"放不下"。
           if (isTRUE(show_pt)) "Sample accessions are labelled."
           else sprintf(paste0("Sample accessions are NOT labelled: %d labels ",
                               "cannot be placed legibly in this panel width."),
                        nrow(pca_df))),
           fig_width = mm(165)),
         x = sprintf("PC1 (%.1f%% variance)", var_explained[1L]),
         y = sprintf("PC2 (%.1f%% variance)", var_explained[2L])) +
    theme_paper(10)
  save_pdf(file.path(res, "pca_plot.pdf"), print(p_pca), width = mm(165), height = mm(152))
  log_info(sprintf("已生成 pca_plot.pdf（PC1=%.1f%%, PC2=%.1f%%）",
                   var_explained[1L], var_explained[2L]))

  # ---- 4. 样本间相关性 ----------------------------------------------------
  pearson  <- stats::cor(expr, method = "pearson")
  spearman <- stats::cor(expr, method = "spearman")

  # 每个样本与其它样本的最低相关系数；低于阈值判为离群
  min_cor <- apply(pearson, 1L, function(x) min(x[x < 1], na.rm = TRUE))
  outlier <- names(min_cor)[min_cor < cfg$thresholds$outlier_cor]
  if (length(outlier) > 0L) {
    log_warn(sprintf("相关性离群样本 (min Pearson < %.2f): %s",
                     cfg$thresholds$outlier_cor, paste(outlier, collapse = ", ")))
  } else {
    log_info(sprintf("无相关性离群样本（全部 min Pearson >= %.2f）",
                     cfg$thresholds$outlier_cor))
  }

  # ---- 4b. 组内 vs 组间相关性 ---------------------------------------------
  #
  # "min Pearson < 阈值 即离群"隐含一个假设：所有样本是同一组织的技术重复。
  # 当分组本身对应一个巨大的整体表达位移（批次效应，或组织成分差异）时，
  # 每个样本都会与另一组的样本低相关，于是**全部样本都被判为离群**。
  # 这是规则失效，不是 9 个样本都坏了 —— 所以必须把组内/组间拆开报告。
  grp_chr <- as.character(groups)
  ns <- ncol(pearson)
  same_grp <- outer(grp_chr, grp_chr, "==")
  upper <- upper.tri(pearson)
  within_r  <- pearson[which(same_grp & upper, arr.ind = TRUE)]
  between_r <- pearson[which(!same_grp & upper, arr.ind = TRUE)]
  mean_within  <- if (length(within_r)  > 0L) mean(within_r)  else NA_real_
  mean_between <- if (length(between_r) > 0L) mean(between_r) else NA_real_

  # 每个样本自身的组内/组间平均相关
  mw_by_sample <- vapply(seq_len(ns), function(i) {
    j <- setdiff(seq_len(ns), i); j <- j[grp_chr[j] == grp_chr[i]]
    if (length(j) == 0L) NA_real_ else mean(pearson[i, j])
  }, numeric(1))
  mb_by_sample <- vapply(seq_len(ns), function(i) {
    j <- setdiff(seq_len(ns), i); j <- j[grp_chr[j] != grp_chr[i]]
    if (length(j) == 0L) NA_real_ else mean(pearson[i, j])
  }, numeric(1))

  # 组间平均相关低于离群阈值 => 分组与一个全局表达位移混杂，
  # tumor-vs-normal 的差异里分不清多少是分组、多少是这个位移
  confounded <- is.finite(mean_within) && is.finite(mean_between) &&
    mean_between < cfg$thresholds$outlier_cor
  log_info(sprintf("组内平均 Pearson = %.3f，组间平均 Pearson = %.3f",
                   mean_within, mean_between))
  if (confounded) {
    log_warn(sprintf(paste0(
      "组间平均相关 %.3f 低于离群阈值 %.2f：分组与一个全局表达位移高度混杂。",
      "此时 tumor-vs-normal 的差异无法与批次/组织成分差异分离，",
      "limma 结果只能作为假设生成，不能当作肿瘤特异事件。"),
      mean_between, cfg$thresholds$outlier_cor))
  }

  # ---- 5. 样本间相关性热图 ------------------------------------------------
  # 相关性全部在 0.9+ 区间，是**单向**量（越高越好），所以用序列色而不是发散色。
  # 用发散色会让人以为 0.5 是"中性"中点，而这里根本没有负相关。
  annotation_col <- data.frame(group = groups, row.names = colnames(expr))
  annotation_colors <- list(group = gp)

  # **格子里那个数字放不放得下，同样要算。**
  # `display_numbers = TRUE` 会给**每一个**格子画一个数：121×121 就是 14641 个
  # 6pt 的 "0.987"。183 mm 宽 ÷ 121 列 = 每格 1.5 mm，而 "0.987" 要 6 mm ——
  # 这**不是"对比度差"，是一个数都读不出来**：所有数字叠在一起，
  # 图上表现为一片灰糊，而"糊了"从图注上完全看不出来。
  # 判据沿用同一个一维模型，可用长度取热图面板宽（约 0.62 × 整幅）：
  #   每个数占 = nchar(number_format) × 0.5 × fontsize + min_gap
  num_fs  <- 6
  num_len <- nchar(sprintf("%.3f", 0)) * 0.5 * num_fs
  show_num <- decide_rownames(ncol(pearson), W_DOUBLE, num_len, "相关性热图格内数值",
                              panel_frac = 0.62, min_gap = 2.5,
                              figure = "correlation_heatmap")
  # 行列名同理：旋转 90°，横向占用的正是它的行高 `fontsize + min_gap`。
  cor_fs   <- 5
  show_cn2 <- decide_rownames(ncol(pearson), W_DOUBLE, cor_fs, "相关性热图样本名",
                              panel_frac = 0.62, min_gap = 2.5,
                              figure = "correlation_heatmap")
  log_info(sprintf(paste0("相关性热图：%d 个样本 -> 格内数值 %s（字号 %gpt，",
                          "每个数占 %.1fpt）、行列名 %s（字号 %gpt）"),
                   ncol(pearson),
                   if (isTRUE(show_num)) "显示" else "隐藏（放不下）", num_fs, num_len,
                   if (isTRUE(show_cn2)) "显示" else "隐藏（放不下）", cor_fs))

  save_pdf(file.path(res, "correlation_heatmap.pdf"), {
    pheatmap::pheatmap(
      pearson,
      annotation_col = annotation_col,
      annotation_row = annotation_col,
      annotation_colors = annotation_colors,
      display_numbers = isTRUE(show_num), number_format = "%.3f",
      fontsize_number = num_fs,
      show_rownames = isTRUE(show_cn2), show_colnames = isTRUE(show_cn2),
      fontsize = cor_fs,
      color = pal_sequential(100),
      border_color = "white", treeheight_row = 18, treeheight_col = 18,
      # **少画了哪一层，标题上要写出来。** pheatmap 没有副标题，
      # 不说的话读者会以为这张热图本来就不带数值。
      main = sprintf("Sample-sample Pearson correlation (%s)%s", cfg$dataset_id,
                     if (isTRUE(show_num)) ""
                     else " - cell values omitted (they do not fit); see correlation_matrix.csv"),
      silent = FALSE
    )
    # pheatmap 的色条固定在图右侧，没有位置参数可调。
    # 它是**细长条**，占的宽度远小于 ggplot 的右侧图例，所以这张图保留右侧。
  }, width = W_DOUBLE, height = mm(165))
  log_info("已生成 correlation_heatmap.pdf")

  # ---- 6. 相关性矩阵落盘 --------------------------------------------------
  cor_df <- data.frame(
    sample = rownames(pearson),
    group = as.character(groups),
    min_pearson = round(min_cor[rownames(pearson)], 4),
    mean_within_group = round(mw_by_sample, 4),
    mean_between_group = round(mb_by_sample, 4),
    outlier = rownames(pearson) %in% outlier,
    pearson, spearman,
    check.names = FALSE
  )
  utils::write.csv(cor_df, file.path(res, "correlation_matrix.csv"), row.names = FALSE)

  write_json(file.path(res, "qc_summary.json"), list(
    dataset_id = cfg$dataset_id,
    genes = nrow(expr),
    samples = ncol(expr),
    pc1_variance = round(var_explained[1L], 2),
    pc2_variance = round(var_explained[2L], 2),
    pc3_variance = if (length(var_explained) >= 3L) round(var_explained[3L], 2) else NA,
    outlier_samples = outlier,
    outlier_threshold = cfg$thresholds$outlier_cor,
    median_min_pearson = round(stats::median(min_cor), 4),
    mean_within_group_pearson = round(mean_within, 4),
    mean_between_group_pearson = round(mean_between, 4),
    group_confounded_with_global_shift = confounded
  ))

  log_info(sprintf("已生成 correlation_matrix.csv（离群样本 %d 个）", length(outlier)))
  invisible(list(pearson = pearson, spearman = spearman, pca = pca, var_explained = var_explained))
}

if (!GEO_ORCHESTRATED()) {
  cfg <- load_config()
  run_02_qc_pca_correlation(cfg)
}
