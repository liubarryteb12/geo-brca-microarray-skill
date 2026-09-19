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

#' 统一的分组配色
group_palette <- function(groups) {
  lv <- sort(unique(groups))
  base <- c("#2E5FA3", "#C1443C", "#3E8E5A", "#8A6BBE", "#D08C34", "#4C9BB5")
  stats::setNames(rep(base, length.out = length(lv)), lv)
}

run_02_qc_pca_correlation <- function(cfg) {
  log_info("=== 步骤 02：QC / PCA / 样本相关性 ===")
  ensure_dirs(cfg)

  res <- cfg$output$results_dir
  expr_raw <- readRDS(file.path(cfg$output$data_dir, "expr_raw.rds"))
  expr <- readRDS(file.path(cfg$output$data_dir, "expr_clean.rds"))
  group <- utils::read.csv(file.path(cfg$output$data_dir, "group.csv"), stringsAsFactors = FALSE)
  groups <- factor(group$group, levels = unique(group$group))

  # ---- 1. 标准化前后箱线图 ------------------------------------------------
  long <- rbind(
    cbind(matrix_to_long(expr_raw), stage = "before"),
    cbind(matrix_to_long(expr), stage = "after")
  )
  long$stage <- factor(long$stage, levels = c("before", "after"))

  p_box <- ggplot(long, aes(x = sample, y = expression, fill = stage)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.25) +
    facet_wrap(~stage, ncol = 2, scales = "free_y") +
    labs(title = "Expression distribution before / after quantile normalization",
         subtitle = sprintf("%s - %d genes x %d samples", cfg$dataset_id, nrow(expr), ncol(expr)),
         x = NULL, y = "log2 expression") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          legend.position = "none")
  save_pdf(file.path(res, "boxplot_before_after.pdf"), print(p_box), width = 10, height = 6)
  log_info("已生成 boxplot_before_after.pdf")

  # ---- 2. 密度曲线 --------------------------------------------------------
  p_density <- ggplot(long, aes(x = expression, colour = sample)) +
    geom_density(linewidth = 0.4) +
    facet_wrap(~stage, ncol = 1, scales = "free_y") +
    labs(title = "Expression density before / after normalization", x = "log2 expression", y = "density") +
    theme_bw(base_size = 9) +
    theme(legend.position = "none")
  save_pdf(file.path(res, "density_plot.pdf"), print(p_density), width = 8, height = 7)
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
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = group, label = sample)) +
    geom_point(size = 3.2) +
    geom_text(vjust = -1, size = 2.4, show.legend = FALSE) +
    scale_colour_manual(values = group_palette(groups)) +
    labs(title = sprintf("PCA of %s", cfg$dataset_id),
         subtitle = sprintf("top %d variable genes", nrow(pca_input)),
         x = sprintf("PC1 (%.1f%% variance)", var_explained[1L]),
         y = sprintf("PC2 (%.1f%% variance)", var_explained[2L])) +
    theme_bw(base_size = 10)
  save_pdf(file.path(res, "pca_plot.pdf"), print(p_pca), width = 7, height = 6)
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
  annotation_col <- data.frame(group = groups, row.names = colnames(expr))
  save_pdf(file.path(res, "correlation_heatmap.pdf"), {
    pheatmap::pheatmap(
      pearson,
      annotation_col = annotation_col,
      annotation_row = annotation_col,
      display_numbers = TRUE, number_format = "%.3f", fontsize_number = 6,
      color = grDevices::colorRampPalette(c("#2E5FA3", "white", "#C1443C"))(100),
      main = sprintf("Sample-sample Pearson correlation (%s)", cfg$dataset_id),
      silent = FALSE
    )
  }, width = 8, height = 7)
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
