# ============================================================================
# 06_wgcna.R — 加权基因共表达网络分析（WGCNA）
# ============================================================================
# 只在 design_mode: cohort 下运行。small_sample 下**记录原因后跳过**，
# 不报错、也不产出一张没有意义的模块图。
#
# 为什么 n < 15 不做：WGCNA 的输入是基因两两相关矩阵，样本量决定相关系数的
# 稳定性。n=6 时观测到 r=0.8 的 95% 置信区间是 [-0.03, +0.98] —— 下界为负，
# 意思是"这 6 个样本根本无法把 r=0.8 和 r=0 区分开"。模块划分建立在这种
# 相关矩阵上，得到的模块不可复现，画出来只会让人误以为发现了结构。
#
# 为什么只用肿瘤组：若把 104 例癌和 17 例正常一起放进去，**第一个模块必然是
# "肿瘤 vs 正常"轴** —— 组织成分差异会盖过肿瘤内部的异质性。而"肿瘤 vs 正常"
# 这件事 step 03 的 DEG 已经答过了，WGCNA 要回答的是另一个问题：
# 在癌组织内部，哪些基因协同变化、这些模块与临床性状如何关联。
# 这个选择会写进 wgcna_status.json，不是隐含假设。
#
# 输出：results/<GSE>/wgcna_modules.csv, wgcna_module_trait.csv,
#       wgcna_soft_power.csv, wgcna_module_sizes.csv,
#       01-06-01-unit1-wgcna-scale-free-fit.pdf, 01-06-02-unit1-wgcna-module-trait-heatmap.pdf, wgcna_status.json
# ============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
})

# ---- bootstrap -------------------------------------------------------------
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

WGCNA_MIN_SAMPLES <- 15L    # 与 00_validate_inputs.R 的 COHORT_MIN_SAMPLES 一致
WGCNA_MIN_GENES   <- 2000L  # 基因太少时共表达网络没有意义

#' 从 clinical.csv 里挑出可作性状的字段
#'
#' 规则是显式的：
#'   * 数值字段：非缺失 >= 80%、取值 >= 2 种 -> 直接用
#'   * 两水平分类字段：非缺失 >= 80% -> 编码成 0/1（按字母序，第一级为 0）
#'   * 三水平及以上：跳过。**不做哑变量展开** —— 一个三水平因子展开成两列
#'     再和模块算相关，得到的是"该水平 vs 其余"的混合对比，解释起来很容易出错。
#'
#' @return list(traits = data.frame, encoded = 命名字符向量, skipped = 命名字符向量)
build_traits <- function(clinical, sample_ids) {
  clin <- clinical[match(sample_ids, clinical$gsm), , drop = FALSE]
  traits <- list(); encoded <- character(0); skipped <- character(0)
  for (k in setdiff(colnames(clinical), "gsm")) {
    v <- clin[[k]]
    present <- !is.na(v) & nzchar(v)
    if (sum(present) < 0.8 * length(v)) { skipped[k] <- "coverage < 80%"; next }
    num <- suppressWarnings(as.numeric(v[present]))
    if (!anyNA(num) && length(unique(num)) >= 2L) {
      out <- rep(NA_real_, length(v)); out[present] <- num
      traits[[k]] <- out; encoded[k] <- "numeric as-is"; next
    }
    lv <- sort(unique(v[present]))
    if (length(lv) == 2L) {
      out <- rep(NA_real_, length(v)); out[present] <- as.numeric(v[present] == lv[2L])
      traits[[k]] <- out
      encoded[k] <- sprintf("binary: %s=0, %s=1", lv[1L], lv[2L])
      next
    }
    skipped[k] <- sprintf("%d levels (only numeric or binary are used)", length(lv))
  }
  list(traits = as.data.frame(traits, stringsAsFactors = FALSE),
       encoded = encoded, skipped = skipped)
}

#' 选软阈值 power
#'
#' 判据：取**最小的**使 scale-free topology R² >= 0.8 的 power。
#' 若没有任何 power 达到 0.8，取 R² 最大的那个并**如实记录** ——
#' 这时网络不是无标度的，模块结果要打折看，不能假装通过。
#'
#' @return list(power, r2, reached, table)
pick_power <- function(datExpr, powers = c(1:10, seq(12, 20, by = 2)), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers,
                                  networkType = "signed", verbose = 0)
  tab <- data.frame(power = sft$fitIndices[, 1],
                    r2 = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
                    slope = sft$fitIndices[, 3],
                    mean_k = sft$fitIndices[, 5],
                    median_k = sft$fitIndices[, 6])
  # **联合判据**（PLAN-T-W1，文献规范）：R2>=0.8 且 mean connectivity >= 1。
  # 只看 R2 会选中 mean_k 趋 0 的 power —— 网络碎成孤点，模块无意义。
  # 不存在两者都满足的 power 时退回 R2 判据，并如实记录。
  ok_both <- which(tab$r2 >= 0.8 & tab$mean_k >= 1)
  # **联合判据**（PLAN-T-W1，文献规范）：R2>=0.8 且 mean connectivity >= 1。
  # 只看 R2 会选中 mean_k 趋 0 的 power —— 网络碎成孤点，模块无意义。
  ok_both <- which(tab$r2 >= 0.8 & tab$mean_k >= 1)
  ok <- which(tab$r2 >= 0.8)
  if (length(ok_both) > 0L) {
    i <- ok_both[1L]
    list(power = as.integer(tab$power[i]), r2 = tab$r2[i], reached = TRUE,
         joint_criterion = TRUE, mean_k = tab$mean_k[i], table = tab)
  } else if (length(ok) > 0L) {
    i <- ok[1L]
    list(power = as.integer(tab$power[i]), r2 = tab$r2[i], reached = TRUE,
         joint_criterion = FALSE, mean_k = tab$mean_k[i], table = tab)
  } else {
    i <- which.max(tab$r2)
    list(power = as.integer(tab$power[i]), r2 = tab$r2[i], reached = FALSE,
         joint_criterion = FALSE, mean_k = tab$mean_k[i], table = tab)
  }
}

#' 软阈值诊断图 —— **两张单图，不是一张双面板图**
#'
#' 两个面板问的是**两个不同的问题**、y 轴含义也不同（R² vs 平均连通性），
#' 所以按仓库的出图约定拆成 `unit1` / `unit2` 两张单图。
#'
#' 原来的双面板图还有一个隐患：`scales = "free_y"` 下左右两边的 y 轴刻度
#' **不一样**，读者很容易把两条曲线按同一把尺子比 —— 而副标题里写着
#' "Left/Right"，一旦面板顺序或布局变了，文字与图就对不上。
#' 拆开之后各自的 y 轴含义直接写在轴上，不需要靠文字指认左右。
#'
#' （原来还专门为"facet_wrap 按字母序排面板会让左右与文字相反"钉过 levels，
#' 拆成单图之后那一类问题从根上没有了。）
make_soft_power_panel <- function(tab, chosen, cfg, panel) {
  if (identical(panel, "scale-free topology fit")) {
    d <- data.frame(power = tab$power, value = tab$r2)
    ylab <- "scale-free topology fit R2"
    title <- sprintf("WGCNA scale-free topology fit - %s", cfg$dataset_id)
    sub <- sprintf(paste0("dashed line = chosen power %d. Target R2 >= 0.8; ",
                          "the chosen power is the smallest one reaching it."), chosen)
  } else {
    d <- data.frame(power = tab$power, value = tab$mean_k)
    ylab <- "mean connectivity"
    title <- sprintf("WGCNA mean connectivity - %s", cfg$dataset_id)
    sub <- sprintf(paste0("dashed line = chosen power %d. Mean connectivity ",
                          "must stay above 0."), chosen)
  }
  ggplot2::ggplot(d, ggplot2::aes(x = power, y = value)) +
    ggplot2::geom_line(colour = PAL$muted, linewidth = 0.4) +
    ggplot2::geom_point(colour = PAL$ink, size = 1.4) +
    ggplot2::geom_vline(xintercept = chosen, linetype = "dashed",
                        linewidth = 0.4, colour = PAL$up) +
    ggplot2::labs(title = title,
                  subtitle = wrap_subtitle(sub, fig_width = W_ONE_HALF),
                  x = "soft-thresholding power", y = ylab) +
    theme_paper(10)
}

#' 模块-性状关联图
#'
#' 格子里写 r 值本身，**不靠颜色单独承载信息** —— 颜色只表示相关方向与强度，
#' 具体数值和显著性必须能直接读出来。
make_module_trait_plot <- function(cor_df, cfg) {
  cor_df$module <- factor(cor_df$module, levels = rev(unique(cor_df$module)))
  cor_df$trait  <- factor(cor_df$trait,  levels = unique(cor_df$trait))
  cor_df$label  <- sprintf("%.2f%s", cor_df$cor,
                           ifelse(cor_df$p_adj < 0.001, "***",
                                  ifelse(cor_df$p_adj < 0.01, "**",
                                         ifelse(cor_df$p_adj < 0.05, "*", ""))))
  ggplot2::ggplot(cor_df, ggplot2::aes(x = trait, y = module, fill = cor)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 2.4, colour = PAL$ink) +
    ggplot2::scale_fill_gradientn(colours = pal_diverging(100), limits = c(-1, 1),
                                  name = "Pearson r") +
    ggplot2::labs(
      title = sprintf("WGCNA module-trait relationships - %s", cfg$dataset_id),
      subtitle = wrap_subtitle(sprintf(
        paste0("%d modules x %d traits on %d tumour samples. Cell text = Pearson r; ",
               "* p_adj<0.05, ** <0.01, *** <0.001 (BH across all %d tests). ",
               "Modules are named by WGCNA colour labels; the colour here encodes r only."),
        length(unique(cor_df$module)), length(unique(cor_df$trait)),
        cor_df$n[1L], nrow(cor_df)), fig_width = W_DOUBLE),
      x = NULL, y = NULL) +
    theme_paper(10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 8),
      panel.grid  = ggplot2::element_blank())
}

run_06_wgcna <- function(cfg) {
  log_info("=== 步骤 06：WGCNA 共表达网络 ===")
  ensure_dirs(cfg)
  res <- cfg$output$results_dir
  dat <- cfg$output$data_dir
  status <- list(step = "wgcna")

  # ---- 门禁 1：设计模式 ---------------------------------------------------
  if (!identical(cfg$design_mode, "cohort")) {
    status$status <- "not_applicable"
    status$reason <- sprintf(
      paste0("design_mode=%s。WGCNA 的输入是基因两两相关矩阵，样本量决定相关系数的稳定性；",
             "n<15 时 r=0.8 的 95%% 置信区间下界为负，模块划分不可复现。"),
      cfg$design_mode)
    write_json(file.path(res, "wgcna_status.json"), status)
    log_warn(sprintf("跳过 WGCNA：%s", status$reason))
    return(invisible(NULL))
  }

  # ---- 门禁 2：包可用 -----------------------------------------------------
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    status$status <- "package_missing"
    status$reason <- "WGCNA 未安装（见 workflow 的 packages 列表）"
    write_json(file.path(res, "wgcna_status.json"), status)
    log_warn(sprintf("跳过 WGCNA：%s", status$reason))
    return(invisible(NULL))
  }

  # ---- 1. 取肿瘤组样本 ----------------------------------------------------
  expr  <- readRDS(file.path(dat, "expr_clean.rds"))
  group <- utils::read.csv(file.path(dat, "group.csv"), stringsAsFactors = FALSE)
  clinical_path <- file.path(dat, "clinical.csv")
  clinical <- if (file.exists(clinical_path)) {
    utils::read.csv(clinical_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    NULL
  }

  tumor_arm <- cfg$contrast[1L]
  keep <- group$gsm[group$group == tumor_arm]
  keep <- intersect(keep, colnames(expr))
  # **验证专用：样本截断开关**（PLAN-T-W1 L2，用户指示用小样本快速验证代码）。
  # 设了 `wgcna_sample_cap` 时只取前 N 个肿瘤样本 —— 让 WGCNA 在 12-20 样本上
  # 跑完（pickSoftThreshold + blockwiseModules 从分钟级降到秒级），
  # 目的是**验证代码与出图设置**，不是产出结论。全样本跑时不设这个字段。
  # 截断是显式的、进状态文件的 —— 不能悄悄少跑样本。
  cap <- cfg$analysis$wgcna_sample_cap
  if (!is.null(cap) && is.numeric(cap) && length(keep) > cap) {
    keep <- keep[seq_len(as.integer(cap))]
    status$sample_cap <- as.integer(cap)
    status$sample_cap_note <- paste0("验证模式：样本截断到 ", cap,
      " 个（wgcna_sample_cap），用于快速验证代码与出图，结果不可当结论")
    log_warn(status$sample_cap_note)
  }
  if (length(keep) < WGCNA_MIN_SAMPLES) {
    status$status <- "not_applicable"
    status$reason <- sprintf("组 %s 只有 %d 个样本（< %d）", tumor_arm, length(keep), WGCNA_MIN_SAMPLES)
    write_json(file.path(res, "wgcna_status.json"), status)
    log_warn(sprintf("跳过 WGCNA：%s", status$reason))
    return(invisible(NULL))
  }

  # ---- 2. 选高变基因 ------------------------------------------------------
  # WGCNA 是全基因两两相关，全基因组跑既慢又没有额外信息：
  # 低变基因的相关基本是噪声，还会把软阈值拟合带偏。
  top_genes <- cfg$analysis$wgcna_top_genes %||% 5000L
  v <- apply(expr[, keep, drop = FALSE], 1L, stats::var)
  v[!is.finite(v)] <- 0
  if (length(v) < WGCNA_MIN_GENES) {
    status$status <- "not_applicable"
    status$reason <- sprintf("表达矩阵只有 %d 个基因（< %d），共表达网络没有意义",
                             length(v), WGCNA_MIN_GENES)
    write_json(file.path(res, "wgcna_status.json"), status)
    log_warn(sprintf("跳过 WGCNA：%s", status$reason))
    return(invisible(NULL))
  }
  n_top <- min(as.integer(top_genes), length(v))
  sel <- names(sort(v, decreasing = TRUE))[seq_len(n_top)]

  # WGCNA 要的是 样本 x 基因
  datExpr <- t(expr[sel, keep, drop = FALSE])
  log_info(sprintf("WGCNA 输入: %d 个样本 x %d 个高变基因（组 %s）",
                   nrow(datExpr), ncol(datExpr), tumor_arm))

  # ---- 3. 样本/基因级 QC --------------------------------------------------
  # goodSamplesGenes 会挑出全 NA、方差为零、以及表达值离群的基因。
  # **不静默丢** —— 丢了多少要记下来，否则"用了多少基因"对不上。
  gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0L) {
      log_warn(sprintf("WGCNA: 剔除 %d 个坏基因（全 NA / 零方差 / 离群）", sum(!gsg$goodGenes)))
    }
    if (sum(!gsg$goodSamples) > 0L) {
      log_warn(sprintf("WGCNA: 剔除 %d 个坏样本", sum(!gsg$goodSamples)))
    }
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  }

  # ---- 3b. 样本层次聚类树（传统 WGCNA 的离群检测，差距 PLAN-T-W1）----
  # 文献标准做法：对样本做 hclust(dist(cor(datExpr), method="average"))
  # 画树并检查是否有树高显著偏离主体的离群样本。
  # **只标记 + 记录，不自动剔除** —— 剔除与否是人工复核节点
  # outlier_removal 的职责（规范 §8），脚本不能替人做决定。
  sample_dist <- 1 - stats::cor(t(datExpr), use = "pairwise.complete.obs")
  sample_tree <- stats::hclust(as.dist(sample_dist), method = "average")
  # TreeHeight 判据（FAQ 推荐）：主体树高的 1.5 倍以上 = 离群候选
  merge_h <- sample_tree$height
  body_h <- if (length(merge_h) >= 2L) stats::median(merge_h) else 0
  outl <- which(merge_h > 1.5 * body_h)
  status$sample_dendrogram <- list(
    n_samples = nrow(datExpr),
    tree_height_max = round(max(merge_h), 3),
    outlier_candidates = if (length(outl) > 0L)
      as.list(rownames(datExpr)[unique(outl)]) else list(),
    note = "候选 = 合并高度 > 1.5 x 中位合并高。是否剔除由人工复核节点 outlier_removal 决定，脚本不自动删。")
  grp_lv <- unique(group$group)
  grp_pal <- stats::setNames(c(PAL$up, PAL$down, PAL$primary, PAL$muted), grp_lv)
  group_col <- unname(grp_pal[group$group[match(rownames(datExpr), group$gsm)]])
  png(file.path(res, "01-06-03-unit1-sample-dendrogram.png"),
      width = W_DOUBLE, height = mm(80), units = "in", res = 300)
  WGCNA::plotDendroAndColors(sample_tree, group_col, "Group (colours = group)",
                             dendroLabels = FALSE, hang = 0.03,
                             addGuide = TRUE, guideHang = 0.05,
                             main = "Sample dendrogram (outlier check)",
                             cex.labels = 0.4)
  dev.off()
  pdf(file.path(res, "01-06-03-unit1-sample-dendrogram.pdf"),
      width = W_DOUBLE, height = mm(80))
  WGCNA::plotDendroAndColors(sample_tree, group_col, "Group (colours = group)",
                             dendroLabels = FALSE, hang = 0.03,
                             addGuide = TRUE, guideHang = 0.05,
                             main = "Sample dendrogram (outlier check)",
                             cex.labels = 0.4)
  dev.off()
  if (length(outl) > 0L) {
    log_warn(sprintf("WGCNA: %d 个离群候选样本（见 01-06-03-unit1-sample-dendrogram）—— 不自动剔除",
                     length(outl)))
  }
  status$n_samples <- nrow(datExpr)
  status$n_genes_input <- ncol(datExpr)
  status$good_samples_genes_all_ok <- gsg$allOK

  seed <- cfg$analysis$seed
  if (!is.null(seed)) set.seed(seed)

  # ---- 4. 软阈值 ----------------------------------------------------------
  pw <- pick_power(datExpr, seed = seed)
  status$soft_power <- pw$power
  status$soft_power_r2 <- round(pw$r2, 4)
  status$soft_power_target_reached <- pw$reached
  if (!pw$reached) {
    # 如实记录：没达到 0.8 时网络不是无标度的，模块结果要打折看
    log_warn(sprintf("WGCNA: 没有任何 power 使 R2 >= 0.8（最大 %.3f，取 power=%d）；网络非无标度，模块结果需谨慎",
                     pw$r2, pw$power))
  } else {
    log_info(sprintf("WGCNA: 软阈值 power=%d（R2=%.3f >= 0.8）", pw$power, pw$r2))
  }
  utils::write.csv(pw$table, file.path(res, "wgcna_soft_power.csv"), row.names = FALSE)
  save_pdf(file.path(res, "01-06-01-unit1-wgcna-scale-free-fit.pdf"),
           print(make_soft_power_panel(pw$table, pw$power, cfg,
                                       "scale-free topology fit")),
           width = W_ONE_HALF, height = mm(80))
  save_pdf(file.path(res, "01-06-01-unit2-wgcna-mean-connectivity.pdf"),
           print(make_soft_power_panel(pw$table, pw$power, cfg,
                                       "mean connectivity")),
           width = W_ONE_HALF, height = mm(80))

  # ---- 5. 模块识别 --------------------------------------------------------
  #
  # **这是本仓库唯一一处 `library()` 调用，必须临时挂载。**
  #
  # `blockwiseModules` 内部用 `do.call(corFnc, ...)` 计算 KME，而 `corFnc` 来自
  # 包内常量 `.corFnc = c("cor", "bicor", "cor")` —— 是个**字符串**，按名字查找。
  # 本仓库不 attach 任何包，于是它解析到 `stats::cor`，而后者没有
  # `weights.x` / `weights.y` / `cosine` 参数。实测报错：
  #   unused arguments (weights.x = NULL, weights.y = NULL, cosine = FALSE)
  #
  # **试过传 `corFnc = WGCNA::cor`，没用。** 读了 WGCNA 1.74 源码：
  # `blockwiseModules` 的形参表里根本没有 `corFnc`（只有 `corType`），
  # 我的参数掉进 `...`；而 KME 那段用的是包内常量，不看 `...`。报错一字不变。
  # 这是包内部按字符串查函数的行为，**从外面没有任何参数能改**。
  #
  # 所以只能在调用期间把 WGCNA 挂到搜索路径上，让它按名字能查到自己的 `cor`。
  # `on.exit` 立刻 detach，遮蔽窗口仅限这一次调用；本脚本其余所有调用
  # （包括下面的 `moduleEigengenes` / `labels2colors`）仍然写全名。
  # 挂载前先记下它是否已经在搜索路径上，避免把调用方原有的状态拆掉。
  wgcna_attached <- "package:WGCNA" %in% search()
  if (!wgcna_attached) {
    # **`library("WGCNA")`，不能写 `library(WGCNA, character.only = TRUE)`。**
    # 后者会去**求值** `WGCNA` 这个符号，而它当然不存在 ——
    # 实测报 `object 'WGCNA' not found`，而那个报错看着像"包没装"，
    # 完全指不到真正的原因。character.only=TRUE 是给"包里存了个变量名"用的，
    # 这里直接给字符串就行。
    suppressPackageStartupMessages(library("WGCNA", warn.conflicts = FALSE))
    on.exit({
      try(detach("package:WGCNA", unload = FALSE), silent = TRUE)
    }, add = TRUE)
  }

  # randomSeed 必须显式传：blockwiseModules 内部有随机初始化，
  # 不传的话同一份输入两轮给出不同模块（见 AGENTS.md 规则 11）。
  min_mod <- as.integer(cfg$analysis$wgcna_min_module_size %||% 30L)
  net <- WGCNA::blockwiseModules(
    datExpr, power = pw$power, networkType = "signed", TOMType = "signed",
    minModuleSize = min_mod, mergeCutHeight = 0.25,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    randomSeed = if (is.null(seed)) 12345L else as.integer(seed),
    verbose = 0)

  module_label <- net$colors
  mod_names <- WGCNA::labels2colors(module_label)
  status$n_modules <- length(unique(module_label))
  status$n_modules_nongrey <- sum(unique(module_label) != 0L)
  status$min_module_size <- min_mod
  # **状态在这里就定成 ok。** 后面每一段（没有 clinical.csv、没有可用性状、
  # 只有 grey 模块）都是"网络建成了、只是关联做不了"，它们会覆写
  # `module_trait` 字段说明原因，但不该把整个步骤降级成"没跑"——
  # 验收项看的是 `status$status`，漏了这一行会让一个真的跑完的 WGCNA
  # 在验收日志里显示成 FAIL。
  status$status <- "ok"
  log_info(sprintf("WGCNA: 识别出 %d 个模块（其中 %d 个非 grey），最小模块 %d 个基因",
                   status$n_modules, status$n_modules_nongrey, min_mod))

  mod_df <- data.frame(gene = colnames(datExpr),
                       module_label = module_label,
                       module = mod_names,
                       stringsAsFactors = FALSE)
  mod_df <- mod_df[order(mod_df$module_label, mod_df$gene), , drop = FALSE]
  utils::write.csv(mod_df, file.path(res, "wgcna_modules.csv"), row.names = FALSE)

  sizes <- as.data.frame(table(module = mod_df$module), stringsAsFactors = FALSE)
  colnames(sizes) <- c("module", "n_genes")
  sizes <- sizes[order(-sizes$n_genes), , drop = FALSE]
  utils::write.csv(sizes, file.path(res, "wgcna_module_sizes.csv"), row.names = FALSE)

  # ---- 6. 模块-性状关联 ---------------------------------------------------
  if (is.null(clinical) || nrow(clinical) == 0L) {
    status$module_trait <- "skipped: 没有 clinical.csv"
    write_json(file.path(res, "wgcna_status.json"), status)
    log_warn("WGCNA: 没有 clinical.csv，跳过模块-性状关联")
    return(invisible(mod_df))
  }

  bt <- build_traits(clinical, rownames(datExpr))
  status$traits_used <- names(bt$traits)
  status$traits_encoding <- as.list(bt$encoded)
  status$traits_skipped <- as.list(bt$skipped)
  log_info(sprintf("WGCNA 性状: %d 个可用（%s）", ncol(bt$traits),
                   paste(names(bt$traits), collapse = ", ")))
  if (length(bt$skipped) > 0L) {
    log_info(sprintf("WGCNA 性状跳过 %d 个: %s", length(bt$skipped),
                     paste(sprintf("%s (%s)", names(bt$skipped), unlist(bt$skipped)),
                           collapse = "; ")))
  }
  if (ncol(bt$traits) == 0L) {
    status$module_trait <- "skipped: 没有可用性状"
    write_json(file.path(res, "wgcna_status.json"), status)
    log_warn("WGCNA: 没有可用性状，跳过模块-性状关联")
    return(invisible(mod_df))
  }

  me <- WGCNA::moduleEigengenes(datExpr, colors = module_label)$eigengenes

  # **kME 落盘**（PLAN-T-W1，文献标准）：kME = cor(gene, ME) 即
  # module membership，hub 基因的排序依据；不落盘等于没算。
  kme <- stats::cor(datExpr, me, use = "pairwise.complete.obs")
  utils::write.csv(data.frame(gene = colnames(datExpr), kme, check.names = FALSE),
                   file.path(res, "wgcna_kme.csv"), row.names = FALSE)
  hub_top <- lapply(colnames(kme), function(m) as.list(utils::head(sort(kme[, m], decreasing = TRUE), 5L)))
  names(hub_top) <- colnames(kme)
  status$kme_top5 <- hub_top

  gsm_plots <- 0L
  # grey 模块是"未分配"，它的特征基因没有生物学含义，不参与关联
  me <- me[, colnames(me) != "ME0", drop = FALSE]
  if (ncol(me) == 0L) {
    status$module_trait <- "skipped: 只有 grey 模块"
    write_json(file.path(res, "wgcna_status.json"), status)
    log_warn("WGCNA: 只有 grey 模块，跳过模块-性状关联")
    return(invisible(mod_df))
  }

  rows <- list()
  for (m in colnames(me)) {
    for (t in colnames(bt$traits)) {
      x <- me[[m]]; y <- bt$traits[[t]]
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 3L) next
      ct <- stats::cor.test(x[ok], y[ok], method = "pearson")
      rows[[length(rows) + 1L]] <- data.frame(
        module = sub("^ME", "", m), trait = t,
        cor = unname(ct$estimate), p = ct$p.value, n = sum(ok),
        stringsAsFactors = FALSE)
    }
  }
  if (length(rows) == 0L) {
    status$module_trait <- "skipped: 没有可算的模块-性状对"
    write_json(file.path(res, "wgcna_status.json"), status)
    return(invisible(mod_df))
  }
  cor_df <- do.call(rbind, rows)
  # **多重检验校正。** 模块 x 性状是几十到上百次检验，不校正的话
  # p<0.05 的格子会有一堆是偶然的，而热图上它们和真信号长得一样。
  cor_df$p_adj <- stats::p.adjust(cor_df$p, method = "BH")
  cor_df <- cor_df[order(cor_df$p_adj), , drop = FALSE]
  rownames(cor_df) <- NULL
  utils::write.csv(cor_df, file.path(res, "wgcna_module_trait.csv"), row.names = FALSE)

  n_sig <- sum(cor_df$p_adj < 0.05)
  status$module_trait <- "ok"
  status$module_trait_tests <- nrow(cor_df)
  status$module_trait_significant <- n_sig
  status$module_trait_correction <- "BH"
  log_info(sprintf("WGCNA 模块-性状: %d 对检验，BH 校正后 %d 对 p_adj < 0.05",
                   nrow(cor_df), n_sig))

  save_pdf(file.path(res, "01-06-02-unit1-wgcna-module-trait-heatmap.pdf"),
           print(make_module_trait_plot(cor_df, cfg)),
           width = W_DOUBLE, height = max(mm(81), 0.32 * length(unique(cor_df$module)) + 1.6))

  # ---- GS-MM 散点（PLAN-T-W1，文献核心工具）----
  # 整段包 tryCatch：这一段有独立的失败模式，不能让它拖垮前面的产物。
  # 错误详情进 status$gsmm_error（可见、可查），不当"没做"处理。
  status$gsmm_error <- tryCatch({ NULL
  # 每个 trait 选 |cor| 最高的 module，画 GS vs MM 散点，
  # hub 候选（top2% both）标基因名；Spearman rho 标在图上。
  DYNAMIC_FIG_BASES_DECL = '04:8'
  GSMM_BASE <- paste0(as.character(1), "-06-04-unit")
  for (t in colnames(bt$traits)) {
    sub <- cor_df[cor_df$trait == t & is.finite(cor_df$cor), , drop = FALSE]
    if (nrow(sub) == 0L) next
    best_m <- sub$module[which.max(abs(sub$cor))]
    me_name <- paste0("ME", best_m)
    if (!me_name %in% colnames(kme)) next
    gs_v <- abs(stats::cor(datExpr, bt$traits[[t]],
                           use = "pairwise.complete.obs"))[, 1L]
    mm_v <- kme[, me_name]
    rho <- suppressWarnings(stats::cor.test(gs_v, mm_v,
                                            )$estimate)
    dd <- data.frame(gs = gs_v, mm = mm_v, gene = colnames(datExpr),
                     stringsAsFactors = FALSE)
    hub_lab <- dd[dd$gs > quantile(dd$gs, 0.98) &
                  dd$mm > quantile(dd$mm, 0.98), , drop = FALSE]
    p_gsmm <- ggplot2::ggplot(dd, ggplot2::aes(x = mm, y = gs)) +
      ggplot2::geom_point(size = 0.8, alpha = 0.5, colour = PAL$primary) +
      ggplot2::labs(
        title = sprintf("GS vs MM - trait %s, module %s", t, best_m),
        subtitle = wrap_subtitle(sprintf(paste0(
          "GS = |cor(gene, trait)|, MM = kME. Spearman rho = %.2f. ",
          "Labelled: top 2% of both = hub candidates.",
          " High rho = module genuinely tracks this trait."),
          as.numeric(rho)), fig_width = W_ONE_HALF),
        x = sprintf("MM (kME, module %s)", best_m),
        y = sprintf("GS (|cor| with %s)", t)) +
      theme_paper(9)
    if (nrow(hub_lab) > 0L) {
      p_gsmm <- p_gsmm + ggrepel::geom_text_repel(
        data = hub_lab, ggplot2::aes(label = gene),
        size = 2.2, colour = PAL$ink,
        min.segment.length = 0)
    }
    save_pdf(file.path(res, paste0(GSMM_BASE, which(unique(cor_df$trait) == t),
                                   "-gs-mm-", t, ".pdf")),
             print(p_gsmm), width = W_ONE_HALF, height = mm(72))
  }  # 闭 for
  NULL
  }, error = function(e) sprintf("%s: %s", paste(class(e), collapse = "/"), conditionMessage(e)))
  if (!is.null(status$gsmm_error)) {
    status$gsmm_status <- "not_available"
    log_warn(sprintf("[WGCNA] GS-MM 散点未出（其余产物不受影响）: %s", status$gsmm_error))
  }

  write_json(file.path(res, "wgcna_status.json"), status)
  log_info(paste0("已生成 wgcna_modules.csv / wgcna_module_trait.csv / wgcna_module_sizes.csv / ",
                  "wgcna_soft_power.csv / wgcna_status.json"))
  invisible(mod_df)
}

if (!GEO_ORCHESTRATED()) {
  cfg <- load_config()
  run_06_wgcna(cfg)
}
