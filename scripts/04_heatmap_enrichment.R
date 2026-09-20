# ============================================================================
# 04_heatmap_enrichment.R — 聚类热图 + preranked GSEA + GO/KEGG 富集
# ============================================================================
# spec 的 deg_heatmap + go_kegg_enrich 步骤。
#
# 热图：top N 显著 DEG（按 adj.P.Val），不足时自动降级；行 Z-score，euclidean + complete。
#
# 富集分两条路，来自 K-Dense `pathway-enrichment` skill 的指引：
#   "a discrete hit list → ORA; a ranked table with per-gene scores → GSEA"
#   "Never threshold a list and then feed it to GSEA"
#
#   A. preranked GSEA（主力）—— 用**完整的 16,487 基因排序表**，不卡阈值，
#      排序指标是 limma 的 moderated t 统计量。弱功效、效应弥散的数据集
#      正是 GSEA 被设计出来处理的场景；卡阈值再跑 ORA 会扔掉排序信息。
#   B. ORA（辅助）—— 按上/下调**分开**跑。ORA 本身方向无关，混在一起
#      "富集到 X 通路"由上调还是下调基因驱动就分不清了。
#
# 两条路的结果都做**基因重叠去冗余**：GO 会返回几十个近义条目，
# 那是 1 个发现重复几十次，不是几十个发现。
#
# 富集为空或 KEGG 接口失败时写空表 + 状态文件，**不终止流程**。
#
# 输出：results/01-04-01-unit1-top50-heatmap.pdf
#       results/01-04-02-unit1-gsea-go-dotplot.pdf  / GSEA_GO_table.csv
#       results/01-04-03-unit1-gsea-kegg-dotplot.pdf / GSEA_KEGG_table.csv
#       results/01-04-04-unit1-go-ora-dotplot.pdf  / GO_table.csv      （含 direction 列）
#       results/01-04-05-unit1-kegg-ora-dotplot.pdf / KEGG_table.csv   （含 direction 列）
#       results/enrichment_status.json
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
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

#' 选出用于热图的基因：优先显著 DEG，不足时降级
#'
#' 返回 `label` 是**给读者看的图标题用词**，必须随实际取到的基因变：
#' 原来图标题写死 "Top DEG heatmap"，于是 GSE64790（0 个基因通过 FDR、
#' 回退到 top 20 by raw P）的图上仍然写着 "Top DEG" —— 而那张图里
#' 一个差异表达基因都没有。`mode` 是给日志和状态文件的技术描述，
#' `label` 是给图上的事实陈述，两者不能互相顶替。
select_heatmap_genes <- function(deg, cfg) {
  top_n <- cfg$analysis$top_heatmap_genes
  sig <- deg[deg$adj.P.Val < cfg$thresholds$adj_p &
             abs(deg$logFC) > cfg$thresholds$log2fc, , drop = FALSE]
  if (nrow(sig) >= top_n) {
    return(list(genes = head(sig$gene[order(sig$adj.P.Val)], top_n),
                mode = sprintf("top %d significant DEG by adj.P", top_n),
                label = sprintf("Top %d significant DEG", top_n)))
  }
  if (nrow(sig) > 0L) {
    return(list(genes = sig$gene[order(sig$adj.P.Val)],
                mode = sprintf("all %d significant DEG (< %d requested)", nrow(sig), top_n),
                label = sprintf("All %d significant DEG", nrow(sig))))
  }
  list(genes = head(deg$gene[order(deg$adj.P.Val)], 20),
       mode = "WARNING: no significant DEG, fell back to top 20 by raw P",
       label = "Top 20 genes by raw P - no gene passes FDR")
}

run_04a_heatmap <- function(cfg) {
  log_info("=== 步骤 04a：聚类热图 ===")
  ensure_dirs(cfg)

  res <- cfg$output$results_dir
  deg <- utils::read.csv(file.path(res, "deg_table.csv"), stringsAsFactors = FALSE)
  expr <- readRDS(file.path(cfg$output$data_dir, "expr_clean.rds"))
  group <- utils::read.csv(file.path(cfg$output$data_dir, "group.csv"), stringsAsFactors = FALSE)
  groups <- factor(group$group, levels = unique(group$group))

  # ---- 1. 聚类热图 --------------------------------------------------------
  pick <- select_heatmap_genes(deg, cfg)
  genes <- intersect(pick$genes, rownames(expr))
  if (length(genes) == 0L) stop("热图基因与表达矩阵无交集")
  log_info(sprintf("热图: %s（%d 个基因）", pick$mode, length(genes)))
  if (grepl("WARNING", pick$mode)) log_warn(pick$mode)

  mat <- row_zscore(expr[genes, , drop = FALSE])
  mat[is.na(mat)] <- 0
  # 行 Z-score 后截断到 ±3，避免极端值压平配色
  mat[mat > 3] <- 3
  mat[mat < -3] <- -3

  annotation_col <- data.frame(group = groups, row.names = colnames(expr))
  annotation_row <- data.frame(direction = ifelse(deg$logFC[match(genes, deg$gene)] > 0, "up", "down"),
                               row.names = genes)
  # 注释条用与火山图/PCA 完全相同的条件色与方向色
  annotation_colors <- list(
    group     = group_palette(groups, cfg),
    direction = c(up = PAL$up, down = PAL$down)
  )

  # **画布高度要长到装得下行名，而不是让行名被挤掉。**
  #
  # 原来每行固定 0.115in，而 `decide_rownames()` 的预算是
  #   budget = floor(height_in * 72 * panel_frac / (fontsize + min_gap))
  # 反解出每行至少需要 `(fontsize + min_gap) / (72 * panel_frac)` ≈ 0.127in。
  # 0.115 < 0.127，所以 50 个基因时预算只有 45 行 → 判定放不下 →
  # **整张不显示行名**。而热图一旦没有基因名，读者就无法知道画的是哪些基因
  # （规则 21 的另一半正是为此要求落盘 `top50_heatmap_genes.csv`）。
  #
  # 宽度受规则 25 约束（183mm），**高度没有上限** —— 所以正确做法是让画布长高，
  # 而不是牺牲标签。50 个基因 ≈ 169mm 高，仍在常规排版范围内（<200mm）。
  #
  # 每行间距**从 decide_rownames 的同一个算式反解**，不写字面量：
  # 两处各写一个数就会漂移，而漂移的后果是"行名又悄悄不显示了"。
  row_fs     <- 5
  min_gap    <- 2.5
  panel_frac <- 0.82   # pheatmap 无副标题、图例是右侧细色条，面板占比高于 ggplot
  # 1.05 的余量：per_row 按等式反解，等式上刚好等于预算，浮点误差会让
  # floor() 少算 1 行 —— 于是 50 个基因又变成"放不下"。余量把这个刀锋推开。
  per_row <- ((row_fs + min_gap) / (72 * panel_frac)) * 1.05
  fig_h   <- max(mm(140), length(genes) * per_row)
  # **行名放不下就整张不显示。** 原来写的是硬编码的 `length(genes) <= 60`，
  # 和画布高度毫无关系 —— 50 个基因时每行只剩约 3.6px 间隙，糊成一片。
  show_rn <- decide_rownames(length(genes), fig_h, row_fs, "热图基因",
                             panel_frac = panel_frac, min_gap = min_gap,
                             figure = "01-04-01-unit1-top50-heatmap")
  log_info(sprintf("热图画布 %.1f mm（%d 个基因，行名 %s）",
                   fig_h * 25.4, length(genes),
                   if (isTRUE(show_rn)) "显示" else "隐藏"))

  # **行聚类自己算，再把同一棵树传给 pheatmap。**
  # 一是为了拿到显示顺序（见下面的 CSV），二是保证表和图的顺序必然一致 ——
  # 若让 pheatmap 内部再算一次，两边就只是"应该一样"。列同理。
  hc_rows <- stats::hclust(stats::dist(mat, method = "euclidean"), method = "complete")
  hc_cols <- stats::hclust(stats::dist(t(mat), method = "euclidean"), method = "complete")

  # **列名放不放得下，同样算出来，不能靠默认。**
  # pheatmap 的列名默认旋转 90°，所以每个列名在**横向**占用的正是它的行高
  # `fontsize + min_gap` —— 与行名是同一个一维模型，只是把"可排布的那一维"
  # 从高度换成宽度，所以直接复用 decide_rownames()（它的 `height_in` 就是
  # "可排布那一维的长度"，这里传宽度）。
  # 原来这里**根本没传 show_colnames**，pheatmap 默认全画 —— 121 个 GSM 号
  # 在 183 mm 里每个只剩约 1.5 mm，糊成一条黑带，而"糊了"从图上完全看不出来。
  col_fs <- 5
  # 0.62 = 热图面板占整幅宽的比例：左侧树状图约 15 mm + direction 注释条约 4 mm，
  # 右侧色条与图例约 45 mm，183 - 64 = 119 mm，119/183 ≈ 0.65，取 0.62 留余量。
  col_panel_frac <- 0.62
  show_cn <- decide_rownames(ncol(mat), W_DOUBLE, col_fs, "热图样本名（列）",
                             panel_frac = col_panel_frac, min_gap = min_gap,
                             figure = "01-04-01-unit1-top50-heatmap")
  log_info(sprintf("热图列名 %s（%d 个样本，画布宽 %.1f mm）",
                   if (isTRUE(show_cn)) "显示" else "隐藏", ncol(mat), W_DOUBLE * 25.4))

  save_pdf(file.path(res, "01-04-01-unit1-top50-heatmap.pdf"), {
    pheatmap::pheatmap(
      mat,
      annotation_col = annotation_col,
      annotation_row = annotation_row,
      annotation_colors = annotation_colors,
      cluster_rows = hc_rows, cluster_cols = hc_cols,
      show_rownames = show_rn, fontsize_row = row_fs,
      show_colnames = show_cn, fontsize_col = col_fs,
      color = pal_diverging(100), border_color = "white",
      breaks = seq(-3, 3, length.out = 101),
      # **方向分布写进标题。**
      # 这张图按 adj.P 取前 N 个（规范如此，规则 9 只要求 ORA 分方向，没要求
      # 热图平衡），实测 GSE42568 取到的 50 个**全是下调** —— 而同一批数据的
      # 火山图里上调 1957 / 下调 1873，几乎对半。两张图并排看，读者会怀疑
      # 热图选错了基因。图本身没说谎（左侧 direction 注释条如实标了全蓝），
      # 但"为什么一个上调都没有"必须有交代，否则就是"内容表达有歧义"。
      # **这里只加说明，不动选择逻辑** —— 改选择等于改规范。
      main = sprintf("%s (row Z-score) - %s [%d up / %d down]", pick$label,
                     cfg$dataset_id,
                     sum(annotation_row$direction == "up"),
                     sum(annotation_row$direction == "down")),
      silent = FALSE
    )
    # pheatmap 的色条固定在图右侧、无位置参数；细长条，占宽有限，保留。
  }, width = W_DOUBLE, height = fig_h)
  log_info("已生成 01-04-01-unit1-top50-heatmap.pdf")

  # **行名一旦不显示，基因身份就只剩这张表能提供。**
  # 而且必须是**显示顺序** —— 热图按聚类重排行，写 `genes` 的原始顺序
  # 会对不上图上的第 N 行。规则 16（布局要落盘）在热图上同样适用：
  # 从 PNG 反推"第 12 行是哪个基因"是猜，落盘之后就是可核对的数据。
  row_order <- hc_rows$order
  disp <- data.frame(
    display_row = seq_along(row_order),
    gene        = genes[row_order],
    direction   = annotation_row$direction[row_order],
    logFC       = deg$logFC[match(genes[row_order], deg$gene)],
    adj_P_Val   = deg$adj.P.Val[match(genes[row_order], deg$gene)],
    stringsAsFactors = FALSE
  )
  utils::write.csv(disp, file.path(res, "top50_heatmap_genes.csv"), row.names = FALSE)
  log_info(sprintf("已生成 top50_heatmap_genes.csv（%d 行，按图上的显示顺序；行名%s）",
                   nrow(disp), if (show_rn) "已显示" else "未显示，靠这张表对照"))

  # **列名同理。** 藏了列名之后，"图上第 N 列是哪个样本"就只剩这张表能回答，
  # 而且必须是**显示顺序**（列也按聚类重排）。顺序取自上面自己算的那棵树，
  # 与图上必然一致 —— 让 pheatmap 内部再算一次的话，两边只是"应该一样"。
  col_order <- hc_cols$order
  samp <- data.frame(
    display_col = seq_along(col_order),
    sample      = colnames(mat)[col_order],
    group       = as.character(groups[col_order]),
    stringsAsFactors = FALSE
  )
  utils::write.csv(samp, file.path(res, "top50_heatmap_samples.csv"), row.names = FALSE)
  log_info(sprintf("已生成 top50_heatmap_samples.csv（%d 列，按图上的显示顺序；列名%s）",
                   nrow(samp), if (show_cn) "已显示" else "未显示，靠这张表对照"))

  invisible(list(genes = genes, mode = pick$mode,
                 show_rownames = show_rn, show_colnames = show_cn))
}

run_04b_enrichment <- function(cfg) {
  log_info("=== 步骤 04b：preranked GSEA + GO / KEGG 富集 ===")
  ensure_dirs(cfg)

  res <- cfg$output$results_dir
  deg <- utils::read.csv(file.path(res, "deg_table.csv"), stringsAsFactors = FALSE)
  expr <- readRDS(file.path(cfg$output$data_dir, "expr_clean.rds"))

  status <- list(go = list(status = "not_run"), kegg = list(status = "not_run"),
                 gsea_go = list(status = "not_run"), gsea_kegg = list(status = "not_run"))

  # 富集需要基因 symbol；01 若只能拿到探针 ID，这里必须跳过而不是硬凑
  feat_path <- file.path(cfg$output$data_dir, "feature_mode.json")
  feat <- if (file.exists(feat_path)) {
    tryCatch(jsonlite::fromJSON(feat_path, simplifyVector = FALSE), error = function(e) NULL)
  } else NULL
  if (!is.null(feat)) status$feature_mode <- feat$mode
  if (!is.null(feat) && !identical(feat$mode, "symbol")) {
    reason <- sprintf("特征无法映射到基因 symbol（%s）；GO/KEGG 需要基因 symbol，已跳过",
                      feat$reason %||% "原因未知")
    log_warn(reason)
    write_empty_enrichment(cfg, status, reason)
    return(invisible(NULL))
  }

  sel <- select_degs(deg, cfg, min_genes = 5L)
  sig_genes <- sel$genes
  status$deg_mode <- sel$mode
  status$deg_reason <- sel$reason
  status$input_genes <- length(sig_genes)
  status$universe <- cfg$enrichment$universe
  if (identical(sel$mode, "ranked_fallback")) {
    log_warn(sprintf("无基因通过 FDR，ORA 改用 raw P 排序前 %d 个基因（假设生成，非显著 DEG）",
                     length(sig_genes)))
  }

  # 背景集：genome = OrgDb 全部基因；detected = 芯片实测基因。
  # 默认 detected —— 背景应当是"本实验可能检出的基因"，用全基因组会让管家类
  # 条目假显著（K-Dense pathway-enrichment 点名的 ORA 头号误导来源）。
  universe_symbols <- NULL
  if (identical(cfg$enrichment$universe, "detected")) {
    universe_symbols <- rownames(expr)
    log_info(sprintf("ORA 背景: 实测基因集（%d 个）", length(universe_symbols)))
  } else {
    log_info("ORA 背景: 全基因组（OrgDb 默认）—— 注意这会让泛化条目显得更显著")
  }

  jac <- cfg$enrichment$redundancy_jaccard
  if (is.null(jac)) jac <- 0.5

  # ==========================================================================
  # A. preranked GSEA —— 弱功效数据集的主力方法
  # ==========================================================================
  #
  # **为什么 GSEA 是主力而不是 ORA：** K-Dense `pathway-enrichment`：
  #   "a discrete hit list → ORA; a ranked table with per-gene scores → GSEA"
  #   "Never threshold a list and then feed it to GSEA"
  #   "Better when effects are broad/subtle or when a hit list would be very
  #    short or very long"（> 2000 个基因的 ORA 列表会失去特异性）
  #
  # 本设计手上是完整的 16,487 基因排序表，效应弥散、无一通过 FDR ——
  # 正是 GSEA 被设计出来处理的场景。用 ORA 需要先卡阈值，
  # 而卡阈值恰好扔掉了 GSEA 所依赖的排序信息。
  #
  # 排序指标用 limma 的 moderated t 统计量，不是 log2FC：
  #   "Rank by the test statistic (sign = direction, magnitude = evidence).
  #    This is more stable than ranking by log2FoldChange, which is noisy for
  #    low-count genes."
  gsea_ranks <- NULL
  {
    gl <- deg$t
    names(gl) <- deg$gene
    gl <- gl[!is.na(gl) & is.finite(gl) & nzchar(names(gl))]
    gl <- gl[!duplicated(names(gl))]
    gl <- sort(gl, decreasing = TRUE)
    if (length(gl) >= 100L) {
      map <- tryCatch(
        clusterProfiler::bitr(names(gl), fromType = "SYMBOL", toType = "ENTREZID",
                              OrgDb = org.Hs.eg.db),
        error = function(e) NULL)
      if (!is.null(map) && nrow(map) > 0L) {
        gl_e <- gl[map$SYMBOL]
        names(gl_e) <- map$ENTREZID
        gl_e <- gl_e[!duplicated(names(gl_e))]
        gsea_ranks <- sort(gl_e, decreasing = TRUE)
        log_info(sprintf("GSEA 排序表: %d 个基因（指标 = limma moderated t，未卡阈值）",
                         length(gsea_ranks)))
      }
    }
  }

  # fgsea 的 p 值是蒙特卡洛估计的（fgseaMultilevel 自适应采样），因此必须固定 RNG。
  #
  # **注意：`gseGO(seed=123)` 这个参数不足以复现。** 实测连续两轮 CI 的
  # GSEA 显著条目数是 1072 和 1095 —— 参数被接受了（日志里没有退回警告），
  # 但结果仍然每次都不同，说明它没有被真正转发到 fgsea 的采样器。
  # 所以在调用**之前**直接设 RNG 状态，不依赖 clusterProfiler 的转发。
  gsea_call <- function(f, ...) {
    seed <- cfg$analysis$seed
    if (is.null(seed)) seed <- 123
    # 两处都设：set.seed 管 R 全局 RNG；显式 seed= 管那些自己开 RNG 流的实现
    set.seed(seed)
    tryCatch(f(seed = seed, ...), error = function(e) {
      log_warn(sprintf("GSEA 传 seed 失败，退回默认置换: %s", conditionMessage(e)))
      set.seed(seed)
      f(...)
    })
  }

  gsea_go_res <- NULL
  if (!is.null(gsea_ranks)) {
    gsea_go_res <- tryCatch(
      gsea_call(clusterProfiler::gseGO, geneList = gsea_ranks, OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID", ont = cfg$enrichment$ont,
                minGSSize = cfg$analysis$gsea_min_set,
                maxGSSize = cfg$analysis$gsea_max_set,
                pvalueCutoff = cfg$enrichment$pvalue_cutoff,
                pAdjustMethod = cfg$enrichment$p_adjust, verbose = FALSE),
      error = function(e) {
        log_warn(sprintf("GSEA (GO) 失败: %s", conditionMessage(e)))
        status$gsea_go <<- list(status = "failed", reason = conditionMessage(e))
        NULL
      })
    if (!is.null(gsea_go_res)) {
      df <- as.data.frame(gsea_go_res)
      if (nrow(df) > 0L) {
        df <- reduce_terms_by_overlap(df, jac, gene_col = "core_enrichment")
        df <- normalise_gene_lists(df)
        utils::write.csv(df, file.path(res, "GSEA_GO_table.csv"), row.names = FALSE)
        n_rep <- length(unique(stats::na.omit(df$representative)))
        log_info(sprintf("GSEA (GO %s): %d 条显著，去冗余后 %d 个代表条目",
                         cfg$enrichment$ont, nrow(df), n_rep))
        status$gsea_go <- list(
          status = "ok", terms = nrow(df), representative_terms = n_rep,
          top = head(df$Description[order(df$p.adjust)], 5),
          top_up = head(df$Description[df$NES > 0][order(df$p.adjust[df$NES > 0])], 3),
          top_down = head(df$Description[df$NES < 0][order(df$p.adjust[df$NES < 0])], 3))
        decide_rownames(min(2 * cfg$enrichment$top_terms, nrow(df)), 6.5, 7,
                        "GSEA GO 点图", panel_frac = 0.75, min_gap = 2.5,
                        figure = "01-04-02-unit1-gsea-go-dotplot")
        save_pdf(file.path(res, "01-04-02-unit1-gsea-go-dotplot.pdf"),
                 print(make_gsea_dotplot(df, cfg,
                         sprintf("GSEA (preranked) GO %s - %s", cfg$enrichment$ont,
                                 cfg$dataset_id))),
                 width = W_DOUBLE, height = mm(165))
      } else {
        status$gsea_go <- list(status = "empty", reason = "no gene set passed the cutoff")
        utils::write.csv(data.frame(), file.path(res, "GSEA_GO_table.csv"), row.names = FALSE)
      }
    }
  } else {
    status$gsea_go <- list(status = "skipped", reason = "排序表不足 100 个基因，GSEA 无意义")
  }

  gsea_kegg_res <- NULL
  if (!is.null(gsea_ranks)) {
    gsea_kegg_res <- tryCatch(
      gsea_call(clusterProfiler::gseKEGG, geneList = gsea_ranks,
                organism = cfg$enrichment$kegg_organism, keyType = "kegg",
                minGSSize = cfg$analysis$gsea_min_set,
                maxGSSize = cfg$analysis$gsea_max_set,
                pvalueCutoff = cfg$enrichment$pvalue_cutoff,
                pAdjustMethod = cfg$enrichment$p_adjust, verbose = FALSE),
      error = function(e) {
        log_warn(sprintf("GSEA (KEGG) 失败（不终止流程）: %s", conditionMessage(e)))
        status$gsea_kegg <<- list(status = "failed", reason = conditionMessage(e))
        NULL
      })
    if (!is.null(gsea_kegg_res)) {
      df <- as.data.frame(gsea_kegg_res)
      if (nrow(df) > 0L) {
        df <- reduce_terms_by_overlap(df, jac, gene_col = "core_enrichment")
        df <- normalise_gene_lists(df)
        utils::write.csv(df, file.path(res, "GSEA_KEGG_table.csv"), row.names = FALSE)
        n_rep <- length(unique(stats::na.omit(df$representative)))
        log_info(sprintf("GSEA (KEGG): %d 条显著，去冗余后 %d 个代表条目", nrow(df), n_rep))
        status$gsea_kegg <- list(
          status = "ok", terms = nrow(df), representative_terms = n_rep,
          top = head(df$Description[order(df$p.adjust)], 5))
        decide_rownames(min(2 * cfg$enrichment$top_terms, nrow(df)), 6.5, 7,
                        "GSEA KEGG 点图", panel_frac = 0.75, min_gap = 2.5,
                        figure = "01-04-03-unit1-gsea-kegg-dotplot")
        save_pdf(file.path(res, "01-04-03-unit1-gsea-kegg-dotplot.pdf"),
                 print(make_gsea_dotplot(df, cfg,
                         sprintf("GSEA (preranked) KEGG - %s", cfg$dataset_id))),
                 width = W_DOUBLE, height = mm(165))
      } else {
        status$gsea_kegg <- list(status = "empty", reason = "no pathway passed the cutoff")
        utils::write.csv(data.frame(), file.path(res, "GSEA_KEGG_table.csv"), row.names = FALSE)
      }
    }
  } else {
    status$gsea_kegg <- list(status = "skipped", reason = "排序表不足 100 个基因，GSEA 无意义")
  }

  # ==========================================================================
  # B. ORA —— 按上/下调**分开**跑
  # ==========================================================================
  #
  # ORA 本身方向无关：混在一起跑，"富集到 X 通路"到底由上调还是下调基因驱动
  # 就分不清了（K-Dense pathway-enrichment："ORA is direction-agnostic unless
  # you split up/down lists"）。肿瘤 vs 正常组织尤其致命 ——
  # 上调的是增殖、下调的是基质/脂肪，混起来会出现"方向相反的通路同时富集"。
  if (length(sig_genes) < 5L) {
    log_warn(sprintf("可用于 ORA 的基因仅 %d 个（< 5），跳过 GO/KEGG", length(sig_genes)))
    write_empty_enrichment(cfg, status,
                           sprintf("只有 %d 个基因可用于富集，不足 5 个", length(sig_genes)))
    return(invisible(NULL))
  }

  # 名字刻意不叫 run_* —— 那个前缀是留给编排器调用的步骤函数的
  # （tools/check_r_syntax.mjs 会检查 run_* 是否被 main_analysis.R 引用）
  ora_for_direction <- function(genes, direction) {
    if (length(genes) < 5L) {
      log_info(sprintf("ORA (%s): 只有 %d 个基因，跳过", direction, length(genes)))
      return(NULL)
    }
    r <- tryCatch({
      go <- clusterProfiler::enrichGO(
        gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
        ont = cfg$enrichment$ont, universe = universe_symbols,
        pAdjustMethod = cfg$enrichment$p_adjust,
        pvalueCutoff = cfg$enrichment$pvalue_cutoff,
        qvalueCutoff = cfg$enrichment$qvalue_cutoff, readable = TRUE)
      go_df <- as.data.frame(go)
      if (nrow(go_df) == 0L) return(NULL)
      go_df$direction <- direction
      go_df$n_input <- length(genes)
      go_df
    }, error = function(e) {
      log_warn(sprintf("ORA GO (%s) 失败: %s", direction, conditionMessage(e)))
      NULL
    })

    k <- tryCatch({
      ez <- clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                                  OrgDb = org.Hs.eg.db)$ENTREZID
      if (length(ez) < 5L) return(NULL)
      ue <- NULL
      if (!is.null(universe_symbols)) {
        ue <- clusterProfiler::bitr(universe_symbols, fromType = "SYMBOL",
                                    toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
      }
      kk <- clusterProfiler::enrichKEGG(
        gene = ez, organism = cfg$enrichment$kegg_organism, keyType = "kegg",
        universe = ue, pAdjustMethod = cfg$enrichment$p_adjust,
        pvalueCutoff = cfg$enrichment$pvalue_cutoff,
        qvalueCutoff = cfg$enrichment$qvalue_cutoff)
      kk_df <- as.data.frame(kk)
      if (nrow(kk_df) == 0L) return(NULL)
      kk_df$direction <- direction
      kk_df$n_input <- length(genes)
      kk_df
    }, error = function(e) {
      log_warn(sprintf("ORA KEGG (%s) 失败（不终止流程）: %s", direction, conditionMessage(e)))
      NULL
    })
    list(go = r, kegg = k)
  }

  ora_up   <- ora_for_direction(sel$up,   "up")
  ora_down <- ora_for_direction(sel$down, "down")
  log_info(sprintf("ORA 输入: 上调 %d 个 / 下调 %d 个基因（方向已分开）",
                   length(sel$up), length(sel$down)))

  combine_ora <- function(pick) {
    parts <- Filter(Negate(is.null), list(pick(ora_up), pick(ora_down)))
    if (length(parts) == 0L) return(NULL)
    cols <- Reduce(intersect, lapply(parts, colnames))
    do.call(rbind, lapply(parts, function(p) p[, cols, drop = FALSE]))
  }

  go_df   <- combine_ora(function(x) if (is.null(x)) NULL else x$go)
  kegg_df <- combine_ora(function(x) if (is.null(x)) NULL else x$kegg)

  # 去冗余 + 落盘
  #
  # `fig_name` 与 `file_base` **分开两个参数**：前者是出图名（按仓库约定
  # 带阶段-模块-图-单元前缀），后者是 CSV 表名前缀（`GO_table.csv`）。
  # 一个参数兼两用会让改图名顺带改掉数据产物的名字 —— 而数据产物
  # 是别的脚本和验收项在引用的。
  emit_ora <- function(df, label, fig_name, file_base, key) {
    if (is.null(df) || nrow(df) == 0L) {
      utils::write.csv(data.frame(), file.path(res, paste0(file_base, "_table.csv")),
                       row.names = FALSE)
      status[[key]] <<- list(status = "empty", reason = "no term passed the cutoff")
      return(invisible(NULL))
    }
    # 去冗余要**按方向分别做** —— 上调与下调的条目本来就不该互相折叠
    out <- do.call(rbind, lapply(split(df, df$direction), function(d) {
      reduce_terms_by_overlap(d, jac, gene_col = "geneID")
    }))
    rownames(out) <- NULL
    out <- normalise_gene_lists(out)
    utils::write.csv(out, file.path(res, paste0(file_base, "_table.csv")), row.names = FALSE)
    n_rep <- length(unique(stats::na.omit(out$representative)))
    log_info(sprintf("%s: %d 条（上/下调分开），去冗余后 %d 个代表条目", label, nrow(out), n_rep))
    status[[key]] <<- list(
      status = "ok", terms = nrow(out), representative_terms = n_rep,
      direction_split = TRUE,
      top_up = head(out$Description[out$direction == "up"][order(out$p.adjust[out$direction == "up"])], 3),
      top_down = head(out$Description[out$direction == "down"][order(out$p.adjust[out$direction == "down"])], 3))
    # 两个面板并排，每个面板 15 行标签 —— 比原来单面板 30 行宽松一倍。
    # 标签仍按量化判据核一遍：放不下就整张不显示，不缩字号硬塞。
    ora_w <- W_DOUBLE; ora_h <- mm(165)   # 10 in = 254 mm，装不进一页
    decide_rownames(min(cfg$enrichment$top_terms, nrow(out)), ora_h, 7,
                    sprintf("%s 点图", label), panel_frac = 0.68, min_gap = 2.5,
                    figure = fig_name)
    save_pdf(file.path(res, paste0(fig_name, ".pdf")),
             print(make_ora_dotplot(out, cfg, sprintf("%s - %s", label, cfg$dataset_id))),
             width = ora_w, height = ora_h)
    invisible(NULL)
  }

  emit_ora(go_df, sprintf("GO %s ORA (up/down split)", cfg$enrichment$ont),
           "01-04-04-unit1-go-ora-dotplot", "GO", "go")
  emit_ora(kegg_df, "KEGG ORA (up/down split)",
           "01-04-05-unit1-kegg-ora-dotplot", "KEGG", "kegg")

  # 排序指标与置换设置要记录 —— 可复现性清单要求
  status$gsea_ranking_metric <- "limma moderated t statistic (sign = direction)"
  status$gsea_engine <- "clusterProfiler::gseGO / gseKEGG (fgsea)"
  status$gsea_seed <- cfg$analysis$seed
  status$gsea_set_size <- c(cfg$analysis$gsea_min_set, cfg$analysis$gsea_max_set)
  status$redundancy_jaccard <- jac
  write_json(file.path(res, "enrichment_status.json"), status)
  log_info("已生成 enrichment_status.json")

  # 标签决策落盘。**放这里而不是随每张图写** —— 一次性写出本轮所有决策，
  # 便于横向比较"哪张图被藏了行名、为什么"。日志里虽然有同样的算式，
  # 但 CI 日志会滚掉，文件不会。
  write_label_decisions(file.path(res, "label_decisions.csv"))
  invisible(NULL)
}

#' ORA 的 dotplot：**按方向分面**
#'
#' 原来 x 轴是方向，两个列标签都写成 `"<arm>\n(up in <arm>)"` —— 于是
#' **整张图上 "down" 这个词一次都不出现**，读者会以为只有上调的富集。
#' 实测就是这样被问的（"怎么只有上调的富集没有下调的"）。下调的点其实画了
#' （连通域数得出来两列都有点），是标注把人骗了。
#'
#' 改成按方向分面：
#'   * 两个面板各有标题与条目数，方向不可能看漏；
#'   * y 轴标签从 30 行降到每面板 15 行，密集问题一并缓解；
#'   * x 轴腾出来放显著性（原来被方向占着，-log10 P 只能塞进颜色）。
make_ora_dotplot <- function(df, cfg, title) {
  n <- cfg$enrichment$top_terms
  arms <- as.character(cfg$contrast)
  dir_name <- c(up = sprintf("up in %s", arms[1L]),
                down = sprintf("down in %s", arms[1L]))

  keep <- do.call(rbind, lapply(c("up", "down"), function(d) {
    sub <- df[df$direction == d, , drop = FALSE]
    if (nrow(sub) == 0L) return(NULL)
    utils::head(sub[order(sub$p.adjust), , drop = FALSE], n)
  }))
  if (is.null(keep) || nrow(keep) == 0L) return(NULL)
  rownames(keep) <- NULL

  # 面板标题带条目数：读者一眼看出两边各有多少条，不会怀疑某边是空的
  counts <- as.integer(table(factor(keep$direction, levels = c("up", "down"))))
  panel_lab <- sprintf("%s\n(%d terms shown)", dir_name[c("up", "down")], counts)
  keep$panel <- factor(panel_lab[match(keep$direction, c("up", "down"))],
                       levels = panel_lab)
  keep$Description <- factor(keep$Description,
                             levels = unique(keep$Description[order(keep$p.adjust,
                                                                    decreasing = TRUE)]))
  ggplot2::ggplot(keep, ggplot2::aes(x = -log10(p.adjust), y = Description)) +
    ggplot2::geom_point(ggplot2::aes(size = Count, colour = direction)) +
    # 颜色仍走方向色（与火山图/PCA/热图注释条同源）；图例关掉，
    # 因为分面标题已经把方向写在脸上了，再放一个图例是重复。
    ggplot2::scale_colour_manual(values = c(up = PAL$up, down = PAL$down),
                                 guide = "none") +
    # **气泡要收一点**（评审 3.4：行距 ≈25px 而气泡 Ø22–26px，相邻相交）。
    # range 上限从 7 降到 5.5，行距不变时相邻气泡不再相切。
    ggplot2::scale_size_continuous(name = "genes", range = c(1.8, 5.5)) +
    # **右侧留余量**（评审 3.3：GSE42568 KEGG 右侧墨迹距边框仅 2px，
    # 最大的点被边框切平）。默认 expand 把最大点顶到面板边上。
    # 注意参数名是 expand（不是 expansion —— expansion() 只是造取值的函数）。
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.10))) +
    ggplot2::facet_wrap(~ panel, scales = "free_y", nrow = 1) +
    ggplot2::labs(title = title,
                  subtitle = wrap_subtitle(sprintf(
                    paste0("UP TO %d per direction - each panel header gives the count ",
                           "actually shown, which is lower when a direction has fewer ",
                           "terms passing the cutoff. ORA itself is direction-agnostic, ",
                           "so up and down are run as separate gene lists. ",
                           "x = significance, size = number of genes in the term."), n),
                    fig_width = W_DOUBLE),
                  x = expression(-log[10] ~ "(adj.P)"), y = NULL) +
    theme_paper(9) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 7),
      strip.text  = ggplot2::element_text(size = 9, face = "bold"),
      panel.spacing = ggplot2::unit(1.2, "lines"))
}

#' GSEA 的 dotplot：x 轴是 NES，方向直接由符号给出
make_gsea_dotplot <- function(df, cfg, title) {
  keep <- utils::head(df[order(df$p.adjust), , drop = FALSE], 2 * cfg$enrichment$top_terms)
  keep$Description <- factor(keep$Description,
                             levels = keep$Description[order(keep$NES)])
  arms <- as.character(cfg$contrast)
  ggplot2::ggplot(keep, ggplot2::aes(x = NES, y = Description)) +
    # **不加红/蓝背景分区。** 加了会在同一张图里出现两种"红"：
    # 背景红表示方向，点的红表示显著性，读者无法区分。
    # 方向由 NES 的符号承载，轴标签已经写明哪边是哪个条件。
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3,
                        colour = PAL$ink) +
    ggplot2::geom_point(ggplot2::aes(size = setSize, colour = -log10(p.adjust))) +
    # **色标范围必须可见、且不能被数据夹死成一条缝**（评审 3.5）：
    # GSE42568 的 GSEA GO 全部条目挤在 -log10 adj.P 7.65–7.75（0.10 宽），
    # gradientn 默认铺满整个渐变 → 几乎全黄、单个青点像离群。
    # 修法：limits 固定 0 到本图最大值的 1.15 倍（而不是 min–max），
    # 点群落在色带的中后段，"这批条目的显著性其实都在同一量级"这件事
    # 从色标上一眼可见；GSE64790（ranked_fallback，量级低）也不受影响。
    scale_colour_seq("-log10\nadj.P",
                     limits = c(0, 1.15 * max(-log10(keep$p.adjust))),
                     breaks = function(x) pretty(x, n = 4)) +
    ggplot2::scale_size_continuous(name = "set size", range = c(1.8, 5.5)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.04, 0.06))) +
    ggplot2::labs(title = title,
                  # 原来写的是 "right = up in tumor, left = up in normal" ——
                  # 两个方向又都写成 "up"，和 ORA 那张图是同一个毛病。
                  subtitle = wrap_subtitle(sprintf(
                    "preranked on the full gene list (no threshold); right = up in %s (= down in %s), left = down in %s",
                    arms[1L], arms[2L], arms[1L]), fig_width = W_DOUBLE),
                  x = "NES (normalized enrichment score)", y = NULL) +
    theme_paper(9) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))
}

#' 无显著基因时统一写空结果
#'
#' 注意 `status$go` / `status$kegg` 只是**覆盖** ORA 两项，
#' `deg_mode` / `deg_reason` / `universe` 等已经填好的字段要原样保留 ——
#' 验收和结论都要引用它们。
write_empty_enrichment <- function(cfg, status, reason) {
  res <- cfg$output$results_dir
  for (f in c("GO_table.csv", "KEGG_table.csv", "GSEA_GO_table.csv", "GSEA_KEGG_table.csv")) {
    utils::write.csv(data.frame(), file.path(res, f), row.names = FALSE)
  }
  status$go <- list(status = "skipped", reason = reason)
  status$kegg <- list(status = "skipped", reason = reason)
  if (is.null(status$gsea_go))  status$gsea_go  <- list(status = "skipped", reason = reason)
  if (is.null(status$gsea_kegg)) status$gsea_kegg <- list(status = "skipped", reason = reason)
  write_json(file.path(res, "enrichment_status.json"), status)
}

#' 通用 dotplot（enrichplot），保留给需要它的调用方
make_dotplot <- function(x, cfg, title) {
  n <- min(cfg$enrichment$top_terms, nrow(as.data.frame(x)))
  p <- enrichplot::dotplot(x, showCategory = n) +
    labs(title = title, subtitle = wrap_subtitle(sprintf(
      "top %d terms, p.adjust < %g", n, cfg$enrichment$pvalue_cutoff), fig_width = W_DOUBLE)) +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_text(size = 7),
          legend.position = "bottom", legend.direction = "horizontal")
  p
}

if (!GEO_ORCHESTRATED()) {
  cfg <- load_config()
  run_04a_heatmap(cfg)
  run_04b_enrichment(cfg)
}
