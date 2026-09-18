# ============================================================================
# 04_heatmap_enrichment.R — 聚类热图 + GO/KEGG 富集
# ============================================================================
# spec 的 deg_heatmap + go_kegg_enrich 步骤。
#
# 热图：top N 显著 DEG（按 adj.P.Val），不足时自动降级；行 Z-score，euclidean + complete。
# 富集：GO BP（enrichGO）+ KEGG（enrichKEGG）。
#       富集为空或 KEGG 接口失败时写空表 + 状态文件，**不终止流程**。
#
# 输出：results/top50_heatmap.pdf
#       results/GO_dotplot.pdf / GO_table.csv
#       results/KEGG_dotplot.pdf / KEGG_table.csv
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
select_heatmap_genes <- function(deg, cfg) {
  top_n <- cfg$analysis$top_heatmap_genes
  sig <- deg[deg$adj.P.Val < cfg$thresholds$adj_p &
             abs(deg$logFC) > cfg$thresholds$log2fc, , drop = FALSE]
  if (nrow(sig) >= top_n) {
    return(list(genes = head(sig$gene[order(sig$adj.P.Val)], top_n),
                mode = sprintf("top %d significant DEG by adj.P", top_n)))
  }
  if (nrow(sig) > 0L) {
    return(list(genes = sig$gene[order(sig$adj.P.Val)],
                mode = sprintf("all %d significant DEG (< %d requested)", nrow(sig), top_n)))
  }
  list(genes = head(deg$gene[order(deg$adj.P.Val)], 20),
       mode = "WARNING: no significant DEG, fell back to top 20 by raw P")
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

  save_pdf(file.path(res, "top50_heatmap.pdf"), {
    pheatmap::pheatmap(
      mat,
      annotation_col = annotation_col,
      annotation_row = annotation_row,
      cluster_rows = TRUE, cluster_cols = TRUE,
      clustering_distance_rows = "euclidean", clustering_method = "complete",
      clustering_distance_cols = "euclidean",
      show_rownames = length(genes) <= 60, fontsize_row = 5,
      color = grDevices::colorRampPalette(c("#2E5FA3", "white", "#C1443C"))(100),
      breaks = seq(-3, 3, length.out = 101),
      main = sprintf("Top DEG heatmap (row Z-score) - %s", cfg$dataset_id),
      silent = FALSE
    )
  }, width = 8, height = max(6, length(genes) * 0.13))
  log_info("已生成 top50_heatmap.pdf")
  invisible(list(genes = genes, mode = pick$mode))
}

run_04b_enrichment <- function(cfg) {
  log_info("=== 步骤 04b：GO / KEGG 富集 ===")
  ensure_dirs(cfg)

  res <- cfg$output$results_dir
  deg <- utils::read.csv(file.path(res, "deg_table.csv"), stringsAsFactors = FALSE)
  expr <- readRDS(file.path(cfg$output$data_dir, "expr_clean.rds"))

  # ---- 富集分析 -----------------------------------------------------------
  status <- list(go = list(status = "not_run"), kegg = list(status = "not_run"))

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
  sig_genes <- deg$gene[deg$adj.P.Val < cfg$thresholds$adj_p &
                        abs(deg$logFC) > cfg$thresholds$log2fc]

  if (length(sig_genes) == 0L) {
    log_warn("没有显著 DEG，跳过 GO/KEGG 富集（写空表）")
    write_empty_enrichment(cfg, status, "no significant DEG to test")
    return(invisible(NULL))
  }
  log_info(sprintf("富集输入基因数: %d", length(sig_genes)))

  # 背景集：genome = OrgDb 全部基因；detected = 芯片实测基因
  universe_symbols <- NULL
  if (identical(cfg$enrichment$universe, "detected")) {
    universe_symbols <- rownames(expr)
    log_info(sprintf("富集背景: 实测基因集（%d 个）", length(universe_symbols)))
  } else {
    log_info("富集背景: 全基因组（OrgDb 默认）")
  }

  # ---- 2a. GO -------------------------------------------------------------
  go_res <- tryCatch({
    clusterProfiler::enrichGO(
      gene = sig_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
      ont = cfg$enrichment$ont,
      universe = universe_symbols,
      pAdjustMethod = cfg$enrichment$p_adjust,
      pvalueCutoff = cfg$enrichment$pvalue_cutoff,
      qvalueCutoff = cfg$enrichment$qvalue_cutoff,
      readable = TRUE
    )
  }, error = function(e) {
    log_warn(sprintf("GO 富集失败: %s", conditionMessage(e)))
    status$go <<- list(status = "failed", reason = conditionMessage(e))
    NULL
  })

  if (!is.null(go_res)) {
    go_df <- as.data.frame(go_res)
    utils::write.csv(go_df, file.path(res, "GO_table.csv"), row.names = FALSE)
    if (nrow(go_df) > 0L) {
      save_pdf(file.path(res, "GO_dotplot.pdf"),
               print(make_dotplot(go_res, cfg, sprintf("GO %s enrichment - %s",
                                                       cfg$enrichment$ont, cfg$dataset_id))),
               width = 9, height = 7)
      log_info(sprintf("GO %s 富集: %d 条通路，已生成 GO_dotplot.pdf",
                       cfg$enrichment$ont, nrow(go_df)))
      status$go <- list(status = "ok", terms = nrow(go_df),
                        top = head(go_df$Description, 5))
    } else {
      log_warn("GO 富集结果为空，写空表")
      status$go <- list(status = "empty", reason = "no term passed the cutoff")
    }
  }

  # ---- 2b. KEGG -----------------------------------------------------------
  # KEGG 走在线 REST API，可能因网络/限流/授权失败 —— 按 spec 不终止流程
  kegg_res <- tryCatch({
    entrez <- clusterProfiler::bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID",
                                    OrgDb = org.Hs.eg.db)$ENTREZID
    universe_entrez <- NULL
    if (!is.null(universe_symbols)) {
      universe_entrez <- clusterProfiler::bitr(universe_symbols, fromType = "SYMBOL",
                                               toType = "ENTREZID",
                                               OrgDb = org.Hs.eg.db)$ENTREZID
    }
    log_info(sprintf("KEGG 输入: %d 个基因映射到 ENTREZ", length(entrez)))
    r <- clusterProfiler::enrichKEGG(
      gene = entrez, organism = cfg$enrichment$kegg_organism, keyType = "kegg",
      universe = universe_entrez,
      pAdjustMethod = cfg$enrichment$p_adjust,
      pvalueCutoff = cfg$enrichment$pvalue_cutoff,
      qvalueCutoff = cfg$enrichment$qvalue_cutoff
    )
    if (!is.null(r) && nrow(as.data.frame(r)) > 0L) {
      r <- clusterProfiler::setReadable(r, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    }
    r
  }, error = function(e) {
    log_warn(sprintf("KEGG 富集失败（不终止流程）: %s", conditionMessage(e)))
    status$kegg <<- list(status = "failed", reason = conditionMessage(e))
    NULL
  })

  if (!is.null(kegg_res)) {
    kegg_df <- as.data.frame(kegg_res)
    utils::write.csv(kegg_df, file.path(res, "KEGG_table.csv"), row.names = FALSE)
    if (nrow(kegg_df) > 0L) {
      save_pdf(file.path(res, "KEGG_dotplot.pdf"),
               print(make_dotplot(kegg_res, cfg, sprintf("KEGG pathway enrichment - %s",
                                                         cfg$dataset_id))),
               width = 9, height = 7)
      log_info(sprintf("KEGG 富集: %d 条通路，已生成 KEGG_dotplot.pdf", nrow(kegg_df)))
      status$kegg <- list(status = "ok", terms = nrow(kegg_df), top = head(kegg_df$Description, 5))
    } else {
      log_warn("KEGG 富集结果为空，写空表")
      status$kegg <- list(status = "empty", reason = "no pathway passed the cutoff")
    }
  }

  if (is.null(go_res) || nrow(as.data.frame(go_res)) == 0L) {
    utils::write.csv(data.frame(), file.path(res, "GO_table.csv"), row.names = FALSE)
  }
  if (is.null(kegg_res) || nrow(as.data.frame(kegg_res)) == 0L) {
    utils::write.csv(data.frame(), file.path(res, "KEGG_table.csv"), row.names = FALSE)
  }

  status$input_genes <- length(sig_genes)
  status$universe <- cfg$enrichment$universe
  write_json(file.path(res, "enrichment_status.json"), status)
  log_info("已生成 enrichment_status.json")
  invisible(NULL)
}

#' 无显著基因时统一写空结果
write_empty_enrichment <- function(cfg, status, reason) {
  res <- cfg$output$results_dir
  utils::write.csv(data.frame(), file.path(res, "GO_table.csv"), row.names = FALSE)
  utils::write.csv(data.frame(), file.path(res, "KEGG_table.csv"), row.names = FALSE)
  status$go <- list(status = "skipped", reason = reason)
  status$kegg <- list(status = "skipped", reason = reason)
  write_json(file.path(res, "enrichment_status.json"), status)
}

#' dotplot，展示 top N 条目
make_dotplot <- function(x, cfg, title) {
  n <- min(cfg$enrichment$top_terms, nrow(as.data.frame(x)))
  p <- enrichplot::dotplot(x, showCategory = n) +
    labs(title = title, subtitle = sprintf("top %d terms, p.adjust < %g", n,
                                           cfg$enrichment$pvalue_cutoff)) +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_text(size = 7))
  p
}

if (!GEO_ORCHESTRATED()) {
  cfg <- load_config()
  run_04a_heatmap(cfg)
  run_04b_enrichment(cfg)
}
