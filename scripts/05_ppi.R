# ============================================================================
# 05_ppi.R — 差异基因互作网络（STRING PPI）
# ============================================================================
# spec 的 ppi_string 步骤。
#
# 主路径：STRINGdb（species=9606, score >= 400），hub 基因 = degree 前 N。
# 回退路径：STRING 网络文件下载失败时，改用表达相关性构建共表达网络
#           （igraph 本地计算），并在 ppi_status.json 中**明确标注为共表达而非 PPI**。
#
# 输出：results/PPI_network.png
#       results/hub_genes.csv
#       results/ppi_edges.csv
#       results/ppi_status.json
# ============================================================================

suppressPackageStartupMessages({
  library(igraph)
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

MIN_GENES_FOR_NETWORK <- 5L

#' 打开 PNG 设备；无头 Linux runner 上优先用 cairo，缺失时退回默认设备
open_png <- function(path, ...) {
  if (isTRUE(capabilities("cairo"))) {
    grDevices::png(path, ..., type = "cairo")
  } else {
    log_warn("R 未编译 cairo 支持，使用默认 PNG 设备")
    grDevices::png(path, ...)
  }
}

run_05_ppi <- function(cfg) {
  log_info("=== 步骤 05：差异基因互作网络 ===")
  ensure_dirs(cfg)

  res <- cfg$output$results_dir
  deg <- utils::read.csv(file.path(res, "deg_table.csv"), stringsAsFactors = FALSE)
  sig <- deg$gene[deg$adj.P.Val < cfg$thresholds$adj_p & abs(deg$logFC) > cfg$thresholds$log2fc]
  sig <- unique(sig[!is.na(sig) & nzchar(sig)])

  status <- list(species = 9606, score_threshold = cfg$thresholds$string_score)

  if (length(sig) < MIN_GENES_FOR_NETWORK) {
    msg <- sprintf("显著 DEG 仅 %d 个（< %d），无法构建有意义的互作网络",
                   length(sig), MIN_GENES_FOR_NETWORK)
    log_warn(msg)
    status$status <- "skipped"
    status$reason <- msg
    status$n_input_genes <- length(sig)
    write_json(file.path(res, "ppi_status.json"), status)
    utils::write.csv(data.frame(), file.path(res, "hub_genes.csv"), row.names = FALSE)
    utils::write.csv(data.frame(), file.path(res, "ppi_edges.csv"), row.names = FALSE)
    return(invisible(NULL))
  }
  status$n_input_genes <- length(sig)
  log_info(sprintf("PPI 输入基因数: %d", length(sig)))

  # 探针 ID 查不了 STRING，直接走共表达回退并写明原因
  feat_path <- file.path(cfg$output$data_dir, "feature_mode.json")
  feat <- if (file.exists(feat_path)) {
    tryCatch(jsonlite::fromJSON(feat_path, simplifyVector = FALSE), error = function(e) NULL)
  } else NULL
  probe_mode <- !is.null(feat) && !identical(feat$mode, "symbol")
  if (probe_mode) {
    log_warn(sprintf("特征为探针 ID 而非基因 symbol（%s），跳过 STRING 查询，直接构建共表达网络",
                     feat$reason %||% "原因未知"))
    status$string_skipped <- "probe-level features cannot be queried against STRING"
  }

  # ---- 主路径：STRINGdb ---------------------------------------------------
  string_ok <- FALSE
  if (!probe_mode && require_pkg("STRINGdb")) {
    string_ok <- tryCatch({
      input_dir <- file.path(cfg$output$data_dir, "string_cache")
      if (!dir.exists(input_dir)) dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
      log_info("连接 STRINGdb（首次运行需下载网络文件，约 100 MB）...")

      sdb <- STRINGdb::STRINGdb$new(
        version = "11.5", species = 9606,
        score_threshold = cfg$thresholds$string_score,
        input_directory = input_dir
      )
      mapped <- sdb$map(data.frame(gene = sig), "gene", removeUnmappedRows = TRUE)
      if (nrow(mapped) < MIN_GENES_FOR_NETWORK) {
        stop(sprintf("STRING 仅映射到 %d 个基因", nrow(mapped)))
      }
      log_info(sprintf("STRING 映射成功: %d / %d 个基因", nrow(mapped), length(sig)))

      edges <- sdb$get_interactions(mapped$STRING_id)
      edges <- edges[edges$score >= cfg$thresholds$string_score, , drop = FALSE]
      if (nrow(edges) == 0L) stop("STRING 未返回任何达到阈值的互作")

      # STRING 的 id 形如 "9606.ENSP00000..."，还原成基因 symbol
      id2sym <- stats::setNames(mapped$gene, mapped$STRING_id)
      edges$from_gene <- id2sym[edges$from]
      edges$to_gene   <- id2sym[edges$to]
      edges <- edges[!is.na(edges$from_gene) & !is.na(edges$to_gene), , drop = FALSE]
      edges <- unique(edges[, c("from_gene", "to_gene", "score")])
      log_info(sprintf("STRING 互作边: %d 条", nrow(edges)))

      g <- igraph::graph_from_data_frame(edges[, c("from_gene", "to_gene")], directed = FALSE)
      write_ppi_outputs(cfg, g, edges, "string_ppi", status)
      string_ok <- TRUE
      TRUE
    }, error = function(e) {
      log_warn(sprintf("STRINGdb 失败，转入共表达回退: %s", conditionMessage(e)))
      status$string_error <<- conditionMessage(e)
      FALSE
    })
  }

  if (isTRUE(string_ok)) return(invisible(NULL))

  # ---- 回退路径：共表达网络 ----------------------------------------------
  log_warn("使用 igraph 共表达网络回退（这是表达相关性，不是物理蛋白互作）")
  expr <- readRDS(file.path(cfg$output$data_dir, "expr_clean.rds"))
  genes <- intersect(sig, rownames(expr))
  if (length(genes) < MIN_GENES_FOR_NETWORK) {
    status$status <- "failed"
    status$reason <- "STRING 失败且可用于共表达分析的基因不足"
    write_json(file.path(res, "ppi_status.json"), status)
    utils::write.csv(data.frame(), file.path(res, "hub_genes.csv"), row.names = FALSE)
    utils::write.csv(data.frame(), file.path(res, "ppi_edges.csv"), row.names = FALSE)
    return(invisible(NULL))
  }

  # 基因数可能上千，先按 DEG 显著性取前 200 个，控制相关矩阵规模
  genes <- head(genes, 200L)
  cor_mat <- stats::cor(t(expr[genes, , drop = FALSE]), method = "spearman")
  cor_mat[is.na(cor_mat)] <- 0
  thr <- 0.9
  adj <- cor_mat
  adj[abs(adj) < thr] <- 0
  diag(adj) <- 0
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)
  g <- igraph::simplify(g)

  if (igraph::ecount(g) == 0L) {
    status$status <- "failed"
    status$reason <- sprintf("STRING 失败，且共表达网络中无 |Spearman r| >= %.2f 的边", thr)
    write_json(file.path(res, "ppi_status.json"), status)
    utils::write.csv(data.frame(), file.path(res, "hub_genes.csv"), row.names = FALSE)
    utils::write.csv(data.frame(), file.path(res, "ppi_edges.csv"), row.names = FALSE)
    return(invisible(NULL))
  }

  el <- igraph::as_data_frame(g, what = "edges")
  names(el)[names(el) == "from"] <- "from_gene"
  names(el)[names(el) == "to"]   <- "to_gene"
  el$score <- round(el$weight * 1000)   # 与 STRING 分数同量纲，便于下游统一处理
  status$coexpression_threshold <- thr
  write_ppi_outputs(cfg, g, el[, c("from_gene", "to_gene", "score")],
                    "coexpression_fallback", status)
  invisible(NULL)
}

#' 统一的网络落盘：图、hub 基因、边表、状态
write_ppi_outputs <- function(cfg, g, edges, method, status) {
  res <- cfg$output$results_dir

  deg_all <- igraph::degree(g)
  hub_n <- min(cfg$analysis$hub_gene_count, length(deg_all))
  hub <- sort(deg_all, decreasing = TRUE)[seq_len(hub_n)]

  hub_df <- data.frame(gene = names(hub), degree = as.integer(hub), method = method)
  utils::write.csv(hub_df, file.path(res, "hub_genes.csv"), row.names = FALSE)
  utils::write.csv(edges, file.path(res, "ppi_edges.csv"), row.names = FALSE)

  # 网络图：节点大小按 degree，颜色按 degree 深浅
  open_png(file.path(res, "PPI_network.png"), width = 1600, height = 1400, res = 150)
  op <- graphics::par(mar = c(1, 1, 3, 1))
  on.exit({ graphics::par(op); grDevices::dev.off() }, add = TRUE)

  igraph::plot.igraph(
    g,
    layout = igraph::layout_with_fr(g),
    vertex.size = pmin(4 + deg_all * 1.5, 18),
    vertex.color = grDevices::colorRampPalette(c("#AFC7E3", "#C1443C"))(max(deg_all) + 1)[deg_all + 1],
    vertex.frame.color = "white",
    vertex.label = ifelse(deg_all >= stats::quantile(deg_all, 0.8), names(deg_all), NA),
    vertex.label.cex = 0.7,
    vertex.label.color = "black",
    edge.color = grDevices::adjustcolor("grey40", alpha.f = 0.4),
    edge.width = 0.8,
    main = sprintf("%s network - %s\nnodes=%d edges=%d",
                   if (identical(method, "string_ppi")) "STRING PPI" else "Co-expression (FALLBACK)",
                   cfg$dataset_id, igraph::vcount(g), igraph::ecount(g))
  )

  status$status <- if (identical(method, "string_ppi")) "ok" else "fallback"
  status$method <- method
  status$nodes <- igraph::vcount(g)
  status$edges <- igraph::ecount(g)
  status$hub_genes <- hub_df$gene
  if (!identical(method, "string_ppi")) {
    status$caveat <- paste(
      "这不是物理蛋白互作网络。STRINGdb 不可用，回退为基于表达谱的共表达网络",
      "(Spearman |r| 阈值见 coexpression_threshold)。不得当作 PPI 证据引用。"
    )
  }
  write_json(file.path(res, "ppi_status.json"), status)

  log_info(sprintf("已生成 PPI_network.png / hub_genes.csv / ppi_edges.csv（%s，%d 节点 %d 边）",
                   method, igraph::vcount(g), igraph::ecount(g)))
  log_info(sprintf("Hub 基因: %s", paste(hub_df$gene, collapse = ", ")))
  invisible(hub_df)
}

if (!GEO_ORCHESTRATED()) {
  cfg <- load_config()
  run_05_ppi(cfg)
}
