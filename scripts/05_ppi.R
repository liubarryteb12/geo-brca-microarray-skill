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

#' 下载并读取 STRING 的 protein.links 文件
#'
#' 刻意不用 STRINGdb::get_interactions()。实测（run 35409119340）：v12.0 上
#' 1321/1423 个基因映射成功、links 文件也下好了，但 get_interactions() **返回 0 行
#' 且不报错**，流程被静默推进到共表达回退分支。自己读文件，行为可预期。
fetch_string_links <- function(version, threshold, cache_dir) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  # score_threshold >= 400 时 STRING 提供预过滤的 min{N} 流式文件，比全量小得多
  use_min <- threshold >= 400
  fname <- if (use_min) {
    sprintf("9606.protein.links.v%s.min%d.txt.gz", version, threshold)
  } else {
    sprintf("9606.protein.links.v%s.txt.gz", version)
  }
  path <- file.path(cache_dir, fname)
  if (!file.exists(path)) {
    url <- if (use_min) {
      sprintf("https://stringdb-downloads.org/download/stream/protein.links.v%s/%s", version, fname)
    } else {
      sprintf("https://stringdb-downloads.org/download/protein.links.v%s/%s", version, fname)
    }
    log_info(sprintf("下载 STRING links 文件: %s", url))
    utils::download.file(url, path, mode = "wb", quiet = TRUE)
  }
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  utils::read.delim(con, sep = " ", stringsAsFactors = FALSE)
}

run_05_ppi <- function(cfg) {
  log_info("=== 步骤 05：差异基因互作网络 ===")
  ensure_dirs(cfg)

  res <- cfg$output$results_dir
  deg <- utils::read.csv(file.path(res, "deg_table.csv"), stringsAsFactors = FALSE)
  # 与富集用同一套挑选逻辑：FDR 不够时退回 raw P 排序前 N 个，并标注模式
  sel <- select_degs(deg, cfg, min_genes = MIN_GENES_FOR_NETWORK)
  sig <- unique(sel$genes[!is.na(sel$genes) & nzchar(sel$genes)])

  status <- list(species = 9606, score_threshold = cfg$thresholds$string_score,
                 deg_mode = sel$mode, deg_reason = sel$reason)
  if (identical(sel$mode, "ranked_fallback")) {
    log_warn(sprintf("无基因通过 FDR，PPI 改用 raw P 排序前 %d 个基因（假设生成，非显著 DEG）",
                     length(sig)))
  }

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

      # STRING 会随年份下线旧版本；按新到旧依次尝试，避免写死一个版本后静默失效
      sdb <- NULL
      last_err <- NULL
      chosen <- NULL
      for (v in c("12.0", "11.5", "11.0")) {
        sdb <- tryCatch(
          STRINGdb::STRINGdb$new(version = v, species = 9606,
                                 score_threshold = cfg$thresholds$string_score,
                                 input_directory = input_dir),
          error = function(e) { last_err <<- conditionMessage(e); NULL }
        )
        if (!is.null(sdb)) { chosen <- v; log_info(sprintf("STRING 数据库版本 %s 可用", v)); break }
        log_warn(sprintf("STRING 版本 %s 不可用: %s", v, last_err))
      }
      if (is.null(sdb)) stop(sprintf("STRING 各版本均不可用: %s", last_err))

      mapped <- sdb$map(data.frame(gene = sig), "gene", removeUnmappedRows = TRUE)
      if (nrow(mapped) < MIN_GENES_FOR_NETWORK) {
        stop(sprintf("STRING 仅映射到 %d 个基因", nrow(mapped)))
      }
      log_info(sprintf("STRING 映射成功: %d / %d 个基因", nrow(mapped), length(sig)))

      # 只保留两端都在输入基因里的互作
      links <- fetch_string_links(chosen, cfg$thresholds$string_score, input_dir)
      keep <- links$protein1 %in% mapped$STRING_id & links$protein2 %in% mapped$STRING_id
      links <- links[keep, , drop = FALSE]
      if (nrow(links) == 0L) {
        stop(sprintf("STRING v%s 的 links 中，没有两端都落在输入基因内的互作", chosen))
      }

      # STRING 的 id 形如 "9606.ENSP00000..."，还原成基因 symbol
      id2sym <- stats::setNames(mapped$gene, mapped$STRING_id)
      edges <- data.frame(
        from_gene = unname(id2sym[links$protein1]),
        to_gene   = unname(id2sym[links$protein2]),
        score     = as.numeric(links$combined_score),
        stringsAsFactors = FALSE
      )
      edges <- unique(edges[!is.na(edges$from_gene) & !is.na(edges$to_gene), , drop = FALSE])
      log_info(sprintf("STRING 互作边: %d 条（阈值 %d，STRING v%s）",
                       nrow(edges), cfg$thresholds$string_score, chosen))

      g <- igraph::graph_from_data_frame(edges[, c("from_gene", "to_gene")], directed = FALSE)
      status$string_version <- chosen
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
  # Fruchterman-Reingold 布局要求权重非负，而 Spearman r 可以为负。
  # 实测：直接把带符号的 r 当权重会让 layout_with_fr 报
  # "Weights must be positive for Fruchterman-Reingold layout" 并中断整个步骤。
  # 所以布局权重取 |r|，有符号的 r 单独留在边表里。
  el$signed_r <- el$weight
  igraph::E(g)$weight <- abs(el$weight)
  el$score <- round(abs(el$signed_r) * 1000)   # 与 STRING 分数同量纲
  status$coexpression_threshold <- thr
  status$note <- "布局权重使用 |Spearman r|；边表保留有符号的 signed_r"
  write_ppi_outputs(cfg, g, el[, c("from_gene", "to_gene", "score", "signed_r")],
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

  # 状态先落盘：出图是最后一步，画不出来也不能丢掉"网络是否构建成功"这个事实。
  # 实测 run 35409119340 就是先出图后写状态，图挂了连 ppi_status.json 都没有。
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

  # 网络图：节点大小按 degree，颜色按 degree 深浅
  png_path <- file.path(res, "PPI_network.png")
  plot_err <- tryCatch({
    open_png(png_path, width = 1600, height = 1400, res = 150)
    graphics::par(mar = c(1, 1, 3, 1))
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
    grDevices::dev.off()
    NULL
  }, error = function(e) {
    if (grDevices::dev.cur() > 1L) try(grDevices::dev.off(), silent = TRUE)
    conditionMessage(e)
  })
  if (is.null(plot_err)) {
    status$plot <- basename(png_path)
  } else {
    status$plot_error <- plot_err
    log_warn(sprintf("网络图绘制失败（网络本身已构建成功，边表与 hub 基因不受影响）: %s", plot_err))
  }

  write_json(file.path(res, "ppi_status.json"), status)

  log_info(sprintf("已生成 hub_genes.csv / ppi_edges.csv（%s，%d 节点 %d 边）",
                   method, igraph::vcount(g), igraph::ecount(g)))
  log_info(sprintf("Hub 基因: %s", paste(hub_df$gene, collapse = ", ")))
  invisible(hub_df)
}

if (!GEO_ORCHESTRATED()) {
  cfg <- load_config()
  run_05_ppi(cfg)
}
