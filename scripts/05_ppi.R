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
      # **必须显式补上边权。** 上面只传了两列，图里没有 weight 属性；
      # 而出图代码要按边强度筛选、要按权重调透明度，拿到 NULL 会直接报错
      # （实测踩过两次：`edf$weight / max(edf$weight)` 得到长度 0 的向量 ->
      # "replacement has 0 rows"；`order(NULL)` -> 全 NA 的边索引 ->
      # "argument 1 is not a vector"）。两次都让 STRING 静默退化成共表达网络。
      igraph::E(g)$weight <- as.numeric(edges$score)
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

  # ---- 网络图 --------------------------------------------------------------
  #
  # **旧版为什么丑：** 390 节点 / 3284 条边全部画出来是一团毛线，
  # 而且按 degree 上色（连续深浅）看不出任何结构 —— 读者只能看到"中间密、边上稀"。
  # 标签取 degree 前 20%，约 78 个，在毛线球上互相压成一团。
  #
  # 改法：
  #   1. 只取**最大连通分量**（零散小碎片对"互作网络"没有信息量）
  #   2. 节点仍超过上限时按 degree 取前 N 个 —— 这才是毛线球的根因
  #   3. **Louvain 社区着色**，用分类色板。有结构可看，而不是一片渐变色
  #   4. 边按权重调透明度，弱边近乎消失，强边浮现
  #   5. 只标注 top hub，用 ggrepel 避免重叠
  #   6. 用 ggplot2 画而不是 plot.igraph：布局坐标只算一次，
  #      PDF 与 PNG 共用同一份坐标，天然可复现
  plot_ppi_network <- function(cfg, g, method, status) {
    res <- cfg$output$results_dir
    png_path <- file.path(res, "PPI_network.png")

    g_full <- g
    # 1. 最大连通分量
    comps <- igraph::components(g_full)
    if (length(comps$csize) > 1L) {
      g <- igraph::induced_subgraph(g_full, which(comps$membership == which.max(comps$csize)))
      log_info(sprintf("网络图：%d 个连通分量，取最大的一个（%d 节点）",
                       length(comps$csize), igraph::vcount(g)))
    }
    # 2. 节点上限
    #
    # **注意：只按 degree 截节点是反效果的。** 实测取 degree 前 200 个节点后，
    # 边数反而从 3284 涨到 4906 —— 因为高 degree 的节点彼此高度互联，
    # 取它们等于取网络最密的核，比全图还乱。
    # 所以真正的密度控制要落在**边**上：按置信度保留最强的若干条。
    max_nodes <- cfg$analysis$ppi_plot_max_nodes
    if (is.null(max_nodes) || max_nodes <= 0) max_nodes <- 150L
    max_edges <- cfg$analysis$ppi_plot_max_edges
    if (is.null(max_edges) || max_edges <= 0) max_edges <- 700L

    n_before <- igraph::vcount(g)
    e_before <- igraph::ecount(g)
    if (n_before > max_nodes) {
      d <- igraph::degree(g)
      keep <- names(sort(d, decreasing = TRUE))[seq_len(max_nodes)]
      g <- igraph::induced_subgraph(g, keep)
    }
    # 3. 按边权保留最强的 max_edges 条，再丢掉因此变成孤立的节点
    #
    # **先确认图真的有边权。** 没有 weight 时 `order(NULL)` 会返回 integer(0)，
    # 再取 `[seq_len(max_edges)]` 得到一整条 NA，传给 subgraph.edges 直接报
    # "argument 1 is not a vector"。宁可跳过滤，也不要静默把整个 STRING 路径
    # 打进回退分支 —— 那会让"PPI 网络"变成共表达网络而不自知。
    w_all <- igraph::E(g)$weight
    has_weight <- !is.null(w_all) && length(w_all) == igraph::ecount(g)
    if (!has_weight) {
      log_warn("图没有边权属性，跳过按强度筛选（保留全部边）")
    } else if (igraph::ecount(g) > max_edges) {
      strong <- order(w_all, decreasing = TRUE)[seq_len(max_edges)]
      g <- igraph::subgraph.edges(g, strong, delete.vertices = TRUE)
    }
    log_info(sprintf("网络图：%d 节点 / %d 边 → 过滤后 %d 节点 / %d 边",
                     n_before, e_before, igraph::vcount(g), igraph::ecount(g)))
    if (igraph::vcount(g) < 2L || igraph::ecount(g) == 0L) {
      status$plot_error <- "过滤后网络为空，跳过绘图"
      log_warn(status$plot_error)
      return(status)
    }

    seed <- cfg$analysis$seed
    if (!is.null(seed)) set.seed(seed)
    # 4. Louvain 社区（随机算法，必须设种子）
    #
    # 注意全部用 igraph:: 前缀 —— 本仓库不 attach 任何包（没有 library() 调用），
    # 裸 V()/E() 会 "could not find function"
    #
    # **只给最大的若干个模块上色，其余归入灰色 "other"。**
    # 实测 129 节点的过滤网络在 resolution=1 下切出 **11 个**社区，
    # 而经计算验证的色盲安全色板只有 4 色（见 common.R 的 pal_categorical）。
    # 11 种颜色必然走插值降级，插出来的颜色没经过验证、彼此也分不开 ——
    # 那样的图看着花花绿绿，实际读不出结构。
    comm <- igraph::cluster_louvain(g)
    memb <- as.character(comm$membership)
    n_comm_all <- length(unique(memb))
    max_mod <- cfg$analysis$ppi_plot_modules
    if (is.null(max_mod) || max_mod <= 0) max_mod <- 4L
    sizes <- sort(table(memb), decreasing = TRUE)
    shown <- names(sizes)[seq_len(min(max_mod, length(sizes)))]
    memb[!(memb %in% shown)] <- "other"
    # 按模块大小定 levels，保证配色稳定（不随 Louvain 的编号跳变）
    lv <- c(shown, if (any(memb == "other")) "other")
    igraph::V(g)$community <- factor(memb, levels = lv)
    n_comm <- length(lv)
    if (n_comm_all > length(shown)) {
      log_info(sprintf("网络图：Louvain 切出 %d 个模块，只给最大的 %d 个上色，其余 %d 个节点归入 other",
                       n_comm_all, length(shown), sum(memb == "other")))
    }

    vdf <- data.frame(
      name = igraph::V(g)$name,
      degree = as.integer(igraph::degree(g)),
      community = igraph::V(g)$community,
      stringsAsFactors = FALSE
    )

    # 5. 布局：**同心圆环**，不是力导向
    #
    # 力导向（Fruchterman-Reingold）在这张图上是失败的：390 节点挤成一团，
    # 中间密到看不出结构、边上又空着，而且**每次运行布局都不一样**（另一个随机源）。
    #
    # 同心圆环把"谁是 hub"直接编码成半径：内圈 = degree 最高的核心，
    # 外圈 = 边缘基因。读者不用找中心，一眼就知道层次。
    #
    # 三条实现要点：
    #   * **每环节点数按半径成比例**（周长 ∝ 半径），否则内圈挤成一坨、外圈稀稀拉拉
    #   * **环内按社区排序**（同一模块占同一角度扇区）—— 这样模块内的边是短弦，
    #     模块间的边才跨圆心，边交叉大幅减少。按 degree 排会让每条边都横穿全图
    #   * 半径等距递增，配一圈很淡的参考圆，让"几圈"这件事看得见
    n_rings <- cfg$analysis$ppi_plot_rings
    if (is.null(n_rings) || n_rings < 1) n_rings <- 3L
    n_rings <- as.integer(min(n_rings, max(1L, floor(nrow(vdf) / 12L))))
    radii <- 1 + 0.72 * (seq_len(n_rings) - 1L)
    # 按半径比例分配每环节点数，再修正取整误差
    sizes_r <- pmax(1L, as.integer(round(nrow(vdf) * radii / sum(radii))))
    while (sum(sizes_r) > nrow(vdf)) sizes_r[which.max(sizes_r)] <- sizes_r[which.max(sizes_r)] - 1L
    while (sum(sizes_r) < nrow(vdf)) sizes_r[n_rings] <- sizes_r[n_rings] + 1L
    # 上面两个循环可能把某一环减到 0，那样 idx[pos:(pos-1)] 会取到倒序下标。
    # 节点数很少时才可能发生（n_rings 已被 n/12 限制过），兜一下。
    if (any(sizes_r < 1L)) {
      sizes_r <- rep(1L, n_rings)
      sizes_r[n_rings] <- nrow(vdf) - (n_rings - 1L)
    }

    # 先按 degree 降序决定"谁在内圈"
    ring_of <- integer(nrow(vdf))
    idx <- order(-vdf$degree, vdf$name)
    pos <- 1L
    for (r in seq_len(n_rings)) {
      ring_of[idx[pos:(pos + sizes_r[r] - 1L)]] <- r
      pos <- pos + sizes_r[r]
    }

    # 环内按社区排（社区之间按模块大小），使同一模块落在同一角度扇区
    comm_levels <- levels(vdf$community)
    comm_rank <- match(vdf$community, comm_levels)
    vdf$x <- NA_real_; vdf$y <- NA_real_
    for (r in seq_len(n_rings)) {
      mem <- which(ring_of == r)
      mem <- mem[order(comm_rank[mem], -vdf$degree[mem], vdf$name[mem])]
      k <- length(mem)
      ang <- 2 * pi * (seq_len(k) - 1L) / k + pi / 2   # 从正上方开始
      vdf$x[mem] <- radii[r] * cos(ang)
      vdf$y[mem] <- radii[r] * sin(ang)
    }
    vdf$ring <- factor(sprintf("ring %d", ring_of),
                       levels = sprintf("ring %d", seq_len(n_rings)))

    # **布局落盘。** 图本身看不出"第 3 环是不是真的在外圈"，光看 PNG 只能靠猜；
    # 写下每个节点的环号与坐标，环结构就是可核对的数据而不是视觉效果。
    # 同时它让"哪些基因属于核心环"成为可引用的结果。
    utils::write.csv(
      data.frame(gene = vdf$name, ring = ring_of, radius = radii[ring_of],
                 angle_deg = round((atan2(vdf$y, vdf$x) * 180 / pi + 360) %% 360, 1),
                 x = round(vdf$x, 4), y = round(vdf$y, 4),
                 degree = vdf$degree, module = as.character(vdf$community),
                 stringsAsFactors = FALSE),
      file.path(res, "ppi_plot_layout.csv"), row.names = FALSE, fileEncoding = "UTF-8")

    edf <- igraph::as_data_frame(g, what = "edges")
    edf$x    <- vdf$x[match(edf$from, vdf$name)]
    edf$y    <- vdf$y[match(edf$from, vdf$name)]
    edf$xend <- vdf$x[match(edf$to,   vdf$name)]
    edf$yend <- vdf$y[match(edf$to,   vdf$name)]
    edf$w    <- if (has_weight) edf$weight / max(edf$weight) else 0.5

    # 参考圆：让"几圈"看得见，颜色压到几乎不可见，不与数据争视觉
    ring_path <- do.call(rbind, lapply(seq_len(n_rings), function(r) {
      th <- seq(0, 2 * pi, length.out = 240)
      data.frame(x = radii[r] * cos(th), y = radii[r] * sin(th),
                 grp = factor(r), stringsAsFactors = FALSE)
    }))

    # 6. 只标注 top hub
    hub_k <- min(20L, nrow(vdf))
    lab <- vdf[order(-vdf$degree)[seq_len(hub_k)], , drop = FALSE]

    # 上色：只有被展示的模块用验证过的分类色，other 用中性灰
    shown_lv <- setdiff(levels(vdf$community), "other")
    comm_cols <- stats::setNames(pal_categorical(length(shown_lv)), shown_lv)
    if ("other" %in% levels(vdf$community)) comm_cols["other"] <- PAL$ns
    n_comm <- length(shown_lv)

    p <- ggplot2::ggplot() +
      ggplot2::geom_path(
        data = ring_path,
        ggplot2::aes(x = x, y = y, group = grp),
        colour = PAL$grid, linewidth = 0.25) +
      ggplot2::geom_segment(
        data = edf,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend, alpha = w),
        colour = PAL$edge, linewidth = 0.22) +
      ggplot2::scale_alpha_continuous(range = c(0.03, 0.5), guide = "none") +
      ggplot2::geom_point(
        data = vdf,
        ggplot2::aes(x = x, y = y, size = degree, fill = community),
        shape = 21, colour = "white", stroke = 0.35) +
      ggplot2::scale_fill_manual(
        values = comm_cols, name = "module",
        labels = stats::setNames(
          c(sprintf("%s (%d)", shown_lv, as.integer(sizes[shown_lv])),
            if ("other" %in% levels(vdf$community))
              sprintf("other (%d)", sum(memb == "other"))),
          c(shown_lv, if ("other" %in% levels(vdf$community)) "other")),
        guide = "legend") +
      ggplot2::scale_size_continuous(name = "degree", range = c(1.5, 6.5),
                                     breaks = pretty(range(vdf$degree), 4)) +
      ggrepel::geom_text_repel(
        data = lab, ggplot2::aes(x = x, y = y, label = name),
        size = 2.4, colour = PAL$ink, fontface = "bold",
        segment.size = 0.2, segment.colour = PAL$muted,
        min.segment.length = 0, max.overlaps = Inf, box.padding = 0.35,
        # 显式播种：不传时 ggrepel 用环境 RNG，位置会随上游随机数消耗量漂移。
        # seed 默认值是 NA（不是 NULL），所以这里要转换。
        seed = if (is.null(seed)) NA else seed) +
      ggplot2::labs(
        title = sprintf("%s network - %s",
                        if (identical(method, "string_ppi")) "STRING PPI" else "Co-expression (FALLBACK)",
                        cfg$dataset_id),
        subtitle = sprintf(paste0("%d nodes / %d edges shown, laid out on %d concentric rings: ",
                                  "inner ring = highest degree (%d hubs labelled). ",
                                  "Nodes ordered by degree, then grouped by module so each module ",
                                  "occupies one angular sector; node colour = Louvain module ",
                                  "(%d found, top %d coloured, rest grey), size = degree, ",
                                  "edge opacity = interaction confidence. ",
                                  "Filtered for readability - full network in ppi_edges.csv (%d nodes, %d edges)."),
                           igraph::vcount(g), igraph::ecount(g), n_rings,
                           nrow(lab), n_comm_all, length(shown),
                           igraph::vcount(g_full), igraph::ecount(g_full)),
        x = NULL, y = NULL) +
      ggplot2::coord_fixed() +
      ggplot2::theme_void(base_size = 10) +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", size = 11),
        plot.subtitle = ggplot2::element_text(colour = PAL$muted, size = 8),
        legend.position = "right",
        plot.margin = ggplot2::margin(8, 8, 8, 8))

    plot_err <- tryCatch({
      save_pdf(file.path(res, "PPI_network.pdf"), print(p), width = 10, height = 8.5)
      NULL
    }, error = function(e) conditionMessage(e))

    if (is.null(plot_err)) {
      status$plot <- "PPI_network.png"
      status$plot_nodes <- igraph::vcount(g)
      status$plot_edges <- igraph::ecount(g)
      status$plot_modules <- n_comm
      status$plot_modules_total <- n_comm_all
      status$plot_rings <- n_rings
      status$plot_ring_sizes <- as.integer(sizes_r)
      status$plot_layout <- "concentric_rings"
      status$plot_filtered <- igraph::vcount(g) < igraph::vcount(g_full)
      status$plot_note <- sprintf(
        paste0("图为可读性做过过滤：最大连通分量 → degree 前 %d 个节点 → 最强的 %d 条边；",
               "布局为 %d 个同心圆环（内圈 = degree 最高），环内按社区排序。",
               "完整网络见 ppi_edges.csv（%d 节点 %d 边）。"),
        max_nodes, max_edges, n_rings, igraph::vcount(g_full), igraph::ecount(g_full))
      log_info(sprintf("已生成 PPI_network.png（%d 节点 / %d 边 / %d 个模块 / %d 个同心环 %s）",
                       igraph::vcount(g), igraph::ecount(g), n_comm, n_rings,
                       paste(sizes_r, collapse = "-")))
    } else {
      status$plot_error <- plot_err
      log_warn(sprintf("网络图绘制失败（网络本身已构建成功，边表与 hub 基因不受影响）: %s", plot_err))
    }
    status
  }

  # **绘图必须与网络构建隔离。** 出图代码里的任何错误如果逃逸出去，会被
  # 调用方的 `tryCatch` 当成"STRING 失败"接住，于是整条路径静默退化成共表达网络 ——
  # 一个画图的 bug 悄悄换掉了分析方法，而日志里只留一句"STRINGdb 失败"。
  # 实测正是如此：`E(g)$weight` 缺失导致出图报错，连续两轮 CI 的 PPI 都变成了
  # 共表达网络，而 ppi_status.json 里看起来"有图、有 hub 基因"，一切正常。
  #
  # 所以这里自己兜住：出图失败只影响图，方法本身（string_ppi）如实记录。
  status <- tryCatch(
    plot_ppi_network(cfg, g, method, status),
    error = function(e) {
      status$plot_error <- conditionMessage(e)
      log_warn(sprintf("网络图绘制失败（网络本身已构建成功，边表与 hub 基因不受影响）: %s",
                       conditionMessage(e)))
      status
    }
  )

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
