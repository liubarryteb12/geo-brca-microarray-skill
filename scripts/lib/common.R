# ============================================================================
# lib/common.R — 共享工具：配置加载、日志、JSON 状态、GEO SOFT 抓取
# ============================================================================
# 被 scripts/00..05 与 main_analysis.R 共同 source()。
# 不依赖任何 Bioconductor 包，只用 base R + jsonlite + yaml。
# ============================================================================

# 编排器标志：main_analysis.R 会设为 TRUE，脚本据此决定是否自动执行 main()
GEO_ORCHESTRATED <- function() isTRUE(getOption("geo.orchestrated", FALSE))

# ---- 日志 ------------------------------------------------------------------

.geo_log <- function(level, ...) {
  msg <- paste0(..., collapse = "")
  cat(sprintf("[%s] %-5s %s\n", format(Sys.time(), "%H:%M:%S"), level, msg))
  flush(stdout())
}

log_info <- function(...) .geo_log("INFO", ...)
log_warn <- function(...) .geo_log("WARN", ...)
log_error <- function(...) .geo_log("ERROR", ...)

# ---- 配置 ------------------------------------------------------------------

#' 解析 `--config <path>` 命令行参数
#' @return 配置文件路径（默认 assets/config.yml）
parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  cfg <- "assets/config.yml"
  i <- 1L
  while (i <= length(argv)) {
    if (argv[i] %in% c("--config", "-c")) {
      if (i == length(argv)) stop("--config 需要一个路径参数")
      cfg <- argv[i + 1L]
      i <- i + 2L
    } else if (argv[i] %in% c("--help", "-h")) {
      cat("用法: Rscript <script>.R [--config assets/config.yml]\n")
      quit(save = "no", status = 0L)
    } else {
      # 允许直接传位置参数作为配置路径
      cfg <- argv[i]
      i <- i + 1L
    }
  }
  cfg
}

#' 读取并校验配置文件
#' @param path 配置 YAML 路径
#' @return 配置 list（已做必填项与取值合法性校验）
load_config <- function(path = parse_args()) {
  if (!file.exists(path)) stop(sprintf("配置文件不存在: %s", path))
  cfg <- yaml::read_yaml(path)

  required <- c("dataset_id", "group_field", "group_values", "contrast")
  missing <- setdiff(required, names(cfg))
  if (length(missing) > 0L) {
    stop(sprintf("配置缺少必填字段: %s", paste(missing, collapse = ", ")))
  }
  if (length(cfg$contrast) != 2L) {
    stop("contrast 必须是长度为 2 的向量: [分子, 分母]")
  }
  unknown <- setdiff(unlist(cfg$contrast), names(cfg$group_values))
  if (length(unknown) > 0L) {
    stop(sprintf("contrast 引用了 group_values 中不存在的组: %s", paste(unknown, collapse = ", ")))
  }
  if (!is.null(cfg$group_values) && length(cfg$group_values) < 2L) {
    stop("group_values 至少需要两个组才能做对比")
  }

  # 默认值，避免下游到处写 is.null 判断
  cfg$platform_id        <- cfg$platform_id %||% ""
  cfg$paired             <- isTRUE(cfg$paired)
  cfg$thresholds         <- utils::modifyList(
    list(adj_p = 0.05, log2fc = 1.0, string_score = 400, outlier_cor = 0.8),
    cfg$thresholds %||% list()
  )
  cfg$analysis           <- utils::modifyList(
    list(top_heatmap_genes = 50, impute_k = 10, pca_top_genes = 2000, hub_gene_count = 10),
    cfg$analysis %||% list()
  )
  cfg$enrichment         <- utils::modifyList(
    list(universe = "genome", p_adjust = "BH", pvalue_cutoff = 0.05,
         qvalue_cutoff = 0.2, top_terms = 15, ont = "BP", kegg_organism = "hsa"),
    cfg$enrichment %||% list()
  )
  cfg$output             <- utils::modifyList(
    list(results_dir = "results", data_dir = "data"),
    cfg$output %||% list()
  )
  cfg
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' 建立输出目录
ensure_dirs <- function(cfg) {
  for (d in c(cfg$output$results_dir, cfg$output$data_dir)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(TRUE)
}

# ---- JSON / 状态 -----------------------------------------------------------

#' 写 JSON（原子写：先写临时文件再改名，避免半截文件）
write_json <- function(path, obj) {
  tmp <- paste0(path, ".tmp")
  writeLines(jsonlite::toJSON(obj, auto_unbox = TRUE, pretty = TRUE, null = "null"), tmp)
  if (file.exists(path)) unlink(path)
  file.rename(tmp, path)
  invisible(path)
}

#' 记录单个步骤的执行结果到 results/state.json
record_step <- function(cfg, id, status, seconds = NA_real_, message = "", required = TRUE) {
  path <- file.path(cfg$output$results_dir, "state.json")
  state <- if (file.exists(path)) {
    tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE), error = function(e) list(steps = list()))
  } else {
    list(steps = list())
  }
  state$steps <- Filter(function(s) !identical(s$id, id), state$steps %||% list())
  state$steps[[length(state$steps) + 1L]] <- list(
    id = id, status = status, required = required,
    seconds = round(seconds, 1), message = message
  )
  failed <- Filter(function(s) identical(s$status, "failed") && isTRUE(s$required), state$steps)
  state$required_failed <- vapply(failed, function(s) s$id, character(1))
  state$optional_failed <- vapply(
    Filter(function(s) identical(s$status, "failed") && !isTRUE(s$required), state$steps),
    function(s) s$id, character(1)
  )
  state$updated <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  write_json(path, state)
  invisible(state)
}

# ---- 统一调色板 -------------------------------------------------------------
#
# **语义固定，所有图共用同一套。** 之前每张图各自挑颜色，同一个含义在不同图里
# 是不同颜色：PCA 把第一组画成蓝色，火山图却把上调（tumor 高）画成红色；
# p 值直方图的参考线是第三种红（firebrick）；dotplot 用 viridis plasma，
# 与红蓝体系毫无关系。读者要重新学一遍每张图的配色。
#
# **这套色值是算出来的，不是挑出来的。** 判据取自 better-colors skill：
#   * 色相相差 15° 以内视为**同一个颜色** —— 承载不同含义的颜色必须拉开 15° 以上
#   * 颜色从不是唯一的语义载体（见各图的 shape / 位置编码）
#   * 报告任何对比度之前先测量，不要估
#
# 实测（`_palette_audit.py`）暴露了三个凭眼睛看不出的问题：
#   1. 旧 p 值直方图用的 `firebrick` 与 up 红**只差 0.4°** —— 按上述判据就是同一个颜色
#   2. 旧分类色板**不是色盲安全的**：protanopia 下 up 红与 cat brown 距离 0.002
#      （几乎完全重合），deuteranopia 下 up 红与 cat green 距离 0.043
#   3. magma 序列色的中段 #F1605D 色相 24.3°，**距 up 红仅 2.9°** ——
#      序列色的中段就是"上调红"，深色点会被误读成上调
#
# 约定：
#   * **方向 / 条件**：红 = 上调 / tumor，蓝 = 下调 / normal
#   * **序列（显著性、强度）**：紫色单色相 ramp。刻意避开红与蓝 ——
#     红色已经被"方向"占了。magma 的色相跨度 171.7°，本 ramp 只有 9.8°，
#     亮度步长也更均匀（最大 0.122 vs 0.157），彩度峰值落在中段。
#   * **分类（模块、多组）**：Okabe-Ito 的色盲安全子集，
#     剔除 orange（距 nominal 仅 10°）与 blue（距 down 仅 13°）。
#     三种色盲下两两最小 OKLab 距离 0.076，阈值 0.05。
PAL <- list(
  up       = "#C1443C",   # 上调 / tumor 一侧      h=27.2
  down     = "#2E5FA3",   # 下调 / normal 一侧     h=257.2
  ns       = "#BFBFBF",   # 不显著（近中性灰，不承载色相语义）
  nominal  = "#C07A1E",   # 名义显著但未过 FDR（h=67.1，距 up 红 40°，白底 3.5:1）
  mid      = "#F5F5F5",   # 发散色中点
  ink      = "#333333",   # 文字 / 参考线（白底 12.6:1）
  grid     = "#E5E5E5"
)

#' 发散色板：低 = 蓝 → 中 = 近白 → 高 = 红
#' 用于 Z-score 热图、logFC 类连续量
pal_diverging <- function(n = 100) {
  grDevices::colorRampPalette(c(PAL$down, PAL$mid, PAL$up))(n)
}

#' 序列色板：紫色单色相 ramp，用于显著性 / 强度
#'
#' **刻意不用 magma。** 实测 magma 中段 #F1605D 的色相是 24.3°，
#' 而 up 红是 27.2° —— 只差 2.9°，按 15° 判据就是同一个颜色。
#' 也就是说 magma 中段的点会被读成"上调"。本 ramp 色相跨度 9.8°，
#' 亮度单调递减且步长均匀，彩度峰值在中段，距 up 红 70°、距 down 蓝 50°。
pal_sequential <- function(n = 256) {
  grDevices::colorRampPalette(
    c("#F4ECF7", "#DDC8E8", "#C09FD6", "#9E74BE", "#7B4E9E", "#5A3279", "#3B1F52"))(n)
}

#' 分类色板（模块 / 多组）：经计算验证的色盲安全四色
#'
#' 这四个不是挑的，是搜出来的（`_palette_search.py`）：在候选池里贪心挑选，
#' 同时满足三条约束 ——
#'   * 与 up 红 / down 蓝 / nominal 橙 / 序列紫 的色相都拉开 15° 以上
#'   * 与它们、以及彼此之间，在 protanopia / deuteranopia / tritanopia 下
#'     OKLab 距离都 >= 0.05
#'   * 白底对比度 >= 2.0（实心圆点仍清晰可见的底线）
#'
#' 被剔除的典型：`#E69F00`（距 nominal 仅 10°）、`#0072B2`（距 down 仅 13°）、
#' `#007A59` `#14606D` `#4A4E8C`（都撞 down）、`#A65A2E`（撞 up）、
#' `#F0E442`（白底仅 1.32:1，小圆点看不见）。
#'
#' **超过 4 个模块时降级**：在四色之间插值。插出来的颜色**没有**经过上面的
#' 验证，可能色盲下不可分。实测 GSE64790 的 PPI 是 4 个模块，正好用满。
pal_categorical <- function(k) {
  base <- c("#009E73", "#56B4E9", "#C9B800", "#6E7B8B")
  if (k <= length(base)) return(base[seq_len(k)])
  grDevices::colorRampPalette(base)(k)
}

#' 按对比的两端给条件上色
#'
#' 两端时 tumor/上调 = 红、normal/下调 = 蓝；超过两类时前两个仍是方向色，
#' 其余从分类色板补（>2 组的回退路径，本数据集用不到）。
#' **这是修掉"PCA 蓝、火山红"那个不一致的地方** —— 条件色必须由
#' `cfg$contrast` 决定，不能由分组出现的先后顺序决定。
pal_condition <- function(arms) {
  arms <- as.character(arms)
  if (length(arms) == 2L) return(stats::setNames(c(PAL$up, PAL$down), arms))
  stats::setNames(c(PAL$up, PAL$down, pal_categorical(max(0L, length(arms) - 2L)))[seq_along(arms)],
                  arms)
}

#' 分组配色：按 `cfg$contrast` 的两端排序后再上色
#'
#' 放在 common.R 而不是 02 里，因为 04（热图注释条）也要用 ——
#' 定义在 02 会让 `Rscript scripts/04_heatmap_enrichment.R` 单独跑时找不到函数。
group_palette <- function(groups, cfg = NULL) {
  lv <- unique(as.character(groups))
  arms <- if (!is.null(cfg)) as.character(cfg$contrast) else character(0)
  ordered <- c(intersect(arms, lv), sort(setdiff(lv, arms)))
  stats::setNames(unname(pal_condition(ordered)), ordered)
}

# ggplot2 比例尺包装，保证同一语义在每张图上写法一致
scale_colour_condition <- function(arms, name = NULL, ...) {
  ggplot2::scale_colour_manual(values = pal_condition(arms), name = name, ...)
}
scale_fill_condition <- function(arms, name = NULL, ...) {
  ggplot2::scale_fill_manual(values = pal_condition(arms), name = name, ...)
}
scale_colour_seq <- function(name, ...) {
  ggplot2::scale_colour_gradientn(colours = pal_sequential(), name = name, ...)
}
scale_fill_seq <- function(name, ...) {
  ggplot2::scale_fill_gradientn(colours = pal_sequential(), name = name, ...)
}
scale_fill_div <- function(name, ...) {
  ggplot2::scale_fill_gradientn(colours = pal_diverging(), name = name, ...)
}

#' 所有图共用的主题，保证字号、网格、留白一致
theme_paper <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = PAL$grid, linewidth = 0.25),
      panel.border     = ggplot2::element_rect(colour = PAL$grid, linewidth = 0.4),
      plot.title       = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle    = ggplot2::element_text(colour = "#666666", size = base_size - 1.5),
      legend.key.size  = ggplot2::unit(0.9, "lines")
    )
}

# ---- 绘图 ------------------------------------------------------------------

#' 把绘图表达式同时写进 PDF 和 PNG，保证设备一定关闭
#'
#' **为什么要同时出 PNG：** 产物是打包成 artifact zip 下载的，PDF 在 zip 里
#' 不能直接预览 —— 拿到 artifact 的人得先解压再找 PDF 阅读器。PNG 可以直接看。
#' PDF 保留是因为它是矢量图，放大不失真；PNG 是为了能一眼看到。
#'
#' **实现要点：** 绘图代码要跑两遍（一次 PDF、一次 PNG），所以必须用
#' `substitute()` 抓住**未求值**的表达式再 `eval()` 两次。
#' 不能写成 `force(expr); force(expr)` —— R 的 promise 有记忆，
#' 第二次 force 直接返回缓存值，**PNG 设备会开了又关、什么都不画**，得到一张空白图。
#'
#' 每次绘图都用一个独立的 `render()` 开关设备，`on.exit` 才精确对应这一次开设备；
#' 若在 `save_pdf` 主体里 `on.exit(add = TRUE)` 两次，退出时会多关一次设备，
#' 可能把调用方的设备一起关掉。
#'
#' PNG 走 `ragg`（若装了）或 `png(type="cairo")`，两者抗锯齿都更好；
#' 都没有时退回默认设备，仍然出图，不因为画质问题中断流程。
#'
#' @param path  PDF 输出路径；同名 `.png` 会写在旁边
#' @param expr  绘图表达式，会在调用者的环境里求值两次
save_pdf <- function(path, expr, width = 8, height = 6, dpi = 150) {
  code <- substitute(expr)
  env  <- parent.frame()

  render <- function(open_dev) {
    open_dev()
    on.exit(grDevices::dev.off(), add = TRUE)
    eval(code, envir = env)
  }

  render(function() grDevices::pdf(path, width = width, height = height))

  png_path <- sub("\\.pdf$", ".png", path)
  tryCatch(
    render(function() {
      if (requireNamespace("ragg", quietly = TRUE)) {
        ragg::agg_png(png_path, width = width, height = height, units = "in", res = dpi)
      } else if (capabilities("cairo")) {
        grDevices::png(png_path, width = width, height = height, units = "in",
                       res = dpi, type = "cairo")
      } else {
        grDevices::png(png_path, width = width * dpi, height = height * dpi, res = dpi)
      }
    }),
    error = function(e) {
      # 出图失败不能中断分析：PDF 已经拿到了，PNG 只是方便预览
      log_warn(sprintf("PNG 输出失败（PDF 已生成）: %s", conditionMessage(e)))
    }
  )
  invisible(path)
}

# ---- GEO SOFT 抓取（纯 base R，不需要 Bioconductor）------------------------

#' 从 GEO SOFT 文本接口抓取元数据
#'
#' 用 `form=text&view=brief` 拿到的是纯文本，比下载完整 SOFT family 小几个数量级，
#' 而且不需要 GEOquery —— 这样前置校验可以在装任何 Bioconductor 包之前跑完。
#'
#' @param acc GSE 或 GSM 编号
#' @param targ "self" 取 series 级，"gsm" 取全部样本级
#' @return 字符向量，每行一条 SOFT 记录
fetch_geo_soft <- function(acc, targ = "self", view = "brief", retries = 3L) {
  url <- sprintf(
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=%s&form=text&view=%s",
    acc, targ, view
  )
  last_err <- NULL
  for (attempt in seq_len(retries)) {
    con <- NULL
    res <- tryCatch({
      con <- url(url, open = "rt")
      lines <- readLines(con, warn = FALSE)
      if (length(lines) == 0L) stop("返回为空")
      lines
    }, error = function(e) {
      last_err <<- conditionMessage(e)
      NULL
    }, finally = {
      if (!is.null(con)) try(close(con), silent = TRUE)
    })
    if (!is.null(res)) return(res)
    log_warn(sprintf("GEO 抓取失败 (第 %d/%d 次): %s", attempt, retries, last_err))
    Sys.sleep(2 * attempt)
  }
  stop(sprintf("无法从 GEO 获取 %s 的元数据: %s", acc, last_err))
}

#' 从 SOFT 行中提取某个键的全部取值
#' @param lines SOFT 行
#' @param key 不含前导 "!" 的键名，例如 "Sample_characteristics_ch1"
soft_values <- function(lines, key) {
  pattern <- paste0("^!", key, " = ")
  hits <- grep(pattern, lines, value = TRUE)
  sub(pattern, "", hits)
}

#' 提取 SOFT 中第一个匹配键的值
soft_value <- function(lines, key, default = NA_character_) {
  v <- soft_values(lines, key)
  if (length(v) == 0L) default else v[[1L]]
}

# ---- 基因/统计小工具 -------------------------------------------------------

#' 按行 Z-score（用于热图）
row_zscore <- function(m) {
  t(scale(t(as.matrix(m))))
}

#' 为下游（富集、PPI）挑选差异基因
#'
#' **为什么需要降级路径：** spec 要求样本数 < 10，而在这个量级上，
#' 对全基因组（约 1.6 万个基因）做 BH 校正几乎不可能有任何基因通过 ——
#' GSE64790（n=6，3 对配对）实测：1,456 个基因 raw P < 0.05 且 |log2FC| > 1，
#' 但最小的 adj.P 也有 0.394。这是样本量本身的限制，不是分析错误。
#'
#' 所以：FDR 显著基因够用时用 FDR；不够时退回「raw P 排序前 N 个（仍要求
#' |log2FC| > 阈值）」。**降级必须被标注**，mode 会写进
#' enrichment_status.json / ppi_status.json，结论里不得把它当成显著差异基因。
#'
#' **返回的 up / down 是分开的。** ORA 本身是方向无关的：把上调和下调基因
#' 混在一个列表里跑，"富集到 X 通路"到底是被上调基因驱动还是被下调基因驱动
#' 就分不清了（K-Dense `pathway-enrichment`：*"ORA is direction-agnostic unless
#' you split up/down lists; GSEA NES sign gives direction"*）。
#' 对肿瘤 vs 正常组织尤其致命 —— 上调的是增殖，下调的是基质/脂肪，
#' 混起来会得到"两条方向相反的通路同时富集"这种无法解释的结果。
select_degs <- function(deg, cfg, min_genes = 5L) {
  padj <- cfg$thresholds$adj_p
  lfc  <- cfg$thresholds$log2fc
  # 样本数只用于把降级原因写清楚；从 03 步的摘要里取，取不到就算了
  sp <- file.path(cfg$output$results_dir, "deg_summary.json")
  n_s <- if (file.exists(sp)) {
    tryCatch(jsonlite::fromJSON(sp)$n_samples, error = function(e) NA_integer_)
  } else NA_integer_
  if (is.null(n_s) || length(n_s) == 0L) n_s <- NA_integer_

  # 按 logFC 符号拆方向，供 ORA 分方向跑
  by_direction <- function(tab) {
    if (is.null(tab) || nrow(tab) == 0L) return(list(up = character(0), down = character(0)))
    list(up   = unique(tab$gene[!is.na(tab$logFC) & tab$logFC > 0]),
         down = unique(tab$gene[!is.na(tab$logFC) & tab$logFC < 0]))
  }

  sig <- deg[!is.na(deg$adj.P.Val) & deg$adj.P.Val < padj & abs(deg$logFC) > lfc, , drop = FALSE]
  if (nrow(sig) >= min_genes) {
    d <- by_direction(sig)
    return(list(
      genes = unique(sig$gene), mode = "fdr", n = nrow(sig), table = sig,
      up = d$up, down = d$down,
      reason = sprintf("adj.P < %g 且 |log2FC| > %g，共 %d 个（上调 %d / 下调 %d）",
                       padj, lfc, nrow(sig), length(d$up), length(d$down))
    ))
  }

  top_n <- cfg$analysis$ranked_fallback_genes
  if (is.null(top_n)) top_n <- 500L
  cand <- deg[!is.na(deg$P.Value) & abs(deg$logFC) > lfc, , drop = FALSE]
  cand <- cand[order(cand$P.Value), , drop = FALSE]
  cand <- utils::head(cand, top_n)
  d <- by_direction(cand)
  list(
    genes = unique(cand$gene), mode = "ranked_fallback", n = nrow(cand), table = cand,
    up = d$up, down = d$down,
    reason = sprintf(paste0(
      "FDR 显著基因仅 %d 个（需 >= %d）。样本数 %s 下对 %d 个基因做 BH 校正过严，",
      "最小的 adj.P 为 %.3f。退回按 raw P 排序、|log2FC| > %g 的前 %d 个基因",
      "（上调 %d / 下调 %d）。",
      "**这是假设生成，不是显著差异基因清单。**"),
      nrow(sig), min_genes, if (is.na(n_s)) "很少" else as.character(n_s),
      nrow(deg), if (nrow(deg)) min(deg$adj.P.Val, na.rm = TRUE) else NA_real_,
      lfc, nrow(cand), length(d$up), length(d$down))
  )
}

#' 按基因重叠把冗余的条目折叠成代表条目
#'
#' GO 会返回大量近义条目 —— "mitotic cell cycle phase transition"、
#' "regulation of mitotic cell cycle phase transition"、
#' "positive regulation of mitotic cell cycle phase transition"……
#' 列 44 条这种东西不是 44 个发现，是 1 个发现重复了 44 次。
#' K-Dense `pathway-enrichment`：*"GO especially returns many near-duplicate
#' terms. Collapse with an enrichment map (term–term similarity), leading-edge
#' overlap, or parent terms, and report representative terms."*
#'
#' 这里用**基因重叠 Jaccard 的单链接聚类**实现（不引入 GOSemSim 这类重依赖）：
#' 两个条目共享的基因占并集的比例 >= 阈值就归为一类，每类保留 adj.P 最小的
#' 那个作代表。返回的 representative 列标出每条属于哪一类、该类有几个成员。
#'
#' @param tab 含 `ID` / `Description` / `p.adjust` 与基因列的富集结果表
#' @param threshold Jaccard 阈值；>= 1 时关闭去冗余
#' @param gene_col 基因列的列名。ORA 是 `geneID`，GSEA 是 `core_enrichment`
#'   （两者都是 `/` 分隔的基因串，但列名不同，所以必须显式传）
reduce_terms_by_overlap <- function(tab, threshold = 0.5, gene_col = "geneID") {
  if (is.null(tab) || nrow(tab) == 0L) return(tab)
  tab$representative <- NA_character_
  tab$cluster_size <- 1L
  if (!all(c("ID", gene_col) %in% colnames(tab))) return(tab)
  if (is.na(threshold) || threshold >= 1) return(tab)

  sets <- strsplit(as.character(tab[[gene_col]]), "/", fixed = TRUE)
  names(sets) <- as.character(tab$ID)
  n <- nrow(tab)

  # 单链接聚类：并查集
  parent <- seq_len(n)
  find <- function(i) { while (parent[i] != i) { parent[i] <<- parent[parent[i]]; i <- parent[i] }; i }
  union <- function(i, j) { ri <- find(i); rj <- find(j); if (ri != rj) parent[ri] <<- rj }

  for (i in seq_len(n - 1L)) {
    for (j in seq(i + 1L, n)) {
      inter <- length(intersect(sets[[i]], sets[[j]]))
      if (inter == 0L) next
      jac <- inter / (length(sets[[i]]) + length(sets[[j]]) - inter)
      if (jac >= threshold) union(i, j)
    }
  }

  roots <- vapply(seq_len(n), find, integer(1))
  for (r in unique(roots)) {
    idx <- which(roots == r)
    # 代表 = 该类里 adj.P 最小的条目
    best <- idx[which.min(tab$p.adjust[idx])]
    tab$representative[idx] <- as.character(tab$ID[best])
    tab$cluster_size[idx] <- length(idx)
  }
  # 按类大小与显著性排序，代表条目排在自己那一类的首位
  tab <- tab[order(tab$representative, tab$p.adjust), , drop = FALSE]
  rownames(tab) <- NULL
  tab
}

#' 把富集结果里的基因列表列排序归一化
#'
#' `enrichGO` / `gseGO` 写出的 `geneID` / `core_enrichment` 是 `/` 分隔的基因串，
#' 其**内部顺序不稳定** —— 同样的基因、同样的 p 值，两轮运行的行顺序可以不同
#' （实测 234 行里 230 行的 `geneID` 排列不同，而 p 值、计数完全一致）。
#' 这对科学结论没有影响，但让产物无法逐字节比对，
#' "两轮结果是否一致" 这种验证就做不了。
#'
#' 按字母排序即可消除这一差异。GSEA 的 `core_enrichment` 是 leading-edge 子集，
#' 顺序本来也没有语义，排序是安全的。
normalise_gene_lists <- function(tab, cols = c("geneID", "core_enrichment")) {
  if (is.null(tab) || nrow(tab) == 0L) return(tab)
  for (cl in intersect(cols, colnames(tab))) {
    tab[[cl]] <- vapply(strsplit(as.character(tab[[cl]]), "/", fixed = TRUE),
                        function(g) paste(sort(g), collapse = "/"), character(1))
  }
  tab
}

#' 安全地调用一个可选包；缺失时返回 NULL 而不是报错
require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    log_warn(sprintf("可选包 %s 未安装，跳过相关分析", pkg))
    return(FALSE)
  }
  TRUE
}
