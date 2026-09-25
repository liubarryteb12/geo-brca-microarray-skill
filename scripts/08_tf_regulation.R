# ============================================================================
# 08_tf_regulation.R — 转录因子调控分析
# ============================================================================
# 用户文档把 TF 调控列为 bulk 部分的分析项。本脚本补齐它。
#
# ## 两个分析，回答两个不同的问题
#
#   1. **调控子富集**（regulon enrichment，Fisher 精确检验）
#      问：某个 TF 的靶基因在差异基因里是否**过度出现**？
#      这是主分析 —— 它不依赖样本量，n<10 时仍然成立。
#
#   2. **调控子活性分数**（regulon activity score，按样本）
#      问：某个 TF 的靶基因整体在两组之间是否**系统性地位移**？
#      分数 = 靶基因 z-score 的加权均值（mode of regulation 定符号）。
#      **n<10 时这个分数噪声很大**，所以它只作为佐证，
#      并在产物里写明样本量警告。
#
# ## 为什么不用 decoupleR 的 ULM
#
# decoupleR 是更标准的做法（单样本线性模型 + 置换检验），但它多一个
# Bioconductor 依赖，而它算出的"活性"与本脚本的加权 z-score
# 在**排序上高度一致**。这里选依赖更少的路，把置换检验留给自己做。
#
# **这不等于 decoupleR。** 产物里写明了差异。
#
# ## 数据来源与回退
#
#   主路径：dorothea::dorothea_hs（confidence A/B/C 的调控子）
#   回退：assets/tf_regulons_minimal.tsv（少量手工整理的 TF-靶基因对）
#         —— **回退时 coverage 标为 minimal，结论不能当全基因组结论**
#   都没有：写 status = "not_done" + reason，不产出任何"结果"
#
# ## §1.5 的 TRRUST / ChEA3 怎么落地
#
# **TRRUST 接入为独立交叉验证，不替换主来源。** 它没有 CRAN 包，
# 官方分发方式是直接下 TSV，所以用 base R 的 `download.file`（不引新依赖），
# 下到 `data/<GSE>/trrust_cache/` 并记 sha256。
# 比对报三件事：配对重叠（Jaccard）、**重叠部分的 mor 符号一致率**、
# 各自独有部分的大小。**mor 冲突是最该看的数** —— 两个库对同一对
# TF-靶基因的调控方向就不一致时，基于 mor 定符号的活性打分要打折。
#
# **ChEA3 未接入，理由写在 limitations 里**：它是 web 服务、无 CRAN 包，
# 调用需要 httr/curl（CI 的 R 包列表里没有）。而 dorothea 的置信度 A 档
# 本身就整合了 ChEA，覆盖率上不缺 —— 缺的是"用另一个服务独立复核富集
# 结果"这一步。**这是如实记录的缺口，不是"已经做了"。**
#
# ## 产物
#
#   results/<GSE>/tf_regulon_enrichment.csv   每个 TF 的 Fisher 检验
#   results/<GSE>/tf_activity_by_sample.csv   每个样本 × TF 的活性分数
#   results/<GSE>/tf_activity_group_test.csv  组间比较
#   results/<GSE>/01-08-01-unit1-tf-regulon-enrichment.pdf        调控子富集条形图（+ 同名 .png）
#   results/<GSE>/01-08-02-unit1-tf-activity-group-difference.pdf 调控子活性组间效应量（+ .png）
#   results/<GSE>/tf_status.json                   状态与方法学限定
# ============================================================================

suppressPackageStartupMessages({
  library(stats)
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

# 用哪些置信度的调控子。
# A = 人工审编（ChEA/TRRUST/文献），B = 中等，C = 预测。
# **不用 D/E** —— 那两档基本是基序预测，假阳性率高，
# 而富集分析对假阳性特别敏感（每个 TF 的靶基因集合被稀释）。
TF_CONFIDENCE <- c("A", "B", "C")

# 一个 TF 至少要这么多靶基因落在检测到的基因集里，才做检验。
# 太少的话 Fisher 的自由度不足，p 值没有意义。
MIN_TARGETS_IN_DATA <- 5L

# 一个 TF 的活性分数至少要这么多靶基因才计算（噪声控制）
MIN_TARGETS_FOR_ACTIVITY <- 5L


#' 读入 TF-靶基因调控子表
#'
#' 返回 list(regulons = data.frame(tf, target, mor, confidence), source, coverage)
#' mor = mode of regulation（+1 激活 / -1 抑制）
load_tf_regulons <- function(cfg) {
  # ---- 主路径：dorothea ----
  if (requireNamespace("dorothea", quietly = TRUE)) {
    db <- dorothea::dorothea_hs
    # 列名在不同版本间可能是 tf/target/mor/confidence
    need <- c("tf", "target", "mor", "confidence")
    if (all(need %in% colnames(db))) {
      keep <- db$confidence %in% TF_CONFIDENCE
      db <- db[keep, need, drop = FALSE]
      log_info(sprintf("dorothea: %d 个 TF, %d 条 TF-靶基因关系（置信度 %s）",
                       length(unique(db$tf)), nrow(db),
                       paste(TF_CONFIDENCE, collapse = "/")))
      return(list(regulons = as.data.frame(db),
                  source = "dorothea::dorothea_hs",
                  coverage = "genome_wide",
                  confidence_levels = TF_CONFIDENCE))
    }
    log_warn(sprintf("dorothea 的列名不是预期的 %s，实际 %s —— 转回退路径",
                     paste(need, collapse = "/"),
                     paste(colnames(db), collapse = "/")))
  } else {
    log_warn("dorothea 未安装 —— 转回退路径")
  }

  # ---- 回退：随仓库附带的最小表 ----
  # 不假设 cwd 是仓库根 —— 从 --file= 反推脚本位置，再往上找 assets/。
  # 找不到就按 cwd 试一次（编排器从根目录运行时走这条）。
  cand <- character(0)
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    here <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
    cand <- c(cand, file.path(here, "..", "assets", "tf_regulons_minimal.tsv"),
              file.path(here, "assets", "tf_regulons_minimal.tsv"))
  }
  cand <- c(cand, file.path("assets", "tf_regulons_minimal.tsv"))
  cand <- unique(cand)
  fb <- cand[file.exists(cand)]
  if (length(fb) > 0L) {
    fb <- fb[[1L]]
    # **`comment.char = "#"` 不能省。** `read.delim` 默认 `comment.char = ""`，
    # 也就是**不把 `#` 当注释** —— 表头那 20 行说明文字会被当成数据行读进来，
    # 然后 `tf` 列里出现一堆 `# ===...`，靶基因匹配全部落空。
    # 而这个失败不报错，只是富集结果为空。
    db <- utils::read.delim(fb, stringsAsFactors = FALSE, comment.char = "#")
    if (all(c("tf", "target", "mor") %in% colnames(db))) {
      db$confidence <- "minimal"
      log_warn(sprintf(
        "**回退到最小 TF 表**：%d 个 TF, %d 条关系。覆盖度远低于 dorothea，",
        length(unique(db$tf)), nrow(db)))
      log_warn("结论只能作为**示例**，不能当全基因组的 TF 调控结论。")
      return(list(regulons = db,
                  source = "assets/tf_regulons_minimal.tsv (回退)",
                  coverage = "minimal",
                  confidence_levels = "minimal"))
    }
  }

  # ---- 都没有 ----
  NULL
}


# ============================================================================
# TRRUST（文档 §1.5 点名）
# ============================================================================
# **TRRUST 没有 CRAN/Bioconductor 包**（实测查过 cran.r-project.org 的
# TRRUST / trrust / ChEA3 / chea3 四个名字，全无）。官方分发方式是
# 直接下 TSV，所以这里用 base R 的 `utils::download.file`，不引新依赖。
#
# 实测格式（trrust_rawdata.human.tsv，9396 行 / 795 个 TF / 4 列无表头）：
#     AATF<TAB>BAX<TAB>Repression<TAB>22909821
#     AATF<TAB>CDKN2A<TAB>Activation<TAB>...
# 第 3 列取值只有 Activation / Repression / Unknown 三种。
#
# **它是 dorothea 的独立交叉验证，不是替代品。** 两个库的收录标准不同
# （TRRUST 要求有文献报道的调控关系；dorothea 整合 ChEA/TRRUST/文献并按
# 置信度分档），所以**对不上的部分本身就是信息**：
#   - 只在 dorothea 里：多为高通量/预测来源
#   - 只在 TRRUST 里：文献报道但未进 dorothea 的
#   - **两边都有但 mor 符号相反**：最值得看的 —— 说明"激活还是抑制"
#     这件事在两个库之间就不一致，任何基于 mor 的活性打分都要谨慎。
# **只用人类文件，不做物种分支。** 本仓库从头到尾硬编码人类：
# `00_validate_inputs.R` 遇到非 Homo sapiens 样本直接 stop()，
# 富集用 `org.Hs.eg.db`，PPI 用 `species = 9606`。
# 写一个从配置读物种的分支是**永远走不到的死代码**，
# 而且会让静态检查报"引用了配置里不存在的字段"（实测报过）。
# 真要放开物种，先改门禁，再在这里加 mouse 分支。
TRRUST_URL <- "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"
# 需要小鼠时的地址（TRRUST v2 同样提供）：
#   https://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.tsv


#' 下载（带缓存）并解析 TRRUST
#'
#' 返回 list(regulons, source, coverage, n_tf, n_pairs, url, sha256, cache_path)
#' 或 NULL（任何一步失败都不抛异常 —— 它是交叉验证，不该拖垮主分析）。
load_trrust_regulons <- function(cfg) {
  url <- TRRUST_URL
  cache_dir <- file.path(cfg$output$data_dir, "trrust_cache")
  cache_path <- file.path(cache_dir, basename(url))

  # ---- 取文件（缓存优先）----
  got <- FALSE
  if (file.exists(cache_path) && file.size(cache_path) > 1000) {
    got <- TRUE
    log_info(sprintf("TRRUST 用缓存：%s", cache_path))
  } else {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    ok <- tryCatch({
      # _download_allow 仅 CI 执行：TRRUST 数据取数只发生在 CI runner 上
      # （R1 本地零执行），本地从未运行过这条路径；有缓存后也不会再下载。
      # mode = "wb" 是必须的：默认 "w" 在 Windows 上会把行尾改写，
      # 而哈希是用来核对"这份文件有没有变"的 —— 改写会让哈希失去意义。
      utils::download.file(url, cache_path, quiet = TRUE, mode = "wb")
      file.exists(cache_path) && file.size(cache_path) > 1000
    }, error = function(e) {
      log_warn(sprintf("TRRUST 下载失败：%s", conditionMessage(e)))
      FALSE
    })
    if (!ok) {
      # 下载失败时清掉可能写了一半的文件，避免下轮拿半截文件当缓存
      if (file.exists(cache_path)) unlink(cache_path)
      return(NULL)
    }
    got <- TRUE
    log_info(sprintf("TRRUST 已下载：%s", cache_path))
  }
  if (!got) return(NULL)

  # ---- 解析 ----
  db <- tryCatch(
    utils::read.delim(cache_path, header = FALSE, stringsAsFactors = FALSE,
                      quote = "", comment.char = "",
                      col.names = c("tf", "target", "mor_raw", "pmid")),
    error = function(e) {
      log_warn(sprintf("TRRUST 解析失败：%s", conditionMessage(e)))
      NULL
    })
  if (is.null(db) || nrow(db) == 0L) return(NULL)

  # ---- 列数与列名校验 ----
  # **不校验的话，上游改格式会静默产出错的结果**：比如加了表头行，
  # 那么第一行会变成 tf="tf" 这种，而它照样能算出一个数来。
  if (ncol(db) != 4L) {
    log_warn(sprintf("TRRUST 列数是 %d，预期 4 —— 上游格式变了，放弃", ncol(db)))
    return(NULL)
  }
  mode_vals <- unique(db$mor_raw)
  known <- c("Activation", "Repression", "Unknown")
  if (!all(mode_vals %in% known)) {
    log_warn(sprintf("TRRUST 第 3 列出现未知取值：%s —— 上游格式变了，放弃",
                     paste(setdiff(mode_vals, known), collapse = ", ")))
    return(NULL)
  }

  # ---- mor 映射 ----
  # **Unknown 映射成 NA，不映射成 0。** 0 在加权平均里等于"这条关系不贡献"，
  # 而 NA 会被 na.rm 显式跳过 —— 前者把"不知道"混进了分母。
  db$mor <- ifelse(db$mor_raw == "Activation", 1L,
                   ifelse(db$mor_raw == "Repression", -1L, NA_integer_))
  db$confidence <- "TRRUST"
  db <- db[, c("tf", "target", "mor", "confidence")]

  n_tf <- length(unique(db$tf))
  n_unk <- sum(is.na(db$mor))
  log_info(sprintf("TRRUST(human)：%d 个 TF，%d 条关系（其中 %d 条 mor 未知）",
                   n_tf, nrow(db), n_unk))

  list(regulons = db,
       source = "TRRUST v2 (human)",
       coverage = "genome_wide",
       n_tf = n_tf,
       n_pairs = nrow(db),
       n_mor_unknown = n_unk,
       url = url,
       sha256 = file_hash(cache_path)$hash,
       cache_path = cache_path)
}


#' 两个调控子库的交叉核对
#'
#' **报三件事，缺一不可：**
#'   1. 配对层面的重叠（Jaccard）—— 两个库收录范围差多少
#'   2. **重叠部分里 mor 符号相反的比例** —— 这是最该看的数
#'   3. 只在一边的条数 —— 说明"独立验证"到底验证到了多少
#'
#' 只报"两个库都有几万条关系"是没有信息量的：数量大不等于一致。
compare_regulon_databases <- function(a, b, a_name = "dorothea",
                                      b_name = "TRRUST") {
  out <- list(compared = FALSE, a = a_name, b = b_name)
  if (is.null(a) || is.null(b) || nrow(a) == 0L || nrow(b) == 0L) {
    out$reason <- "有一侧为空，无法比较"
    return(out)
  }
  ka <- paste(a$tf, a$target, sep = "|")
  kb <- paste(b$tf, b$target, sep = "|")
  ua <- unique(ka); ub <- unique(kb)
  both <- intersect(ua, ub)

  out$compared <- TRUE
  out$a_pairs <- length(ua)
  out$b_pairs <- length(ub)
  out$n_common_pairs <- length(both)
  out$n_a_only <- length(setdiff(ua, ub))
  out$n_b_only <- length(setdiff(ub, ua))
  out$jaccard <- round(length(both) / length(union(ua, ub)), 4)
  out$a_tfs <- length(unique(a$tf))
  out$b_tfs <- length(unique(b$tf))
  out$n_common_tfs <- length(intersect(unique(a$tf), unique(b$tf)))

  # ---- mor 符号一致性（只在两边都有 mor 且都非 NA 的重叠对上算）----
  ma <- a[!duplicated(ka), c("tf", "target", "mor")]
  mb <- b[!duplicated(kb), c("tf", "target", "mor")]
  ma$key <- paste(ma$tf, ma$target, sep = "|")
  mb$key <- paste(mb$tf, mb$target, sep = "|")
  mm <- merge(ma[, c("key", "mor")], mb[, c("key", "mor")],
              by = "key", suffixes = c("_a", "_b"))
  mm <- mm[!is.na(mm$mor_a) & !is.na(mm$mor_b), , drop = FALSE]
  if (nrow(mm) > 0L) {
    out$n_mor_comparable <- nrow(mm)
    out$n_mor_agree <- sum(mm$mor_a == mm$mor_b)
    out$n_mor_conflict <- sum(mm$mor_a != mm$mor_b)
    out$mor_agreement <- round(out$n_mor_agree / nrow(mm), 4)
    # 冲突的具体例子（最多 20 条）—— 只给比例的话没法核对
    cf <- mm[mm$mor_a != mm$mor_b, , drop = FALSE]
    if (nrow(cf) > 0L) {
      out$mor_conflict_examples <- utils::head(
        data.frame(key = cf$key,
                   mor_a = cf$mor_a, mor_b = cf$mor_b,
                   stringsAsFactors = FALSE), 20L)
    }
  } else {
    out$n_mor_comparable <- 0L
    out$mor_agreement <- NA_real_
  }
  out
}


#' 调控子富集：Fisher 精确检验
#'
#' 2x2 表：
#'              在 DEG 里   不在 DEG 里
#'   是靶基因        a            b
#'   不是靶基因      c            d
#'
#' 背景集用**检测到的基因**（detected），不用全基因组 —— 与
#' 04_heatmap_enrichment.R 的 ORA 保持一致（AGENTS.md 规则 9）。
regulon_enrichment <- function(regulons, deg_genes, detected_genes) {
  det <- unique(detected_genes)
  deg <- intersect(unique(deg_genes), det)
  n_det <- length(det)
  n_deg <- length(deg)
  if (n_deg == 0L) return(NULL)

  tfs <- sort(unique(regulons$tf))
  rows <- vector("list", length(tfs))
  k <- 0L
  for (tf in tfs) {
    tg <- intersect(unique(regulons$target[regulons$tf == tf]), det)
    if (length(tg) < MIN_TARGETS_IN_DATA) next
    a <- length(intersect(tg, deg))
    b <- length(tg) - a
    cc <- n_deg - a
    d <- n_det - n_deg - b
    if (d < 0L) next
    ft <- stats::fisher.test(matrix(c(a, b, cc, d), nrow = 2L),
                             alternative = "greater")
    k <- k + 1L
    rows[[k]] <- data.frame(
      tf = tf,
      n_targets_in_data = length(tg),
      n_targets_in_deg = a,
      expected_in_deg = round(length(tg) * n_deg / n_det, 3),
      odds_ratio = round(unname(ft$estimate), 4),
      p_value = ft$p.value,
      stringsAsFactors = FALSE)
  }
  if (k == 0L) return(NULL)
  out <- do.call(rbind, rows[seq_len(k)])
  out$p_adj_bh <- stats::p.adjust(out$p_value, method = "BH")
  out <- out[order(out$p_value), , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' 调控子活性分数（按样本）
#'
#' 分数 = 该 TF 靶基因 z-score 的加权均值，符号由 mor 决定。
#' 这是经典的 regulon score，**不是 decoupleR 的 ULM**。
regulon_activity <- function(expr, regulons, detected_genes) {
  # 按基因（行）z-score：让不同基因可比
  z <- t(scale(t(expr)))
  z[!is.finite(z)] <- 0
  det <- intersect(rownames(z), detected_genes)
  z <- z[det, , drop = FALSE]

  tfs <- sort(unique(regulons$tf))
  scores <- list()
  for (tf in tfs) {
    sub <- regulons[regulons$tf == tf, , drop = FALSE]
    sub <- sub[sub$target %in% rownames(z), , drop = FALSE]
    if (nrow(sub) < MIN_TARGETS_FOR_ACTIVITY) next
    m <- z[sub$target, , drop = FALSE]
    # mor 定符号；缺失 mor 时按 +1（激活）处理并记录
    mor <- ifelse(is.na(sub$mor), 1, sign(sub$mor))
    scores[[tf]] <- as.numeric(crossprod(mor, m) / sum(abs(mor)))
  }
  if (length(scores) == 0L) return(NULL)
  out <- as.data.frame(do.call(cbind, scores), stringsAsFactors = FALSE)
  colnames(out) <- names(scores)
  rownames(out) <- colnames(z)
  out
}


#' 组间比较：每个 TF 的活性分数在两组间是否有差异
#'
#' **n<10 时这里只报效应量，不把 p 值当结论。**
#' 用 Wilcoxon（不假设正态）与组间均值差（效应量）两个都报。
group_test <- function(act, group) {
  grp <- as.character(group)
  lv <- sort(unique(grp))
  if (length(lv) != 2L) return(NULL)
  rows <- vector("list", ncol(act))
  for (j in seq_len(ncol(act))) {
    v <- act[[j]]
    g1 <- v[grp == lv[1]]
    g2 <- v[grp == lv[2]]
    wt <- suppressWarnings(stats::wilcox.test(g1, g2, exact = FALSE))
    rows[[j]] <- data.frame(
      tf = colnames(act)[j],
      mean_group1 = round(mean(g1), 4),
      mean_group2 = round(mean(g2), 4),
      # 效应量：Cohen's d（小样本下也只是参考）
      cohens_d = round((mean(g2) - mean(g1)) /
                         sqrt((stats::var(g1) + stats::var(g2)) / 2), 4),
      delta_mean = round(mean(g2) - mean(g1), 4),
      p_value = wt$p.value,
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out$p_adj_bh <- stats::p.adjust(out$p_value, method = "BH")
  out$group1 <- lv[1]
  out$group2 <- lv[2]
  out[order(out$p_value), , drop = FALSE]
}


run_08_tf_regulation <- function(cfg) {
  set.seed(cfg$analysis$seed)
  res <- cfg$output$results_dir
  dat <- cfg$output$data_dir
  status_path <- file.path(res, "tf_status.json")

  # ---- 输入 ---------------------------------------------------------------
  deg_path <- file.path(res, "deg_table.csv")
  expr_path <- file.path(dat, "expr_clean.rds")
  if (!file.exists(deg_path) || !file.exists(expr_path)) {
    write_json(status_path, list(
      status = "not_configured",
      reason = paste("缺少输入：",
                     if (!file.exists(deg_path)) "deg_table.csv" else "",
                     if (!file.exists(expr_path)) "expr_clean.rds" else ""),
      note = "TF 调控分析需要 limma 的差异结果与清洗后的表达矩阵"))
    log_warn("TF 调控：缺输入，记 not_configured")
    return(invisible(NULL))
  }

  deg <- utils::read.csv(deg_path, stringsAsFactors = FALSE)
  expr <- readRDS(expr_path)
  group <- utils::read.csv(file.path(dat, "group.csv"),
                           stringsAsFactors = FALSE)

  # ---- 调控子 -------------------------------------------------------------
  tfdb <- load_tf_regulons(cfg)
  if (is.null(tfdb)) {
    write_json(status_path, list(
      status = "not_done",
      reason = paste(
        "既没有 dorothea 包，也没有 assets/tf_regulons_minimal.tsv。",
        "**没有 TF-靶基因关系就无法做调控分析** —— 不用代理指标冒充。",
        "装 dorothea：BiocManager::install('dorothea')"),
      method = NULL))
    log_warn("TF 调控：没有调控子来源，记 not_done（不用代理指标冒充）")
    return(invisible(NULL))
  }

  detected <- rownames(expr)

  # ---- 分析 1：调控子富集 -------------------------------------------------
  # 用**显著 DEG**；没有显著基因时用最显著的 N 个（与 deg_mode 的措辞约束一致）
  deg_sig <- deg$gene[!is.na(deg$adj.P.Val) & deg$adj.P.Val < 0.05]
  deg_mode <- if (length(deg_sig) > 0L) "significant" else "ranked_fallback"
  if (deg_mode == "ranked_fallback") {
    n_fb <- min(200L, nrow(deg))
    deg_sig <- deg$gene[seq_len(n_fb)]
    log_warn(sprintf(
      "没有 FDR<0.05 的基因 —— 回退到最显著的 %d 个（deg_mode=ranked_fallback）。",
      n_fb))
    log_warn("**措辞约束**：只能说『在最显著的 N 个基因里富集到……』，")
    log_warn("不能说『显著差异基因富集到……』（AGENTS.md 规则 3 / 设计文档 §2.10）")
  }

  enr <- regulon_enrichment(tfdb$regulons, deg_sig, detected)
  if (!is.null(enr)) {
    utils::write.csv(enr, file.path(res, "tf_regulon_enrichment.csv"),
                     row.names = FALSE)
    n_sig <- sum(enr$p_adj_bh < 0.05, na.rm = TRUE)
    log_info(sprintf("调控子富集: %d 个 TF 可检验，BH<0.05 的 %d 个",
                     nrow(enr), n_sig))
    if (nrow(enr) > 0L) {
      log_info(sprintf("  最显著: %s", paste(
        sprintf("%s(p=%.2g)", utils::head(enr$tf, 5),
                utils::head(enr$p_value, 5)), collapse = ", ")))
    }
  } else {
    log_warn("调控子富集：没有 TF 满足最小靶基因数要求")
  }

  # ---- 分析 2：调控子活性 -------------------------------------------------
  act <- regulon_activity(expr, tfdb$regulons, detected)
  act_test <- NULL
  if (!is.null(act)) {
    utils::write.csv(data.frame(sample = rownames(act), act,
                                check.names = FALSE),
                     file.path(res, "tf_activity_by_sample.csv"),
                     row.names = FALSE)
    log_info(sprintf("调控子活性: %d 个样本 x %d 个 TF", nrow(act), ncol(act)))

    # 分组向量与表达矩阵的列对齐。
    # **列名是 GSM 编号，group.csv 里的对应列叫 `gsm`**（不是 `sample`）。
    # 实测用 `group$sample` 时 121 个样本**全部**匹配失败 —— 而代码只是
    # 打了一行 warning 就把 act 清空，组间比较静默变成"什么都没做"。
    # 这正是本仓库反复防的那类失败：产物文件照样产出，内容是空的。
    if (!"gsm" %in% colnames(group)) {
      stop(sprintf("group.csv 里没有 gsm 列（实际列：%s）",
                   paste(colnames(group), collapse = "/")))
    }
    idx <- match(rownames(act), group$gsm)
    g <- group$group[idx]
    n_miss <- sum(is.na(g))
    if (n_miss == nrow(act)) {
      # **全部匹配不上就停下。** 一个 0 样本的"组间比较"不该产出文件。
      stop(sprintf(
        paste0("调控子活性：%d 个样本在 group.csv 里一个都没匹配上。\n",
               "  表达矩阵列名例：%s\n",
               "  group.csv$gsm 例：%s\n",
               "  列名口径不一致 —— 这种情况必须停下来，",
               "否则会产出一份空的组间比较而看起来一切正常。"),
        nrow(act),
        paste(utils::head(rownames(act), 3L), collapse = ", "),
        paste(utils::head(group$gsm, 3L), collapse = ", ")))
    }
    if (n_miss > 0L) {
      log_warn(sprintf("有 %d/%d 个样本在 group.csv 里找不到分组，已剔除",
                       n_miss, nrow(act)))
      keep <- !is.na(g)
      act <- act[keep, , drop = FALSE]
      g <- g[keep]
    }
    act_test <- group_test(act, g)
    if (!is.null(act_test)) {
      utils::write.csv(act_test, file.path(res, "tf_activity_group_test.csv"),
                       row.names = FALSE)
      log_info(sprintf("组间比较: %d 个 TF，BH<0.05 的 %d 个",
                       nrow(act_test), sum(act_test$p_adj_bh < 0.05, na.rm = TRUE)))
    } else {
      log_warn(sprintf(
        "组间比较跳过：分组水平数不是 2（实际 %d 个）—— 写了 tf_status 说明",
        length(unique(g))))
    }
  } else {
    log_warn("调控子活性：没有 TF 满足最小靶基因数要求")
  }

  # ---- 出图 ---------------------------------------------------------------
  # **绘图单独兜住，不让画图错误影响方法本身的记录**（AGENTS.md 规则 14）。
  # 两张图分开写，不用 patchwork —— 少一个依赖，两张图本来也回答两个问题。
  #
  # `figs_written` 收集**实际落盘成功**的图名。写 status 时按它填
  # `figures`，而不是按"代码走到过这个分支"。上一版就是因为按分支填了
  # `figure_written = TRUE` 的兄弟逻辑不严谨（实际报的是 FALSE 但验收
  # 没看它），图没了而验收全绿。
  figs_written <- character(0)
  plot_err <- NULL
  plot_ok <- tryCatch({
    n_show <- min(20L, if (!is.null(enr)) nrow(enr) else 0L)
    if (n_show >= 1L) {
      top <- utils::head(enr, n_show)
      top$tf <- factor(top$tf, levels = rev(top$tf))
      p1 <- ggplot2::ggplot(top, ggplot2::aes(x = .data$tf,
                                              y = -log10(.data$p_value))) +
        ggplot2::geom_col(fill = PAL$up, width = 0.7) +
        ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                            colour = PAL$muted, linewidth = 0.4) +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = sprintf("TF regulon enrichment among DEGs - %s",
                          cfg$dataset_id),
          subtitle = wrap_subtitle(sprintf(
            "Fisher exact test, background = %d detected genes; %s. Dashed line = p 0.05 (uncorrected)",
            length(detected),
            if (identical(deg_mode, "significant"))
              sprintf("%d FDR<0.05 genes", length(deg_sig))
            else sprintf("ranked_fallback: top %d genes (no FDR<0.05)",
                         length(deg_sig)))),
          x = NULL, y = expression(-log[10](italic(p)))) +
        theme_paper()
      # 走仓库自己的 save_pdf()（同时写 PDF 与 PNG，自己管设备与 dpi）。
      # **不要用 ggsave + cfg$analysis$figure_dpi** —— 那个字段不存在，
      # 传 NULL 给 dpi 会报 "`dpi` must be a single number or string"，
      # 被下面的 tryCatch 接住，于是图静默消失、status 里 figure_written=false。
      save_pdf(file.path(res, "01-08-01-unit1-tf-regulon-enrichment.pdf"),
               print(p1), width = W_DOUBLE, height = mm(152))
      figs_written <- c(figs_written, "01-08-01-unit1-tf-regulon-enrichment.pdf")
      if (!is.null(act_test) && nrow(act_test) >= 2L) {
        tt <- utils::head(act_test[order(-abs(act_test$cohens_d)), ], 20L)
        # **按 |d| 单调排序**（评审 3.8：原 factor 按 tt 行序，红蓝交替横跳，
        # |d| 排序完全不可见）。coord_flip 后 levels 顺序 = 图上自下而上，
        # 所以 levels 给 cohens_d 升序、图上自然自上而下 |d| 递减。
        tt$tf <- factor(tt$tf, levels = tt$tf[order(tt$cohens_d)])
        tt$direction <- ifelse(tt$delta_mean > 0,
                               sprintf("higher in %s", tt$group2[1]),
                               sprintf("higher in %s", tt$group1[1]))
        # 方向色板：红 = 在 group2 高，蓝 = 在 group1 高。
        # **同一个红在火山图里也是"上调"** —— 一个颜色一个含义。
        lv2 <- sprintf("higher in %s", tt$group2[1])
        lv1 <- sprintf("higher in %s", tt$group1[1])
        p2 <- ggplot2::ggplot(tt, ggplot2::aes(x = .data$tf,
                                               y = .data$cohens_d,
                                               fill = .data$direction)) +
          ggplot2::geom_col(width = 0.7) +
          # **条上标数值**（评审同一条：效应量要能读出值，不只靠长度）
          ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$cohens_d)),
                             hjust = ifelse(tt$cohens_d >= 0, -0.15, 1.15),
                             size = 2.3, colour = PAL$ink) +
          ggplot2::coord_flip() +
          # 数值标注伸到条外，横轴要留余量
          ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0.10, 0.16))) +
          ggplot2::scale_fill_manual(values = stats::setNames(
            c(PAL$up, PAL$down), c(lv2, lv1))) +
          ggplot2::labs(
            title = sprintf("TF regulon activity, group difference - %s",
                            cfg$dataset_id),
            subtitle = wrap_subtitle(sprintf(
              paste0("Cohen's d (%s vs %s), min n=%d per group. Activity = ",
                     "mor-weighted mean of target-gene z-scores per sample ",
                     "(dorothea regulons). n is small - read effect size, not p"),
              act_test$group2[1], act_test$group1[1], min(table(g)))),
            x = NULL, y = "Cohen's d", fill = NULL) +
          theme_paper()
        save_pdf(file.path(res, "01-08-02-unit1-tf-activity-group-difference.pdf"),
                 print(p2), width = W_DOUBLE, height = mm(152))
        figs_written <- c(figs_written, "01-08-02-unit1-tf-activity-group-difference.pdf")
      }
      TRUE
    } else {
      FALSE
    }
  }, error = function(e) {
    plot_err <<- conditionMessage(e)
    log_warn(sprintf("TF 图绘制失败（方法本身不受影响）: %s", plot_err))
    FALSE
  })

  # 再核一遍：文件真的在磁盘上才算数。
  # **`save_pdf` 自己会兜住 PNG 失败**（PDF 拿到就不中断），所以
  # "没报错"不等于"图在"。以文件系统为准。
  figs_written <- figs_written[
    file.exists(file.path(res, figs_written))]
  figs_missing <- setdiff(
    c("01-08-01-unit1-tf-regulon-enrichment.pdf",
      if (!is.null(act_test) && nrow(act_test) >= 2L)
        "01-08-02-unit1-tf-activity-group-difference.pdf"),
    figs_written)
  if (length(figs_missing) > 0L) {
    log_warn(sprintf("TF 图缺失: %s", paste(figs_missing, collapse = ", ")))
  }
  plot_ok <- length(figs_written) > 0L && length(figs_missing) == 0L

  # ---- 状态 ---------------------------------------------------------------
  n_sig_enr <- if (!is.null(enr)) sum(enr$p_adj_bh < 0.05, na.rm = TRUE) else 0L
  n_sig_act <- if (!is.null(act_test)) {
    sum(act_test$p_adj_bh < 0.05, na.rm = TRUE)
  } else 0L
  n_per_group <- if (exists("g") && length(g) > 0L) {
    as.list(table(g))
  } else {
    NULL
  }

  # ---- §1.5 TRRUST 交叉验证 -------------------------------------------------
  # 失败不抛异常：它是交叉验证，不该拖垮主分析。但**必须留下为什么没做**。
  trrust <- load_trrust_regulons(cfg)
  trrust_cmp <- if (is.null(trrust)) {
    list(compared = FALSE,
         reason = "TRRUST 不可用（下载失败或格式变了）—— 详见上面的 WARN")
  } else {
    compare_regulon_databases(tfdb$regulons, trrust$regulons,
                              a_name = tfdb$source, b_name = trrust$source)
  }
  if (isTRUE(trrust_cmp$compared)) {
    log_info(sprintf(
      "TRRUST vs dorothea：共有 %d 对（Jaccard %.4f），仅 dorothea %d，仅 TRRUST %d",
      trrust_cmp$n_common_pairs, trrust_cmp$jaccard,
      trrust_cmp$n_a_only, trrust_cmp$n_b_only))
    if (!is.na(trrust_cmp$mor_agreement)) {
      log_info(sprintf(
        "  mor 符号一致率 %.4f（%d/%d 可比），**冲突 %d 条**",
        trrust_cmp$mor_agreement, trrust_cmp$n_mor_agree,
        trrust_cmp$n_mor_comparable, trrust_cmp$n_mor_conflict))
    }
  } else {
    log_warn(sprintf("TRRUST 交叉验证未完成：%s", trrust_cmp$reason))
  }

  write_json(status_path, list(
    status = "ok",
    regulon_source = tfdb$source,
    coverage = tfdb$coverage,
    confidence_levels = tfdb$confidence_levels,
    trrust = if (is.null(trrust)) {
      list(status = "not_used",
           reason = "下载失败或上游格式变了（见 WARN 日志）")
    } else {
      list(status = "ok",
           source = trrust$source,
           n_tf = trrust$n_tf,
           n_pairs = trrust$n_pairs,
           n_mor_unknown = trrust$n_mor_unknown,
           url = trrust$url,
           sha256 = trrust$sha256,
           is_primary = FALSE,
           role = "独立交叉验证，不替换主来源")
    },
    trrust_vs_dorothea = trrust_cmp,
    n_tfs_total = length(unique(tfdb$regulons$tf)),
    n_interactions = nrow(tfdb$regulons),
    deg_mode = deg_mode,
    n_genes_in_deg_set = length(deg_sig),
    n_genes_detected = length(detected),
    n_tfs_testable = if (!is.null(enr)) nrow(enr) else 0L,
    n_tfs_enriched_bh05 = n_sig_enr,
    n_tfs_activity_tested = if (!is.null(act_test)) nrow(act_test) else 0L,
    n_tfs_activity_bh05 = n_sig_act,
    n_per_group = n_per_group,
    # 按**实际落盘的文件**填，不按"代码走到过这个分支"
    figure_written = plot_ok,
    figures = as.list(figs_written),
    figures_missing = as.list(figs_missing),
    plot_error = plot_err,
    method = paste(
      "调控子富集 = Fisher 精确检验（背景集 = 检测到的基因）；",
      "调控子活性 = 靶基因 z-score 的 mor 加权均值"),
    not_decoupler = paste(
      "**这不是 decoupleR 的 ULM。** decoupleR 用单样本线性模型 + 置换检验；",
      "本脚本用加权 z-score。两者在**排序上**通常一致，但 p 值与绝对值不同。",
      "要 decoupleR 的结果请另装该包。"),
    limitations = c(
      sprintf("**每组 n 很小**（%s）。调控子活性分数的组间比较只应看效应量，",
              if (!is.null(n_per_group))
                paste(names(n_per_group), unlist(n_per_group),
                      sep = "=", collapse = ", ") else "未知"),
      "      p 值不可作为结论 —— 小样本下 Wilcoxon 的功效极低。",
      sprintf("调控子覆盖度 = %s。", tfdb$coverage),
      if (identical(tfdb$coverage, "minimal"))
        "     **回退表只覆盖少量 TF，结论只能作为示例。**"
      else "     但仍不含低置信度（D/E）的预测调控子。",
      "TF 调控是**关联**，不是因果。表达相关不等于该 TF 在调控这些靶基因。",
      "bulk 数据里 TF 的 mRNA 水平与其**蛋白活性**经常不相关",
      "（翻译后修饰、核转位都不反映在 mRNA 上）—— 这是本分析的根本限制。",
      "靶基因集合之间大量重叠（一个基因受多个 TF 调控），",
      "所以各 TF 的 p 值**不独立**，BH 校正偏保守。",
      if (isTRUE(trrust_cmp$compared) && !is.na(trrust_cmp$mor_agreement))
        sprintf(paste0("**两个库对同一对 TF-靶基因的调控方向就不完全一致**：",
                       "dorothea 与 TRRUST 重叠 %d 对里，mor 符号一致率仅 %.1f%%，",
                       "冲突 %d 条。任何基于 mor 定符号的活性打分都要按此打折。"),
                trrust_cmp$n_mor_comparable,
                100 * trrust_cmp$mor_agreement,
                trrust_cmp$n_mor_conflict)
      else
        "TRRUST 交叉验证未完成，调控方向没有第二个库背书。",
      "**ChEA3 未接入**：它是 Ma'ayan Lab 的 web 服务，无 CRAN/Bioconductor 包，",
      "      调用需要 httr/curl（本仓库 CI 的 R 包列表里没有）。",
      "      而 dorothea 的置信度 A 档本身就整合了 ChEA —— 覆盖率上不缺，",
      "      缺的是「用另一个服务独立复核富集结果」这一步。已如实记在 trrust 字段旁。"),
    outputs = c("tf_regulon_enrichment.csv", "tf_activity_by_sample.csv",
                "tf_activity_group_test.csv", figs_written)))

  log_info("TF 调控分析完成")
  invisible(NULL)
}


if (!GEO_ORCHESTRATED()) {
  cfg <- parse_args()
  run_08_tf_regulation(cfg)
}
