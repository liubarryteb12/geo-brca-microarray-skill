# ============================================================================
# 07_lasso.R — LASSO-Cox 预后签名 + 外部验证
# ============================================================================
# 只在 design_mode: cohort 且 config 显式声明了随访终点时运行。
#
# **终点必须由 config 指定，不自动配对。** 实测 GSE20685 的字段是
# `event_death` 和 `follow_up_duration (years)` —— 名字里没有任何共同词，
# 靠关键词配对（time/duration/survival）在别的数据集上一定会配错，
# 而配错不会报错，只会安静地算出一个错的 C-index。
#
# **EPV 门禁。** 通行判据是每个入选变量至少 10 个事件（events per variable）。
# 实测 GSE42568 的 OS 只有 35 个事件 -> 签名超过 3 个基因就开始过拟合。
# 这里不"禁止"超过 EPV 的签名（文献里到处都是），但会：
#   * 把上限算出来写进 status 和日志；
#   * 用 lambda.1se（而不是 lambda.min）压变量数；
#   * **报外部验证的 C-index**，因为只有它不受过拟合影响。
#
# 输出：results/<GSE>/lasso_coefficients.csv, lasso_risk_scores.csv,
#       lasso_cv_curve.csv, lasso_stability.csv,
#       01-07-01-unit1-lasso-km-training.pdf（+ unit2 验证集，有外部队列时）,
#       lasso_validation.csv（有外部队列时）, lasso_status.json
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

  # 外部验证要复用 step 00/01 的辅助函数：`split_soft_samples()`（00）、
  # `fetch_platform_annotation()` / `map_features_to_symbols()` / `collapse_to_symbol()`（01）。
  # 编排器（main_analysis.R）已经把这两个脚本 source 进来了，所以只在缺失时补。
  # 不补的话 `Rscript scripts/07_lasso.R` 能跑完训练集、然后在外部验证那一步
  # 报 "could not find function"——而那正好是最不该失败的地方。
  if (!exists("collapse_to_symbol", mode = "function")) {
    for (f in c("00_validate_inputs.R", "01_download_clean.R")) {
      p <- file.path(here, f)
      if (file.exists(p)) source(p)
    }
  }
})

EPV_MIN <- 10   # 每个变量至少 10 个事件

#' Harrell's C-index，定义写死在代码里
#'
#' **不用 `survival::concordance()` 的公式接口。** 实测它在
#' `Surv(time, event) ~ risk` 下给出的是 1 - Harrell C：训练集报 **0.121**，
#' 而同一个模型 `cv.glmnet(type.measure="C")` 报 **0.793** —— 两者正好互补。
#'
#' 这个方向约定藏在函数的默认参数里，读调用点看不出来；而错的 C-index
#' 又**看着像个正常数字**（0.121 完全可能是一个真的很差但合理的模型），
#' 所以它不会引起任何怀疑。这正是本仓库一直在防的那类失败。
#'
#' 定义（Harrell 1982，与 `cv.glmnet` 的 `type.measure="C"` 同口径）：
#'   * 可比对的一对 (i, j)：生存时间短的那个**发生了事件**；
#'   * 一致：生存时间短的那个风险分**更高**（Cox 的线性预测子是 log 风险比，
#'     越大越危险）；
#'   * 风险分打平算 0.5；
#'   * C = (一致 + 0.5 x 打平) / 可比对总数。
#'
#' @return 数值；没有可比对样本对时返回 NA
harrell_c <- function(time, event, risk) {
  n <- length(time)
  if (n < 2L) return(NA_real_)
  conc <- 0; tied <- 0; cmp <- 0
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      # 只有"较早发生事件"的那一对才可比对；两个都删失、或较晚的那个先事件，都不可比
      if (time[i] < time[j] && event[i] == 1L) {
        a <- risk[i]; b <- risk[j]
      } else if (time[j] < time[i] && event[j] == 1L) {
        a <- risk[j]; b <- risk[i]
      } else {
        next
      }
      cmp <- cmp + 1
      if (a > b) conc <- conc + 1 else if (a == b) tied <- tied + 1
    }
  }
  if (cmp == 0L) return(NA_real_)
  (conc + 0.5 * tied) / cmp
}

#' 从 clinical.csv 读出随访时间与事件
#'
#' 列名由 config 显式给出；事件取值也显式给出（`event_value`）。
#' 把"什么算事件"写死成 `== 1` 会在事件编码成 "dead"/"recurred" 时静默
#' 把所有样本判成删失，C-index 变成 0.5 而没有任何报错。
#'
#' @return list(time, event, n_events, n_usable) 或 NULL（并记录原因）
read_survival <- function(cfg) {
  dat <- cfg$output$data_dir
  path <- file.path(dat, "clinical.csv")
  if (!file.exists(path)) return(list(err = "clinical.csv 不存在（step 00 未产出）"))
  clin <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  sv <- cfg$survival
  tc <- sv$time_column; ec <- sv$event_column
  if (is.null(tc) || is.null(ec) || !nzchar(tc) || !nzchar(ec)) {
    return(list(err = "config 未指定 survival.time_column / event_column"))
  }
  if (!tc %in% colnames(clin)) {
    return(list(err = sprintf("clinical.csv 里没有时间列「%s」；实际列: %s",
                              tc, paste(utils::head(colnames(clin), 20), collapse = ", "))))
  }
  if (!ec %in% colnames(clin)) {
    return(list(err = sprintf("clinical.csv 里没有事件列「%s」", ec)))
  }

  time  <- suppressWarnings(as.numeric(clin[[tc]]))
  ev_raw <- trimws(as.character(clin[[ec]]))
  ev_levels <- as.character(unlist(sv$event_value %||% list("1")))
  event <- ifelse(ev_raw %in% ev_levels, 1L,
                  ifelse(is.na(ev_raw) | !nzchar(ev_raw), NA_integer_, 0L))

  ok <- is.finite(time) & time > 0 & !is.na(event)
  if (sum(ok) == 0L) {
    return(list(err = sprintf(
      "时间/事件列没有一个可用样本。事件列「%s」的实际取值: %s；config 声明的事件取值: %s",
      ec, paste(utils::head(unique(ev_raw), 8), collapse = ", "),
      paste(ev_levels, collapse = ", "))))
  }
  list(gsm = clin$gsm[ok], time = time[ok], event = event[ok],
       n_usable = sum(ok), n_events = sum(event[ok] == 1L),
       time_column = tc, event_column = ec, event_levels = ev_levels)
}

#' KM 曲线（不依赖 survminer，自己用 survfit + ggplot 画）
#'
#' 风险分组用**训练集的中位数**切，验证集也用它 —— 若在验证集里重新取中位数，
#' 两个队列的"高风险"就不是同一个定义，KM 图看起来能对上而实际上不可比。
make_km_plot <- function(df, cutoff, title, cfg, time_unit = "days") {
  df$stratum <- factor(ifelse(df$risk > cutoff, "high risk", "low risk"),
                       levels = c("low risk", "high risk"))
  fit <- survival::survfit(survival::Surv(time, event) ~ stratum, data = df)
  # 手写 KM 阶梯：survfit 对象里没有现成的 ggplot 接口
  s <- summary(fit)
  steps <- data.frame(
    time    = s$time,
    surv    = s$surv,
    stratum = sub("^stratum=", "", as.character(s$strata)),
    stringsAsFactors = FALSE)
  # 每条曲线的起点 (0, 1)
  starts <- do.call(rbind, lapply(unique(steps$stratum), function(g) {
    data.frame(time = 0, surv = 1, stratum = g, stringsAsFactors = FALSE)
  }))
  steps <- rbind(starts, steps)
  steps <- steps[order(steps$stratum, steps$time), , drop = FALSE]

  # 删失点：从 survfit 的 n.censor 里取，标在曲线上
  cens <- data.frame(time = s$time, n.censor = s$n.censor,
                     stratum = sub("^stratum=", "", as.character(s$strata)),
                     stringsAsFactors = FALSE)
  cens <- cens[cens$n.censor > 0L, , drop = FALSE]
  if (nrow(cens) > 0L) {
    # 找到每个删失时刻对应的生存概率
    cens$surv <- vapply(seq_len(nrow(cens)), function(i) {
      sub <- steps[steps$stratum == cens$stratum[i] & steps$time <= cens$time[i], , drop = FALSE]
      if (nrow(sub) == 0L) NA_real_ else sub$surv[nrow(sub)]
    }, numeric(1))
    cens <- cens[is.finite(cens$surv), , drop = FALSE]
  }

  # log-rank p
  lr <- survival::survdiff(survival::Surv(time, event) ~ stratum, data = df)
  p_lr <- stats::pchisq(lr$chisq, df = length(lr$n) - 1L, lower.tail = FALSE)
  n_hi <- sum(df$stratum == "high risk"); n_lo <- sum(df$stratum == "low risk")

  p <- ggplot2::ggplot(steps, ggplot2::aes(x = time, y = surv, colour = stratum)) +
    ggplot2::geom_step(linewidth = 0.5) +
    ggplot2::scale_colour_manual(
      values = c("low risk" = PAL$down, "high risk" = PAL$up),
      labels = c("low risk" = sprintf("low risk (n=%d)", n_lo),
                 "high risk" = sprintf("high risk (n=%d)", n_hi)),
      name = NULL) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      title = title,
      subtitle = wrap_subtitle(sprintf(
        paste0("median split at risk score = %.3f (cut-off fixed on the TRAINING set). ",
        "log-rank p = %.3g. HR per +1 risk score = %.2f (95%%CI %.2f-%.2f). ",
        "Censored observations are tick marks; "),
        cutoff, p_lr, hr_est, hr_lo, hr_hi), fig_width = W_DOUBLE),
      # **x 轴必须带时间单位**（评审 3.6："time" 无单位 —— 训练队列是天、
      # 验证队列 GSE20685 是年，同一个数字差 365 倍）。
      x = sprintf("Time since diagnosis (%s)", time_unit),
      y = "survival probability") +
    theme_paper(10) +
    ggplot2::theme(legend.position = "bottom")
  if (nrow(cens) > 0L) {
    p <- p + ggplot2::geom_point(data = cens, shape = 124, size = 1.6,
                                 colour = PAL$ink, show.legend = FALSE)
  }
  p
}

#' 下载并清洗外部验证队列的表达矩阵
#'
#' **不走 step 00/01 的完整流程。** 那两个步骤的门禁是"两组设计、每组 >= 10"，
#' 而预后验证队列（如 GSE20685）是**单一队列、全部是癌**，没有对照组 ——
#' 拿它去跑两组门禁只会得到"分组失败"。
#' 这里只需要"表达矩阵 + 临床终点"，所以单独走一条轻量路径，
#' 但探针映射与基因折叠复用 step 01 的函数，保证两个队列的基因空间一致。
#'
#' @return list(expr, clinical) 或 list(err = ...)
load_external_cohort <- function(cfg, gse, platform_id) {
  cache_dir <- file.path(cfg$output$data_dir, paste0("external_", gse))
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  # **每一步都打标签。** 第一次实跑时这一步报了 "argument lengths differ" ——
  # 一条没有任何上下文的 R 错误，既不说哪个函数、也不说哪个对象，
  # 而候选有 getGEO / 探针映射 / 基因折叠 / 临床解析四段。
  # 加标签之后同一个错误会写成 "collapse_to_symbol: argument lengths differ"，
  # 一眼就知道去哪看。这类"错误信息本身没信息量"的情况不值得再猜第二次。
  stage <- function(label, expr) {
    tryCatch(expr, error = function(e)
      stop(sprintf("%s: %s", label, conditionMessage(e)), call. = FALSE))
  }

  eset <- stage("getGEO", {
    e <- tryCatch(
      GEOquery::getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE, destdir = cache_dir, AnnotGPL = FALSE),
      error = function(e) NULL)
    if (is.null(e)) stop(sprintf("GEO 下载 %s 失败", gse))
    if (is.list(e)) e <- e[[1L]]
    e
  })

  expr <- stage("exprs", Biobase::exprs(eset))
  pd <- stage("pData", Biobase::pData(eset))
  gsm <- if ("geo_accession" %in% colnames(pd)) as.character(pd$geo_accession) else rownames(pd)
  if (length(gsm) != ncol(expr)) {
    return(list(err = sprintf("样本 ID 数 %d 与表达矩阵列数 %d 不一致", length(gsm), ncol(expr))))
  }
  colnames(expr) <- gsm
  if (stats::median(expr, na.rm = TRUE) > 50) {
    expr[expr < 0] <- NA
    expr <- log2(expr + 1)
  }

  gpl_id <- Biobase::annotation(eset)
  if (is.null(gpl_id) || !nzchar(gpl_id) || identical(gpl_id, "NA")) {
    gpl_id <- unique(as.character(pd$platform_id))[1L]
  }
  fdata <- stage("fetch_platform_annotation", fetch_platform_annotation(gpl_id, cache_dir))
  if (is.null(fdata) || nrow(fdata) == 0L) {
    return(list(err = sprintf("平台注释 %s 为空", gpl_id)))
  }
  if ("ID" %in% colnames(fdata)) {
    idx <- match(rownames(expr), as.character(fdata$ID))
    fdata <- fdata[idx, , drop = FALSE]
    rownames(fdata) <- rownames(expr)
  } else if (nrow(fdata) != nrow(expr)) {
    # 没有 ID 列时必须按行号对齐 —— 行数不同还继续走，
    # map_features_to_symbols 会拿到一个长度不对的注释表，
    # 而它的 coverage() 用 length(ids) 当分母，长度不对只会算出 >1 的覆盖率，
    # 不会报错。宁可在这里停下。
    return(list(err = sprintf(
      "平台注释没有 ID 列且行数 %d != 探针数 %d，无法对齐", nrow(fdata), nrow(expr))))
  }
  # **`map_features_to_symbols()` 返回的是 list，不是向量。**
  # 第一次实跑时我把它当向量用了，于是 `collapse_to_symbol(expr, <list>)`
  # 报了一句没有任何上下文的 "argument lengths differ"。
  # 加了长度检查后它变成了 "探针->基因映射长度 6 != 探针数 54627" ——
  # 那个 6 正是 list 的元素个数（mode / symbols / method / coverage /
  # mapped_genes / reason），一眼就看出是类型用错了。
  # 正确用法见 01_download_clean.R:370。
  mapping <- stage("map_features_to_symbols", map_features_to_symbols(rownames(expr), fdata))
  if (!identical(mapping$mode, "symbol")) {
    return(list(err = sprintf("平台注释无法映射到基因 symbol：%s", mapping$reason)))
  }
  if (length(mapping$symbols) != nrow(expr)) {
    return(list(err = sprintf("探针->基因映射长度 %d != 探针数 %d",
                              length(mapping$symbols), nrow(expr))))
  }
  expr <- stage("collapse_to_symbol", collapse_to_symbol(expr, mapping$symbols))
  log_info(sprintf("%s 基因折叠后: %d 个基因（映射途径 %s，覆盖率 %.1f%%）",
                   gse, nrow(expr), mapping$method, 100 * mapping$coverage))

  # 临床：复用与 step 00 完全相同的解析器
  gsm_lines <- stage("fetch_geo_soft", fetch_geo_soft(gse, targ = "gsm"))
  blocks <- split_soft_samples(gsm_lines)
  if (length(blocks) == 0L) return(list(err = sprintf("%s 的 SOFT 里没有样本块", gse)))
  block_gsm <- vapply(blocks, function(b) soft_value(b, "Sample_geo_accession"), character(1))
  clinical <- stage("clinical_table", clinical_table(blocks, block_gsm))
  list(expr = expr, clinical = clinical, platform = gpl_id)
}

run_07_lasso <- function(cfg) {
  log_info("=== 步骤 07：LASSO-Cox 预后签名 ===")
  ensure_dirs(cfg)
  res <- cfg$output$results_dir
  status <- list(step = "lasso")

  if (!identical(cfg$design_mode, "cohort")) {
    status$status <- "not_applicable"
    status$reason <- sprintf("design_mode=%s；预后建模需要队列级样本量", cfg$design_mode)
    write_json(file.path(res, "lasso_status.json"), status)
    log_warn(sprintf("跳过 LASSO：%s", status$reason))
    return(invisible(NULL))
  }
  sv <- cfg$survival
  if (is.null(sv) || !isTRUE(sv$enabled)) {
    status$status <- "not_configured"
    status$reason <- "config 没有启用 survival 段（终点必须显式声明，不自动配对）"
    write_json(file.path(res, "lasso_status.json"), status)
    log_warn(sprintf("跳过 LASSO：%s", status$reason))
    return(invisible(NULL))
  }
  for (pkg in c("glmnet", "survival")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      status$status <- "package_missing"
      status$reason <- sprintf("%s 未安装", pkg)
      write_json(file.path(res, "lasso_status.json"), status)
      log_warn(sprintf("跳过 LASSO：%s", status$reason))
      return(invisible(NULL))
    }
  }

  # ---- 1. 终点 ------------------------------------------------------------
  surv <- read_survival(cfg)
  if (!is.null(surv$err)) {
    status$status <- "endpoint_error"
    status$reason <- surv$err
    write_json(file.path(res, "lasso_status.json"), status)
    log_warn(sprintf("跳过 LASSO：%s", surv$err))
    return(invisible(NULL))
  }
  status$time_column <- surv$time_column
  status$event_column <- surv$event_column
  status$event_levels <- as.list(surv$event_levels)
  status$n_usable <- surv$n_usable
  status$n_events <- surv$n_events

  min_events <- as.integer(sv$min_events %||% 20L)
  if (surv$n_events < min_events) {
    status$status <- "too_few_events"
    status$reason <- sprintf("只有 %d 个事件（< %d），Cox 模型无法给出可用估计",
                             surv$n_events, min_events)
    write_json(file.path(res, "lasso_status.json"), status)
    log_warn(sprintf("跳过 LASSO：%s", status$reason))
    return(invisible(NULL))
  }

  # **EPV 门禁**：算出来、写下来，但不假装它能被绕过
  epv_cap <- floor(surv$n_events / EPV_MIN)
  status$epv_min <- EPV_MIN
  status$epv_max_variables <- epv_cap
  log_info(sprintf("EPV: %d 个事件 / %d = 签名最多 %d 个基因（超过即过拟合）",
                   surv$n_events, EPV_MIN, epv_cap))

  # ---- 2. 表达矩阵（只用肿瘤组）-------------------------------------------
  expr  <- readRDS(file.path(cfg$output$data_dir, "expr_clean.rds"))
  group <- utils::read.csv(file.path(cfg$output$data_dir, "group.csv"), stringsAsFactors = FALSE)
  tumor <- group$gsm[group$group == cfg$contrast[1L]]
  common <- Reduce(intersect, list(surv$gsm, tumor, colnames(expr)))
  if (length(common) < 20L) {
    status$status <- "too_few_samples"
    status$reason <- sprintf("同时有表达、属于 %s、且有随访的样本只有 %d 个",
                             cfg$contrast[1L], length(common))
    write_json(file.path(res, "lasso_status.json"), status)
    log_warn(sprintf("跳过 LASSO：%s", status$reason))
    return(invisible(NULL))
  }
  idx <- match(common, surv$gsm)
  time <- surv$time[idx]; event <- surv$event[idx]
  status$n_samples <- length(common)
  status$n_events_used <- sum(event == 1L)
  log_info(sprintf("建模样本: %d 个（%s 组且有随访），其中 %d 个事件",
                   length(common), cfg$contrast[1L], sum(event == 1L)))

  # ---- 3. 候选基因：FDR 显著 DEG ------------------------------------------
  deg <- utils::read.csv(file.path(res, "deg_table.csv"), stringsAsFactors = FALSE)
  sig <- deg[deg$adj.P.Val < cfg$thresholds$adj_p &
             abs(deg$logFC) > cfg$thresholds$log2fc, , drop = FALSE]
  cand <- intersect(sig$gene, rownames(expr))
  status$deg_mode <- if (nrow(sig) > 0L) "fdr" else "ranked_fallback"
  if (length(cand) < 5L) {
    # 与 step 04/05 一致：FDR 为空时退回排序表前 N，但**明确标注**
    n_fb <- as.integer(cfg$analysis$ranked_fallback_genes %||% 500L)
    cand <- intersect(utils::head(deg$gene[order(deg$P.Value)], n_fb), rownames(expr))
    status$deg_mode <- "ranked_fallback"
    log_warn(sprintf("FDR 显著基因为空，LASSO 候选退回 raw P 前 %d 个（结论只能作假设生成）", n_fb))
  }
  status$n_candidates <- length(cand)
  log_info(sprintf("LASSO 候选基因: %d 个（deg_mode=%s）", length(cand), status$deg_mode))

  # ---- 4. 重复 CV ---------------------------------------------------------
  seed <- cfg$analysis$seed
  if (!is.null(seed)) set.seed(seed)
  x <- t(expr[cand, common, drop = FALSE])
  # 标准化：glmnet 默认 standardize=TRUE，但**验证队列必须用训练集的均值方差**，
  # 所以这里自己算并存下来，避免验证时又按验证集的尺度标准化。
  x_center <- colMeans(x); x_scale <- apply(x, 2L, stats::sd)
  x_scale[!is.finite(x_scale) | x_scale == 0] <- 1
  x <- scale(x, center = x_center, scale = x_scale)
  y <- survival::Surv(time, event)

  folds  <- as.integer(sv$cv_folds %||% 10L)
  repeats <- as.integer(sv$cv_repeats %||% 5L)
  foldid_list <- lapply(seq_len(repeats), function(r) {
    # 显式 foldid：不传的话 cv.glmnet 自己抽，重复之间不可比
    sample(rep(seq_len(folds), length.out = length(time)))
  })

  cvs <- lapply(foldid_list, function(fid) {
    glmnet::cv.glmnet(x, y, family = "cox", alpha = 1, foldid = fid, type.measure = "C")
  })
  status$cv_folds <- folds
  status$cv_repeats <- repeats

  # 稳定性：每轮 lambda.1se 选中的基因数与集合
  stab <- do.call(rbind, lapply(seq_along(cvs), function(i) {
    b <- as.matrix(stats::coef(cvs[[i]], s = "lambda.1se"))
    nz <- rownames(b)[b[, 1L] != 0]
    data.frame(repeat_id = i, lambda_1se = cvs[[i]]$lambda.1se, n_selected = length(nz),
               genes = paste(nz, collapse = " | "), stringsAsFactors = FALSE)
  }))
  utils::write.csv(stab, file.path(res, "lasso_stability.csv"), row.names = FALSE)
  status$n_selected_per_repeat <- as.list(stab$n_selected)
  status$n_selected_median <- stats::median(stab$n_selected)

  # 基因入选频率：比"某一轮选了哪些"更能说明稳定性
  all_sel <- unlist(strsplit(stab$genes, " \\| "))
  all_sel <- all_sel[nzchar(all_sel)]
  freq <- sort(table(all_sel), decreasing = TRUE)
  utils::write.csv(data.frame(gene = names(freq), n_repeats = as.integer(freq),
                              stringsAsFactors = FALSE),
                   file.path(res, "lasso_selection_frequency.csv"), row.names = FALSE)

  # 最终模型用第一轮的 foldid 对应的 lambda.1se（显式、可复现）
  final_cv <- cvs[[1L]]
  fit <- glmnet::glmnet(x, y, family = "cox", alpha = 1, lambda = final_cv$lambda.1se)
  beta <- as.matrix(stats::coef(fit))
  coef_df <- data.frame(gene = rownames(beta), coef = beta[, 1L], stringsAsFactors = FALSE)
  coef_df <- coef_df[coef_df$coef != 0, , drop = FALSE]
  coef_df <- coef_df[order(-abs(coef_df$coef)), , drop = FALSE]
  rownames(coef_df) <- NULL
  utils::write.csv(coef_df, file.path(res, "lasso_coefficients.csv"), row.names = FALSE)

  status$lambda_1se <- final_cv$lambda.1se
  status$n_signature_genes <- nrow(coef_df)
  status$exceeds_epv <- nrow(coef_df) > epv_cap

  # **签名稳不稳定，是比"选出了哪些基因"更重要的一件事。**
  # 实测 GSE42568 的 5 轮重复 CV 选出 [16, 3, 3, 22, 3] 个基因 ——
  # 同一个数据集、同一份代码、只换 foldid，基因数差 7 倍。
  # 这时候报"最终签名有 16 个基因"而不报这个分布，等于把不稳定性藏起来。
  status$signature_stable <- length(unique(stab$genes)) == 1L
  if (!status$signature_stable) {
    log_warn(sprintf(
      "签名不稳定：%d 轮重复 CV 各选出 %s 个基因（不同 foldid）—— 最终签名只是其中一轮，不是稳定解",
      repeats, paste(stab$n_selected, collapse = "/")))
  }

  if (nrow(coef_df) == 0L) {
    status$status <- "empty_signature"
    status$reason <- "lambda.1se 下所有系数被压为 0，没有可用签名"
    write_json(file.path(res, "lasso_status.json"), status)
    log_warn(sprintf("LASSO：%s", status$reason))
    return(invisible(NULL))
  }
  if (status$exceeds_epv) {
    log_warn(sprintf("签名 %d 个基因 > EPV 上限 %d —— 训练集 C-index 会偏高，以外部验证为准",
                     nrow(coef_df), epv_cap))
  } else {
    log_info(sprintf("签名 %d 个基因（EPV 上限 %d）", nrow(coef_df), epv_cap))
  }

  # ---- 4b. EPV 合规的替代模型 --------------------------------------------
  #
  # **lambda.1se 不保证满足 EPV。** 实测它给了 16 个基因配 35 个事件
  # （EPV 2.2），远低于通行判据的 10。文献里这种签名到处都是，
  # 但它在独立队列上几乎必然缩水。
  #
  # 所以再给一个**显式满足 EPV 的版本**：取"非零系数 <= epv_cap"的
  # **最大 lambda**（正则化最强的那一端）。两个都落盘、都报 C-index，
  # 让读者自己看代价：少要基因换来多少外部验证性能。
  cap_info <- NULL
  idx_cap <- which(final_cv$nzero <= epv_cap)
  if (length(idx_cap) > 0L && nrow(coef_df) > epv_cap) {
    lam_cap <- final_cv$lambda[max(idx_cap)]   # lambda 递减，max(idx) = 最强正则
    fit_cap <- glmnet::glmnet(x, y, family = "cox", alpha = 1, lambda = lam_cap)
    b_cap <- as.matrix(stats::coef(fit_cap))
    coef_cap <- data.frame(gene = rownames(b_cap), coef = b_cap[, 1L], stringsAsFactors = FALSE)
    coef_cap <- coef_cap[coef_cap$coef != 0, , drop = FALSE]
    coef_cap <- coef_cap[order(-abs(coef_cap$coef)), , drop = FALSE]
    rownames(coef_cap) <- NULL
    if (nrow(coef_cap) > 0L) {
      utils::write.csv(coef_cap, file.path(res, "lasso_coefficients_epv.csv"), row.names = FALSE)
      risk_cap <- as.numeric(predict(fit_cap, newx = x, type = "link"))
      # **基因名要留一份 character 向量给计算用。**
      # `cap_info$genes` 是给 JSON 看的 list，拿它去 subset 矩阵会报
      # "invalid subscript type 'list'" —— 而那个报错发生在外部验证那一段，
      # 看起来像"验证队列有问题"，其实是这里类型存错了。
      cap_genes <- as.character(coef_cap$gene)
      cap_info <- list(lambda = lam_cap, n_genes = nrow(coef_cap),
                       cindex_train = harrell_c(time, event, risk_cap),
                       genes = as.list(cap_genes))
      status$epv_model <- cap_info
      log_info(sprintf("EPV 合规模型: lambda=%.4f, %d 个基因, 训练集 C-index %.3f（EPV>=10 判据）",
                       lam_cap, nrow(coef_cap), cap_info$cindex_train))
    }
  }

  # ---- 4c. 比例风险（PH）假设检验 ----------------------------------------
  #
  # **这一步原来完全没有，而且不是"漏了一行"那么简单。**
  # 全流程用的是 `glmnet(family="cox")` + 手写的 `harrell_c`，
  # **从来没有拟合过一个标准的多因素 Cox 模型** —— 所以 `cox.zph` 根本没有对象可跑。
  # 补这一步等于先把 coxph 拟合出来。
  #
  # PH 是 Cox 模型的核心前提：协变量的 HR 不随时间变化。违反了并不让 C-index 失效，
  # 但**风险分的含义会变成"平均效应"** —— "某基因早期有害、晚期保护"这种情形
  # 会被平均掉，而报告里读起来和真正的恒定效应一模一样。
  #
  # **只在 EPV 合规模型上跑。** 16 基因配 35 个事件时 coxph 的 16 个参数估不准，
  # cox.zph 的输出只是噪声 —— 拿噪声下"PH 成立"的结论比不跑更糟。
  if (!is.null(cap_info)) {
    ph <- tryCatch({
      d <- data.frame(time = time, event = event, stringsAsFactors = FALSE)
      for (g in cap_genes) d[[g]] <- as.numeric(x[, g])
      fit_cox <- survival::coxph(survival::Surv(time, event) ~ ., data = d)
      z <- survival::cox.zph(fit_cox)
      tab <- as.data.frame(z$table)
      tab$term <- rownames(tab)
      rownames(tab) <- NULL
      list(ok = TRUE, table = tab, global_p = unname(z$table["GLOBAL", "p"]))
    }, error = function(e) list(ok = FALSE, reason = conditionMessage(e)))

    if (isTRUE(ph$ok)) {
      utils::write.csv(ph$table, file.path(res, "cox_zph.csv"), row.names = FALSE)
      bad <- setdiff(ph$table$term[ph$table$p < 0.05], "GLOBAL")
      status$ph_assumption <- list(
        status = "ok",
        test = "survival::cox.zph（scaled Schoenfeld 残差）",
        model = "EPV-compliant signature",
        n_genes = length(cap_genes),
        global_p = ph$global_p,
        violated_terms = as.list(bad),
        verdict = if (ph$global_p < 0.05) "violated" else "not_violated")
      if (ph$global_p < 0.05) {
        log_warn(sprintf(
          "PH 假设被拒绝（全局 p=%.4g）：风险分的效应随时间变化，HR 只能当平均效应读",
          ph$global_p))
      } else {
        log_info(sprintf("PH 假设未被拒绝（全局 p=%.3f，Schoenfeld 残差）", ph$global_p))
      }
      if (length(bad) > 0L) {
        log_warn(sprintf("个别协变量 PH 不成立: %s", paste(bad, collapse = ", ")))
      }
    } else {
      status$ph_assumption <- list(status = "failed", reason = ph$reason)
      log_warn(sprintf("cox.zph 失败: %s", ph$reason))
    }
  } else {
    status$ph_assumption <- list(
      status = "not_applicable",
      reason = sprintf(
        "没有 EPV 合规模型（%d 个事件 -> 上限 %d 个变量）；%d 基因签名配 %d 个事件时 coxph 估不准，cox.zph 只会给出噪声",
        surv$n_events, epv_cap, nrow(coef_df), surv$n_events))
    log_warn("跳过 PH 假设检验：没有 EPV 合规模型可作为检验对象")
  }

  # ---- 5. 训练集风险分与 C-index ------------------------------------------
  risk_train <- as.numeric(predict(fit, newx = x, type = "link"))
  names(risk_train) <- common
  c_train <- harrell_c(time, event, risk_train)
  status$cindex_train <- c_train

  # 交叉验证得到的 C-index 是**乐观程度更小**的那个，一并报出来
  # 用 which.min(abs(...)) 而不是 which(lambda == lambda.1se)：
  # 后者在浮点不完全相等时返回 integer(0)，status 里就会出现一个空的
  # cindex_cv，而空值在 JSON 里看着像"没算"，不像"算错了"。
  status$cindex_cv <- unname(final_cv$cvm[which.min(abs(final_cv$lambda - final_cv$lambda.1se))])
  log_info(sprintf("C-index: 训练集 %.3f | 交叉验证 %.3f",
                   status$cindex_train, status$cindex_cv))

  # ---- 6. 外部验证 --------------------------------------------------------
  val_gse <- sv$validation_dataset
  cutoff <- stats::median(risk_train)
  risk_df <- data.frame(gsm = common, time = time, event = event, risk = risk_train,
                        set = "training", stringsAsFactors = FALSE)

  if (!is.null(val_gse) && nzchar(val_gse)) {
    val <- tryCatch(load_external_cohort(cfg, val_gse, sv$validation_platform_id),
                    error = function(e) list(err = conditionMessage(e)))
    if (!is.null(val$err)) {
      status$validation <- "failed"
      status$validation_error <- val$err
      log_warn(sprintf("外部验证失败（签名本身有效）: %s", val$err))
    } else {
      vexpr <- val$expr; vclin <- val$clinical
      vs <- read_survival_from_table(vclin, sv, val_gse)
      if (!is.null(vs$err)) {
        status$validation <- "endpoint_error"
        status$validation_error <- vs$err
        log_warn(sprintf("外部验证：%s", vs$err))
      } else {
        genes <- coef_df$gene
        miss <- setdiff(genes, rownames(vexpr))
        have <- intersect(genes, rownames(vexpr))
        vcommon <- Reduce(intersect, list(vs$gsm, colnames(vexpr)))
        status$validation_dataset <- val_gse
        status$validation_n_samples <- length(vcommon)
        status$validation_n_events <- sum(vs$event[match(vcommon, vs$gsm)] == 1L)
        status$validation_genes_found <- length(have)
        status$validation_genes_missing <- as.list(miss)

        if (length(have) == 0L || length(vcommon) < 10L) {
          status$validation <- "unusable"
          log_warn("外部验证：可用样本或基因不足")
        } else {
          # **用训练集的均值方差标准化验证集**，否则两个队列的风险分不在同一尺度上，
          # 中位数切点也就失去意义。
          vx <- t(vexpr[have, vcommon, drop = FALSE])
          vx <- scale(vx, center = x_center[have], scale = x_scale[have])
          vx[!is.finite(vx)] <- 0
          b <- coef_df$coef[match(have, coef_df$gene)]
          risk_val <- as.numeric(vx %*% b)
          vi <- match(vcommon, vs$gsm)
          c_val <- harrell_c(vs$time[vi], vs$event[vi], risk_val)
          status$validation <- "ok"
          status$cindex_validation <- c_val
          log_info(sprintf("外部验证 %s: n=%d, %d 个事件, C-index = %.3f（训练集 %.3f）",
                           val_gse, length(vcommon), status$validation_n_events,
                           c_val, status$cindex_train))

          # EPV 合规模型也在同一个外部队列上打分 —— 只有这样才能回答
          # "少要 13 个基因换来多少外部性能"。
          if (!is.null(cap_info)) {
            have_cap <- intersect(cap_genes, rownames(vexpr))
            if (length(have_cap) > 0L) {
              vxc <- t(vexpr[have_cap, vcommon, drop = FALSE])
              vxc <- scale(vxc, center = x_center[have_cap], scale = x_scale[have_cap])
              vxc[!is.finite(vxc)] <- 0
              bc <- coef_cap$coef[match(have_cap, coef_cap$gene)]
              risk_cap_val <- as.numeric(vxc %*% bc)
              c_cap <- harrell_c(vs$time[vi], vs$event[vi], risk_cap_val)
              status$epv_model$cindex_validation <- c_cap
              status$epv_model$genes_found_in_validation <- length(have_cap)
              log_info(sprintf("EPV 合规模型外部验证 C-index = %.3f（%d/%d 个基因在验证队列中找到）",
                               c_cap, length(have_cap), length(cap_info$genes)))
            }
          }
          risk_df <- rbind(risk_df,
                           data.frame(gsm = vcommon, time = vs$time[vi], event = vs$event[vi],
                                      risk = risk_val, set = val_gse, stringsAsFactors = FALSE))
        }
      }
    }
  } else {
    status$validation <- "not_configured"
    log_warn("未配置外部验证队列；训练集 C-index 不能作为签名性能的证据")
  }

  utils::write.csv(risk_df, file.path(res, "lasso_risk_scores.csv"), row.names = FALSE)

  # **签名基因的训练集表达矩阵落盘**（供 10 的 risk-plot 层3 热图读取；
  # 10 不读表达矩阵本体，只读这张小表）。矩阵 = EPV 合规模型的基因 x 训练样本。
  if (!is.null(cap_info)) {
    sig_in <- intersect(cap_genes, rownames(expr))
    if (length(sig_in) >= 1L) {
      utils::write.csv(expr[sig_in, common, drop = FALSE],
                       file.path(res, "signature_expr.csv"), row.names = TRUE)
    }
  }
  utils::write.csv(data.frame(lambda = final_cv$lambda, cvm = final_cv$cvm,
                              cvsd = final_cv$cvsd, nzero = final_cv$nzero),
                   file.path(res, "lasso_cv_curve.csv"), row.names = FALSE)

  # KM：训练集与验证集**各出一张单图**（`unit1` / `unit2` 共用一个图号），
  # 切点都用训练集的中位数（跨队列可比）。
  #
  # **原来是把两张图 `print()` 进同一个 `save_pdf`。** 那会写出一个
  # **两页的 PDF**，而 PNG 只留得下其中一页 —— 读者拿到 PDF 看到两张、
  # 拿到 PNG 只看到一张，而两边的文件名是同一个。
  # 拆成单图之后两边一致（这也是仓库出图约定里"尽量出单图"的一条）。
  save_pdf(file.path(res, "01-07-01-unit1-lasso-km-training.pdf"), {
    print(make_km_plot(risk_df[risk_df$set == "training", , drop = FALSE], cutoff,
                       sprintf("LASSO-Cox risk groups (training) - %s", cfg$dataset_id), cfg,
                       time_unit = "days"))
  }, width = W_DOUBLE, height = mm(70))
  if (any(risk_df$set != "training")) {
    vset <- unique(risk_df$set[risk_df$set != "training"])[1L]
    save_pdf(file.path(res, "01-07-01-unit2-lasso-km-validation.pdf"), {
      print(make_km_plot(risk_df[risk_df$set == vset, , drop = FALSE], cutoff,
                         sprintf("LASSO-Cox risk groups (external validation: %s)", vset), cfg,
                         time_unit = "years"))
    }, width = W_DOUBLE, height = mm(70))
  }

  # ---- 8. 推荐哪一个签名 --------------------------------------------------
  #
  # **EPV 超标不做硬拒绝，而是明确"推荐哪一个"。**
  # 我最初的想法是把超标直接拒掉，想清楚后改了 —— 那是错的：
  # 16 基因配 35 个事件是 glmnet 在 lambda.1se 下的**真实输出**，
  # 把它删掉恰好掩盖了"这个方法在这个样本量下就会产出过拟合签名"这个事实本身，
  # 而那是这份分析最该说出来的话。
  #
  # 真正的问题是**两个签名都报 C-index 却不说哪个该用**，读者会拿错那个看
  # （16 基因训练集 0.879 看着比 3 基因的 0.806 漂亮得多，而外部验证是反的）。
  # 所以加机器可读的推荐字段，而不是指望读者自己判断。
  if (!isTRUE(status$exceeds_epv)) {
    status$signature_recommended <- "lasso_1se"
    status$recommendation_reason <- sprintf(
      "lambda.1se 签名 %d 个基因配 %d 个事件，EPV = %.1f >= %d，本身已合规",
      nrow(coef_df), surv$n_events, surv$n_events / max(nrow(coef_df), 1L), EPV_MIN)
  } else if (!is.null(cap_info)) {
    status$signature_recommended <- "epv_compliant"
    status$recommendation_reason <- sprintf(
      "lambda.1se 签名 %d 个基因配 %d 个事件，EPV = %.1f < %d（超标）；推荐 EPV 合规的 %d 基因版本。两个都落盘并各报三个 C-index，代价可自行核对",
      nrow(coef_df), surv$n_events, surv$n_events / max(nrow(coef_df), 1L), EPV_MIN,
      cap_info$n_genes)
  } else {
    status$signature_recommended <- "none"
    status$recommendation_reason <- sprintf(
      "lambda.1se 签名 %d 个基因配 %d 个事件，EPV = %.1f < %d，且没有可用的 EPV 合规替代（正则化最强的 lambda 下仍有 %d 个非零系数）。这份签名不应作为预后模型使用",
      nrow(coef_df), surv$n_events, surv$n_events / max(nrow(coef_df), 1L), EPV_MIN, epv_cap)
  }
  log_info(sprintf("推荐签名: %s —— %s", status$signature_recommended,
                   status$recommendation_reason))

  status$status <- "ok"
  write_json(file.path(res, "lasso_status.json"), status)
  log_info(paste0("已生成 lasso_coefficients.csv / lasso_risk_scores.csv / lasso_stability.csv / ",
                  "lasso_cv_curve.csv / 01-07-01-unit1-lasso-km-training.pdf",
                  "（+ unit2 验证集）/ lasso_status.json"))
  invisible(coef_df)
}

#' 在任意临床表上按 config 的列名取终点（外部队列用）
read_survival_from_table <- function(clin, sv, gse) {
  tc <- sv$validation_time_column %||% sv$time_column
  ec <- sv$validation_event_column %||% sv$event_column
  if (!tc %in% colnames(clin)) {
    return(list(err = sprintf("%s 没有时间列「%s」；实际列: %s", gse, tc,
                              paste(utils::head(colnames(clin), 20), collapse = ", "))))
  }
  if (!ec %in% colnames(clin)) {
    return(list(err = sprintf("%s 没有事件列「%s」", gse, ec)))
  }
  time <- suppressWarnings(as.numeric(clin[[tc]]))
  ev_raw <- trimws(as.character(clin[[ec]]))
  lv <- as.character(unlist(sv$validation_event_value %||% sv$event_value %||% list("1")))
  event <- ifelse(ev_raw %in% lv, 1L,
                  ifelse(is.na(ev_raw) | !nzchar(ev_raw), NA_integer_, 0L))
  ok <- is.finite(time) & time > 0 & !is.na(event)
  if (sum(ok) == 0L) {
    return(list(err = sprintf("%s 的时间/事件列没有可用样本（事件列取值: %s）", gse,
                              paste(utils::head(unique(ev_raw), 8), collapse = ", "))))
  }
  list(gsm = clin$gsm[ok], time = time[ok], event = event[ok],
       n_events = sum(event[ok] == 1L), time_column = tc, event_column = ec)
}

if (!GEO_ORCHESTRATED()) {
  cfg <- load_config()
  run_07_lasso(cfg)
}
