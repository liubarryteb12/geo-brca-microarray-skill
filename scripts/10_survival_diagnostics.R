# ============================================================================
# 10_survival_diagnostics.R — 时间依赖 AUC 与校准（§1.6）
#
# ## 为什么单独一个脚本，而不是塞进 07_lasso.R
#
# §1.6 点名了 `timeROC` 与 `rms`，而 07_lasso.R 只报了 Harrell C-index。
# **C-index 是一个数，它不告诉你模型在哪个时间段有用。**
# 一个签名完全可能整体 C-index 0.70、而 5 年 AUC 只有 0.55 —— 前者会被
# 读成"模型有效"，后者才是临床上真正关心的（5 年生存是乳腺癌的决策点）。
#
# 放进 07 的风险是**改坏已经验证过的数值**：07 有 700 多行，任何一行
# 动到 `risk_train` 都会让 C-index 变，而 C-index 变了不容易察觉。
# 所以这一步**只读** `lasso_risk_scores.csv`（07 的产物），一行都不改 07。
#
# ## 两个诊断，各答一个问题
#
# | 工具 | 问题 | 产物 |
# |---|---|---|
# | `timeROC` | **在哪个时间点**模型有区分度？ | `time_roc.csv`（各时间点 AUC + CI） |
# | `rms`     | 预测的**绝对风险**准不准？ | `calibration.csv` + `rms_validate.csv` |
#
# **区分度（AUC）和校准（calibration）是两件事，缺一不可。**
# 一个模型可以把所有人排序排得很准（AUC 高），同时把每个人的风险都
# 高估一倍（校准差）—— 而后者直接影响"要不要化疗"这类决策。
# 只报 AUC 等于只报了一半。
#
# ## 时间单位不跨队列共享 —— 这是最容易错的地方
#
# 主队列 `overall survival time_days` 是**天**，验证队列
# `follow_up_duration (years)` 是**年**。同一个数字 1825 在两个队列里
# 差 365 倍。所以时间点是**每个队列各自按事件时间分位数算的**，
# 不是全局常量。共用会把验证队列的 AUC 算在一个只有极少数人随访到的点上。
#
# ## 产物
#
#   results/<GSE>/time_roc.csv            各队列、各时间点的 AUC 与 CI
#   results/<GSE>/calibration.csv         预测风险 vs KM 观测生存（分位分组）
#   results/<GSE>/rms_validate.csv        乐观校正后的 Dxy / 斜率 / R2
#   results/<GSE>/survival_diagnostics_status.json
#   figures/time_roc.png / calibration.png
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

RISK_FILE <- "lasso_risk_scores.csv"

# 一个时间点至少要有这么多事件才报 AUC。**这是硬门槛，不是可调参数。**
# 实测：时间依赖 AUC 在事件数少时方差极大，3 个事件能给出 AUC 0.95 或 0.30。
# 报一个只有 4 个事件的 5 年 AUC，读者会把它当成"5 年区分度"。
MIN_EVENTS_AT_HORIZON <- 10L

# 时间点取事件时间分位数。**按队列各自算** —— 单位不同（天 vs 年）。
HORIZON_QUANTILES <- c(0.25, 0.50, 0.75)

# rms 自助法的重抽样次数。200 是 rms 文档的默认量级；
# 它决定乐观校正的稳定性，和运行时间。
RMS_BOOTSTRAP <- 200L

# 校准曲线的分组数。分位分组，每组人数尽量相等。
CALIBRATION_GROUPS <- 4L


#' 从事件时间分位数挑可用时间点
#'
#' @return list(horizons, note) —— horizons 可能为空（事件太少）
pick_horizons <- function(time, event) {
  ev_times <- time[event == 1L]
  if (length(ev_times) < MIN_EVENTS_AT_HORIZON) {
    return(list(horizons = numeric(0),
                note = sprintf("事件只有 %d 个，不足 %d 个，不报时间依赖 AUC",
                               length(ev_times), MIN_EVENTS_AT_HORIZON)))
  }
  cand <- unique(stats::quantile(ev_times, HORIZON_QUANTILES,
                                 names = FALSE, type = 1))
  cand <- sort(cand[cand > 0])
  # 每个时间点后面还剩多少事件 —— 这是 AUC 能不能报的判据
  keep <- vapply(cand, function(h) sum(ev_times > h) >= MIN_EVENTS_AT_HORIZON,
                 logical(1))
  dropped <- cand[!keep]
  cand <- cand[keep]
  note <- if (length(dropped) > 0L) {
    sprintf("丢弃时间点 %s（其后事件不足 %d 个）",
            paste(format(dropped, digits = 4), collapse = ", "),
            MIN_EVENTS_AT_HORIZON)
  } else {
    "全部分位时间点都保留了"
  }
  list(horizons = cand, note = note)
}


#' 一个队列的时间依赖 AUC
#'
#' **名字不叫 `run_*`** —— `check_r_syntax.mjs` 把 `run_` 前缀当作步骤函数，
#' 会报"定义了但从未被引用"（它只扫编排器里的调用，看不见内部调用）。
#'
#' `timeROC` 的 `iid = TRUE` 才会给 `inference`（标准误），
#' 不给 SE 就只能报点估计 —— 而 AUC 的点估计单独看没有意义。
compute_time_roc <- function(df, label) {
  time <- df$time; event <- as.integer(df$event); risk <- df$risk
  ph <- pick_horizons(time, event)
  out <- list(label = label, horizons = ph$horizons, note = ph$note,
              rows = list(), error = NULL)
  if (length(ph$horizons) == 0L) return(out)

  res <- tryCatch(
    timeROC::timeROC(T = time, delta = event, marker = risk, cause = 1,
                     times = ph$horizons, iid = TRUE),
    error = function(e) e)
  if (inherits(res, "error")) {
    out$error <- conditionMessage(res)
    return(out)
  }

  auc <- as.numeric(res$AUC)
  # inference 是 2 x n 矩阵：第 1 行标准误，第 2 行是……（不同版本不同）
  # 所以**按名字取**，不按位置取 —— 位置在不同 timeROC 版本间变过。
  se <- NULL
  if (!is.null(res$inference)) {
    inf <- res$inference
    if (is.matrix(inf) && "SE" %in% rownames(inf)) {
      se <- as.numeric(inf["SE", ])
    } else if (is.matrix(inf) && nrow(inf) >= 1L) {
      se <- as.numeric(inf[1L, ])
    } else if (is.numeric(inf)) {
      se <- as.numeric(inf)
    }
  }
  for (i in seq_along(ph$horizons)) {
    lo <- hi <- NA_real_
    if (!is.null(se) && length(se) >= i && is.finite(se[[i]])) {
      lo <- auc[[i]] - 1.96 * se[[i]]
      hi <- auc[[i]] + 1.96 * se[[i]]
    }
    out$rows[[length(out$rows) + 1L]] <- list(
      set = label, horizon = ph$horizons[[i]], auc = auc[[i]],
      se = if (is.null(se)) NA_real_ else se[[i]],
      ci_low = lo, ci_high = hi,
      n_at_risk = sum(time >= ph$horizons[[i]]),
      n_events_after = sum(time[event == 1L] > ph$horizons[[i]]))
  }
  out
}


#' 一个队列的校准（预测风险 vs KM 观测生存）
#'
#' **用 KM 算观测生存，不用"有没有事件"的二元比例。**
#' 到时间点 t 时删失的人不是"没发生事件"，直接算比例会系统性低估风险 ——
#' 而低估的方向恰好让校准曲线看起来"预测偏高"，是会被误读成模型保守的。
calibrate_set <- function(df, label, horizon) {
  time <- df$time; event <- as.integer(df$event); risk <- df$risk
  n <- length(time)
  if (n < 20L) return(list(rows = list(), note = "样本不足 20，不做校准"))
  if (sum(time[event == 1L] > horizon) < MIN_EVENTS_AT_HORIZON) {
    return(list(rows = list(),
                note = sprintf("时间点 %g 之后事件不足，不做校准", horizon)))
  }

  # Cox 只用来把 risk 变成绝对风险；区分度已经由 timeROC 答过了
  fit <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ risk),
    error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(rows = list(), note = paste("Cox 拟合失败:", conditionMessage(fit))))
  }

  # 基线生存 -> 每个个体的预测生存
  sf <- survival::survfit(fit, newdata = data.frame(risk = risk))
  # survfit 在 newdata 下给矩阵：行 = 时间点，列 = 个体
  idx <- findInterval(horizon, sf$time)
  pred <- if (idx <= 0L) rep(1, n) else as.numeric(sf$surv[idx, ])
  pred <- pmin(pmax(pred, 1e-6), 1 - 1e-6)

  # 按预测风险分位分组
  g <- min(CALIBRATION_GROUPS, max(2L, floor(n / 10)))
  brk <- unique(stats::quantile(pred, seq(0, 1, length.out = g + 1L),
                                names = FALSE))
  if (length(brk) < 3L) {
    return(list(rows = list(), note = "预测风险几乎无变异，无法分组"))
  }
  grp <- cut(pred, breaks = brk, include.lowest = TRUE, labels = FALSE)

  rows <- list()
  for (k in sort(unique(grp[!is.na(grp)]))) {
    sel <- which(grp == k)
    km <- tryCatch(
      survival::survfit(survival::Surv(time[sel], event[sel]) ~ 1),
      error = function(e) e)
    if (inherits(km, "error")) next
    j <- findInterval(horizon, km$time)
    obs <- if (j <= 0L) 1 else as.numeric(km$surv[j])
    rows[[length(rows) + 1L]] <- list(
      set = label, horizon = horizon, group = k, n = length(sel),
      predicted = mean(pred[sel]), observed = obs,
      n_events = sum(event[sel] == 1L))
  }
  list(rows = rows, note = sprintf("%d 组，Cox 基线生存 + KM 观测", g))
}


#' rms 乐观校正（Dxy / 校准斜率 / R2）
#'
#' **为什么要乐观校正：** 在训练集上算校准斜率，得到的必然接近 1 ——
#' 因为模型就是在这批人身上拟合的。自助法估计"这个 1 里有多少是乐观偏差"，
#' 减掉之后才是真实斜率。不校正的斜率报出来是自我循环。
rms_validate_set <- function(df, label) {
  time <- df$time; event <- as.integer(df$event); risk <- df$risk
  if (sum(event) < MIN_EVENTS_AT_HORIZON) {
    return(list(row = NULL, note = "事件太少，不做 rms 校正"))
  }
  d <- data.frame(time = time, event = event, risk = risk)
  res <- tryCatch({
    fit <- rms::cph(survival::Surv(time, event) ~ risk, data = d,
                    x = TRUE, y = TRUE, surv = TRUE)
    v <- rms::validate(fit, method = "boot", B = RMS_BOOTSTRAP)
    list(fit = fit, v = v)
  }, error = function(e) e)
  if (inherits(res, "error")) {
    return(list(row = NULL, note = paste("rms 失败:", conditionMessage(res))))
  }
  v <- res$v
  # validate 的表：行 = 指标，列 = index.orig / training / test / optimism /
  # index.corrected。**取 index.corrected，不是 index.orig** ——
  # 后者就是上面说的自我循环。
  getm <- function(nm) {
    if (!nm %in% rownames(v)) return(NA_real_)
    if ("index.corrected" %in% colnames(v)) {
      as.numeric(v[nm, "index.corrected"])
    } else {
      as.numeric(v[nm, ncol(v)])
    }
  }
  geto <- function(nm) {
    if (!nm %in% rownames(v)) return(NA_real_)
    as.numeric(v[nm, "index.orig"])
  }
  list(row = list(set = label, n = length(time), n_events = sum(event),
                  Dxy_corrected = getm("Dxy"), Dxy_orig = geto("Dxy"),
                  slope_corrected = getm("Slope"), slope_orig = geto("Slope"),
                  R2_corrected = getm("R2"), R2_orig = geto("R2"),
                  B = RMS_BOOTSTRAP),
       note = sprintf("B=%d 自助法，取 index.corrected", RMS_BOOTSTRAP))
}


#' 主入口
run_10_survival_diagnostics <- function(cfg) {
  res <- cfg$output$results_dir
  # 本仓库的图落在 results_dir，**没有独立的 figures_dir**
  # （`cfg$output` 只有 results_dir / data_dir 两个字段）。
  fig <- res
  ensure_dirs(cfg)

  # ---- 临时挂载 survival（本仓库第二处 library()，理由见下）----------------
  #
  # **`timeROC::timeROC()` 不加这一步会直接报 `could not find function "Surv"`。**
  # 实测（run 35483432635）：两个队列的 timeROC 全部失败，而状态文件当时
  # 写的是"事件数不足" —— 一次崩溃伪装成了"这一步不适用"。
  #
  # **根因是 timeROC 用了未声明的依赖。** 查过 CRAN 上 timeROC 0.4.1 的
  # DESCRIPTION 与 NAMESPACE：
  #
  #   Depends:  R (>= 2.10)                    <- 没有 survival
  #   Imports:  pec (>= 2.4.4), mvtnorm        <- 没有 survival
  #   Suggests: survival, timereg              <- survival 在这里
  #   NAMESPACE: import(pec); import(mvtnorm)
  #              **没有任何 importFrom(survival, ...)**
  #
  # 也就是说 `timeROC` 内部按名字调用 `Surv`，却把 `survival` 只写成
  # `Suggests`。正常用法（用户先 `library(survival)` 再 `library(timeROC)`）
  # 能跑通，是因为 `survival` 恰好在搜索路径上。
  #
  # **`pkg::fun()` 不 attach 任何东西**（连 `Depends` 也不 attach，更别说
  # `Suggests`），所以本仓库的调用约定下它必然失败。
  # 从外面没有参数能改 —— 只能在调用期间把 `survival` 放上搜索路径。
  #
  # 本仓库禁止 `library()`（AGENTS「禁止」段），唯一例外是规则 23 的
  # `blockwiseModules`。这是第二处，同样必须：
  #   1. 先记 `"package:survival" %in% search()`，避免拆掉调用方原有状态
  #   2. `on.exit` 立刻 detach（`unload = FALSE` —— 别的包可能还在用它的
  #      命名空间，卸载会引发难以定位的副作用）
  # 脚本其余所有调用仍然写全名（`survival::coxph` / `survival::survfit` 等）。
  surv_attached <- "package:survival" %in% search()
  if (!surv_attached) {
    suppressPackageStartupMessages(library(survival))
    on.exit({
      if ("package:survival" %in% search()) {
        detach("package:survival", unload = FALSE, character.only = TRUE)
      }
    }, add = TRUE)
    log_info("临时挂载 survival（timeROC 把它写在 Depends 里，:: 不会 attach）")
  }

  status <- list(
    dataset_id = cfg$dataset_id,
    status = "ok",
    method = "timeROC（时间依赖 AUC）+ rms::cph/validate（乐观校正）",
    risk_file = RISK_FILE,
    min_events_at_horizon = MIN_EVENTS_AT_HORIZON,
    horizon_quantiles = HORIZON_QUANTILES,
    rms_bootstrap = RMS_BOOTSTRAP,
    limitations = c(
      "**时间依赖 AUC 与校准都在同一批样本上算的**（训练集）。它们描述的是模型在训练数据上的行为，不是外部性能 —— 外部验证看 07_lasso.R 的 validation 段。",
      "**校准用 Cox 基线生存 + KM 观测**，不是 rms::calibrate 的分组自助法曲线；后者需要 datadist 全局选项，会把状态带进后续调用。斜率与 Dxy 仍走 rms::validate。",
      "**每个队列的时间点各自按事件时间分位数取**（主队列单位是天、验证队列是年），跨队列共享同一组数字会把 AUC 算在几乎没人随访到的点上。",
      "**时间点后事件不足 10 个就不报 AUC**。事件少时时间依赖 AUC 的方差极大，点估计会被读成结论。",
      "校准分组数固定为 4，样本少时组内人数少，KM 阶梯的跳跃会直接体现在校准曲线上。",
      "**这是回顾性数据的内部诊断，不是临床决策依据**；EPV 限制见 07_lasso.R 与 config 注释。"
    )
  )

  rp <- file.path(res, RISK_FILE)
  if (!file.exists(rp)) {
    status$status <- "not_configured"
    status$reason <- sprintf("%s 不存在（07_lasso.R 未产出或未配置随访终点）", RISK_FILE)
    log_warn(status$reason)
    write_json(file.path(res, "survival_diagnostics_status.json"), status)
    return(status)
  }

  risk_df <- utils::read.csv(rp, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("time", "event", "risk", "set")
  miss <- setdiff(need, colnames(risk_df))
  if (length(miss) > 0L) {
    status$status <- "schema_error"
    status$reason <- sprintf("%s 缺列: %s；实际列: %s", RISK_FILE,
                             paste(miss, collapse = ", "),
                             paste(colnames(risk_df), collapse = ", "))
    log_error(status$reason)
    write_json(file.path(res, "survival_diagnostics_status.json"), status)
    return(status)
  }

  sets <- unique(as.character(risk_df$set))
  log_info(sprintf("风险分数表: %d 行, 队列: %s", nrow(risk_df),
                   paste(sets, collapse = ", ")))

  roc_all <- list(); roc_notes <- list()
  for (s in sets) {
    d <- risk_df[risk_df$set == s, , drop = FALSE]
    d <- d[is.finite(d$time) & d$time > 0 & !is.na(d$event) & is.finite(d$risk), ,
           drop = FALSE]
    if (nrow(d) < 20L) {
      roc_notes[[s]] <- sprintf("n=%d，不足 20，跳过", nrow(d))
      next
    }
    # rms 的自助法要用 RNG —— 紧挨着 set.seed（规则 11）
    set.seed(cfg$analysis$seed)
    r <- compute_time_roc(d, s)
    roc_notes[[s]] <- r$note
    if (!is.null(r$error)) {
      log_warn(sprintf("%s 的 timeROC 失败: %s", s, r$error))
      roc_notes[[s]] <- paste(r$note, "| timeROC 失败:", r$error)
    }
    if (length(r$rows) > 0L) roc_all <- c(roc_all, r$rows)
    log_info(sprintf("%s: %d 个时间点 AUC", s, length(r$rows)))
    for (row in r$rows) {
      log_info(sprintf("  t=%g  AUC=%.4f  95%%CI=[%.4f, %.4f]  风险集 %d 事件 %d",
                       row$horizon, row$auc, row$ci_low, row$ci_high,
                       row$n_at_risk, row$n_events_after))
    }
  }

  if (length(roc_all) == 0L) {
    # **崩溃和"事件不够"必须长得不一样。**
    #
    # 实测踩过（run 35483432635）：timeROC 因为 `Surv` 找不到而全部失败，
    # 而这里无条件写 `too_few_events` + "事件数不足" —— 一次**代码崩溃**
    # 伪装成了"这一步不适用"，验收照常 PASS。这正是 AGENTS 规则 24
    # （可选步骤失败不等于这一步不适用）说的那个坑，我自己又踩了一次。
    #
    # 所以先看有没有 error：有 error 就是 failed，理由里带上错误原文。
    errs <- roc_notes[!vapply(roc_notes, function(x) !grepl("失败:", x), logical(1))]
    status$roc_errors <- as.list(errs)
    if (length(errs) > 0L) {
      status$status <- "failed"
      status$reason <- sprintf(
        "timeROC 在 %d 个队列上失败（**不是事件不足**）：%s",
        length(errs), paste(unlist(errs), collapse = " | "))
      log_error(status$reason)
    } else {
      status$status <- "too_few_events"
      status$reason <- "所有队列的事件数都不足以算时间依赖 AUC（timeROC 未报错）"
      log_warn(status$reason)
    }
    status$notes <- roc_notes
    write_json(file.path(res, "survival_diagnostics_status.json"), status)
    return(status)
  }

  roc_df <- do.call(rbind, lapply(roc_all, function(r) {
    data.frame(set = r$set, horizon = r$horizon, auc = r$auc, se = r$se,
               ci_low = r$ci_low, ci_high = r$ci_high,
               n_at_risk = r$n_at_risk, n_events_after = r$n_events_after,
               stringsAsFactors = FALSE)
  }))
  utils::write.csv(roc_df, file.path(res, "time_roc.csv"), row.names = FALSE)

  # ---- 校准：每个队列取它的中位分位时间点 -------------------------------
  cal_all <- list(); cal_notes <- list()
  for (s in sets) {
    d <- risk_df[risk_df$set == s, , drop = FALSE]
    d <- d[is.finite(d$time) & d$time > 0 & !is.na(d$event) & is.finite(d$risk), ,
           drop = FALSE]
    sub <- roc_df[roc_df$set == s, , drop = FALSE]
    if (nrow(d) < 20L || nrow(sub) == 0L) {
      cal_notes[[s]] <- "样本或时间点不足"
      next
    }
    h <- sub$horizon[[which.min(abs(sub$horizon - stats::median(sub$horizon)))]]
    c1 <- calibrate_set(d, s, h)
    cal_notes[[s]] <- c1$note
    if (length(c1$rows) > 0L) cal_all <- c(cal_all, c1$rows)
    log_info(sprintf("%s: 校准用时间点 %g，%d 组", s, h, length(c1$rows)))
  }

  cal_df <- NULL
  if (length(cal_all) > 0L) {
    cal_df <- do.call(rbind, lapply(cal_all, function(r) {
      data.frame(set = r$set, horizon = r$horizon, group = r$group, n = r$n,
                 predicted = r$predicted, observed = r$observed,
                 n_events = r$n_events, stringsAsFactors = FALSE)
    }))
    utils::write.csv(cal_df, file.path(res, "calibration.csv"), row.names = FALSE)
  }

  # ---- rms 乐观校正 -----------------------------------------------------
  rms_rows <- list(); rms_notes <- list()
  for (s in sets) {
    d <- risk_df[risk_df$set == s, , drop = FALSE]
    d <- d[is.finite(d$time) & d$time > 0 & !is.na(d$event) & is.finite(d$risk), ,
           drop = FALSE]
    set.seed(cfg$analysis$seed)          # 自助法用 RNG，紧挨着设
    rv <- rms_validate_set(d, s)
    rms_notes[[s]] <- rv$note
    if (!is.null(rv$row)) {
      rms_rows[[length(rms_rows) + 1L]] <- rv$row
      log_info(sprintf("%s: Dxy 校正 %.4f（原始 %.4f）斜率校正 %.4f（原始 %.4f）",
                       s, rv$row$Dxy_corrected, rv$row$Dxy_orig,
                       rv$row$slope_corrected, rv$row$slope_orig))
    } else {
      log_warn(sprintf("%s 的 rms 校正未完成: %s", s, rv$note))
    }
  }
  if (length(rms_rows) > 0L) {
    utils::write.csv(do.call(rbind, lapply(rms_rows, function(r) {
      data.frame(set = r$set, n = r$n, n_events = r$n_events,
                 Dxy_corrected = r$Dxy_corrected, Dxy_orig = r$Dxy_orig,
                 slope_corrected = r$slope_corrected, slope_orig = r$slope_orig,
                 R2_corrected = r$R2_corrected, R2_orig = r$R2_orig, B = r$B,
                 stringsAsFactors = FALSE)
    })), file.path(res, "rms_validate.csv"), row.names = FALSE)
  }

  # ---- 图 1：时间依赖 AUC -------------------------------------------------
  fig_ok <- tryCatch({
    sets_u <- unique(roc_df$set)
    cols <- stats::setNames(pal_categorical(length(sets_u)), sets_u)
    p <- ggplot2::ggplot(roc_df,
        ggplot2::aes(x = horizon, y = auc, colour = set)) +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
                          colour = PAL$muted, linewidth = 0.35) +
      ggplot2::geom_line(linewidth = 0.6) +
      ggplot2::geom_point(size = 1.6) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_low, ymax = ci_high),
                             width = 0, linewidth = 0.35) +
      ggplot2::scale_colour_manual(values = cols) +
      ggplot2::labs(
        x = "Time since diagnosis (native units per cohort)",
        y = "Time-dependent AUC",
        colour = NULL,
        title = "Time-dependent AUC of the LASSO-Cox risk score",
        subtitle = wrap_subtitle(paste0(
          "Dashed line = no discrimination (AUC 0.5). Horizons are the 25th/50th/75th ",
          "percentiles of EVENT times within each cohort, so the two cohorts' x-axes ",
          "are not on a shared scale (primary cohort in days, validation in years). ",
          "Horizons with fewer than ", MIN_EVENTS_AT_HORIZON,
          " events remaining are not shown."), W_ONE_HALF)) +
      theme_paper() +
      ggplot2::theme(legend.position = "bottom")
    save_pdf(file.path(fig, "time_roc.pdf"), print(p),
             width = W_ONE_HALF, height = mm(80))
    TRUE
  }, error = function(e) {
    log_warn(paste("time_roc 图失败:", conditionMessage(e)))
    FALSE
  })

  # ---- 图 2：校准 ---------------------------------------------------------
  cal_ok <- FALSE
  if (!is.null(cal_df)) {
    cal_ok <- tryCatch({
      sets_u <- unique(cal_df$set)
      cols <- stats::setNames(pal_categorical(length(sets_u)), sets_u)
      p <- ggplot2::ggplot(cal_df,
          ggplot2::aes(x = predicted, y = observed, colour = set)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                             colour = PAL$muted, linewidth = 0.35) +
        ggplot2::geom_line(linewidth = 0.6) +
        ggplot2::geom_point(ggplot2::aes(size = n)) +
        ggplot2::scale_colour_manual(values = cols) +
        ggplot2::scale_size_continuous(name = "n per group", range = c(1.2, 3.6)) +
        ggplot2::labs(
          x = "Mean predicted survival (Cox baseline)",
          y = "Observed survival (Kaplan-Meier)",
          colour = NULL,
          title = "Calibration of the LASSO-Cox risk score",
          subtitle = wrap_subtitle(paste0(
            "Dashed line = perfect calibration. Groups are quartiles of predicted ",
            "risk. Observed values are Kaplan-Meier estimates at the cohort's median ",
            "reported horizon, not the crude event proportion - censored patients ",
            "are not 'event-free'."), W_ONE_HALF)) +
        theme_paper() +
        ggplot2::theme(legend.position = "bottom")
      save_pdf(file.path(fig, "calibration.pdf"), print(p),
               width = W_ONE_HALF, height = mm(80))
      TRUE
    }, error = function(e) {
      log_warn(paste("calibration 图失败:", conditionMessage(e)))
      FALSE
    })
  }

  status$n_time_points <- nrow(roc_df)
  status$sets <- sets
  status$roc_notes <- roc_notes
  # 有队列失败时**不能报 ok** —— 部分失败也是失败，理由要带上错误原文
  status$roc_errors <- as.list(
    roc_notes[!vapply(roc_notes, function(x) !grepl("失败:", x), logical(1))])
  if (length(status$roc_errors) > 0L) {
    status$status <- "partial_error"
    status$reason <- sprintf(
      "%d 个队列的 timeROC 失败：%s",
      length(status$roc_errors), paste(unlist(status$roc_errors), collapse = " | "))
    log_error(status$reason)
  }
  status$calibration_notes <- cal_notes
  status$rms_notes <- rms_notes
  status$n_calibration_rows <- if (is.null(cal_df)) 0L else nrow(cal_df)
  status$n_rms_rows <- length(rms_rows)
  status$figures <- c(if (fig_ok) "time_roc.png", if (cal_ok) "calibration.png")
  status$outputs <- c("time_roc.csv", if (!is.null(cal_df)) "calibration.csv",
                      if (length(rms_rows) > 0L) "rms_validate.csv")

  write_json(file.path(res, "survival_diagnostics_status.json"), status)
  log_info(sprintf("完成: %d 个时间点, %d 行校准, %d 个队列的 rms 校正",
                   nrow(roc_df), status$n_calibration_rows, length(rms_rows)))
  status
}


if (!GEO_ORCHESTRATED()) {
  cfg <- load_config(parse_args())
  # R 侧没有 set_seed() 包装函数（那是 Python 侧的东西）—— 直接 set.seed
  set.seed(cfg$analysis$seed)
  run_10_survival_diagnostics(cfg)
}
