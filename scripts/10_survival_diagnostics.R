# ============================================================================

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
# **区分度（AUC）和校准（calibration 图）是两件事，缺一不可。**
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
#   figures/01-10-01-unit1-time-roc.png / 01-10-02-unit1-calibration.png
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
  # ---- 标准误：`res$inference` 是 **list**，不是 matrix --------------------
  #
  # 实测踩过（run 35483927979）：第一版按 matrix 取（`rownames == "SE"` /
  # `inf[1, ]`），三个分支全不匹配，`se` 一直是 NULL → 所有 CI 报成 `[NA, NA]`，
  # 而 AUC 本身是对的 —— **点估计看着正常，区间全丢**。
  #
  # 读 timeROC 0.4.1 源码（`R/timeROC_3.R` 末尾）确认结构：
  #
  #   inference <- list(mat_iid_rep_2 = mat_iid_rep,      # <-> AUC_2
  #                     mat_iid_rep_1 = mat_iid_rep_star, # <-> AUC_1
  #                     vect_sd_1     = vetc_sestar,      # <-> AUC_1 的 SE
  #                     vect_sd_2     = vetc_se,          # <-> AUC_2 的 SE
  #                     vect_iid_comp_time = ...)
  #
  # 无竞争风险时（本仓库 `cause = 1` 且只有一种事件）返回 `ipcwsurvivalROC`，
  # 此时 `AUC = AUC_1` → **SE 取 `vect_sd_1`**。这不是猜的：包自己的
  # `confint.ipcwsurvivalROC()` 第一行就是
  # `se <- object$inference$vect_sd_1[!is.na(object$AUC)]`。
  #
  # **不用 `confint()` 而是自己算点估计区间**，因为它为了**同时置信带**
  # 要跑 `n.sim = 2000` 次 `rnorm()` —— 那是随机过程，会引入一个新的
  # 需要设种子的来源（AGENTS 规则 11）。点估计区间是闭式的，
  # 用 `qnorm(0.975)` 算，这一步因此保持确定性。
  se <- NULL
  inf <- res$inference
  if (is.list(inf) && !is.null(inf$vect_sd_1)) {
    se <- as.numeric(inf$vect_sd_1)
  } else if (is.matrix(inf) && "SE" %in% rownames(inf)) {
    se <- as.numeric(inf["SE", ])
  } else if (is.numeric(inf)) {
    se <- as.numeric(inf)
  }
  z <- stats::qnorm(0.975)
  for (i in seq_along(ph$horizons)) {
    lo <- hi <- NA_real_
    if (!is.null(se) && length(se) >= i && is.finite(se[[i]])) {
      lo <- auc[[i]] - z * se[[i]]
      hi <- auc[[i]] + z * se[[i]]
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

    # **天和年不能画在同一根线性轴上。**
    #
    # 原来两个队列共用一个 "Time since diagnosis" 线性轴，而主队列是**天**、
    # 验证队列是**年**（规则 29 第 1 条：时间单位不跨队列共享）。实测后果：
    # GSE20685 的点全部塌缩在 x≈0（3–4 个气泡叠成一个团、误差棒互相堆叠），
    # 而训练队列的点散在 x≈430 与 x≈758，中间留下约 40% 的空白。
    # 更要紧的是**横向位置没有意义**，而图看起来像是在做跨队列的时间比较。
    #
    # 改成按队列分面、每个面板各自一根 x 轴。分析层一点没动 ——
    # 时间点本来就是各队列按自己的事件时间分位数算的（见 compute_time_roc）。
    #
    # 单位只存在于 config 的列名里（`overall survival time_days` /
    # `follow_up_duration (years)`）：risk_df 的列名被归一成
    # time/event/risk/set，原始列名已经丢了。所以单位在这里现取，
    # **不写进 time_roc.csv**（那张表已经落盘，不动它）。
    unit_of <- function(col) {
      if (is.null(col) || !nzchar(col)) return(NA_character_)
      m <- regmatches(col, regexpr("\\((days|weeks|months|years)\\)\\s*$",
                                   col, ignore.case = TRUE))
      if (length(m) == 1L) return(gsub("[()]", "", m))
      m <- regmatches(col, regexpr("_(days|weeks|months|years)\\s*$",
                                   col, ignore.case = TRUE))
      if (length(m) == 1L) return(sub("^_", "", m))
      NA_character_
    }
    val_ds <- as.character(cfg$survival$validation_dataset %||% "")
    roc_df$unit <- vapply(roc_df$set, function(s) {
      if (nzchar(val_ds) && identical(s, val_ds)) {
        unit_of(cfg$survival$validation_time_column)
      } else {
        unit_of(cfg$survival$time_column)
      }
    }, character(1))
    roc_df$panel <- ifelse(is.na(roc_df$unit), roc_df$set,
                           sprintf("%s (%s)", roc_df$set, roc_df$unit))
    roc_df$panel <- factor(roc_df$panel, levels = unique(roc_df$panel))
    if (any(is.na(roc_df$unit))) {
      log_warn("time_roc: 有队列的时间单位没能从 config 列名解析出来，面板标题只写队列名")
    }

    # **图例标签内嵌均值 AUC**（差距清单 #4，SRC-3/4 双篇惯例）:
    # 每个队列的平均 time-dependent AUC 跟在队列名后面，读者不用翻 CSV。
    auc_mean <- tapply(roc_df$auc, roc_df$set, mean, na.rm = TRUE)
    cols <- stats::setNames(
      lapply(names(cols), function(s) cols[[s]]),
      sprintf("%s (mean AUC = %.2f)", names(cols), auc_mean[names(cols)]))
    p <- ggplot2::ggplot(roc_df,
        ggplot2::aes(x = horizon, y = auc, colour = set)) +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
                          colour = PAL$muted, linewidth = 0.35) +
      ggplot2::geom_line(linewidth = 0.6) +
      ggplot2::geom_point(size = 1.6) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_low, ymax = ci_high),
                             width = 0, linewidth = 0.35) +
      ggplot2::scale_colour_manual(values = cols) +
      # free_x：每个队列一根自己的 x 轴，天和年因此永远不会落在同一尺度上
      ggplot2::facet_wrap(~ panel, scales = "free_x") +
      ggplot2::labs(
        x = "Time since diagnosis (each panel in its cohort's native unit)",
        y = "Time-dependent AUC",
        colour = NULL,
        title = "Time-dependent AUC of the LASSO-Cox risk score",
        subtitle = wrap_subtitle(paste0(
          "Dashed line = no discrimination (AUC 0.5). Horizons are the 25th/50th/75th ",
          "percentiles of EVENT times within each cohort, so each cohort is drawn on ",
          "its OWN x-axis (native unit in the panel title) — days and years are never ",
          "placed on a shared scale, and the panels must not be read across. ",
          "Horizons with fewer than ", MIN_EVENTS_AT_HORIZON,
          " events remaining are not shown."), W_DOUBLE)) +
      theme_paper() +
      ggplot2::theme(legend.position = "bottom")
    save_pdf(file.path(fig, "01-10-01-unit1-time-roc.pdf"), print(p),
             width = W_DOUBLE, height = mm(80))
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
      # **队列要能被认出来**（评审 3.8：训练队列只叫 "training"，读者不知道
      # 是哪个队列）。主队列用 dataset_id 标注，其余保留原名。
      set_lab <- stats::setNames(
        vapply(sets_u, function(s) {
          if (identical(s, "training")) sprintf("%s (training cohort)", cfg$dataset_id)
          else s
        }, character(1)), sets_u)
      cal_df$set_lab <- factor(set_lab[as.character(cal_df$set)], levels = unname(set_lab))
      p <- ggplot2::ggplot(cal_df,
          ggplot2::aes(x = predicted, y = observed, colour = set_lab)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                             colour = PAL$muted, linewidth = 0.35) +
        ggplot2::geom_line(linewidth = 0.6) +
        ggplot2::geom_point(ggplot2::aes(size = n)) +
        ggplot2::scale_colour_manual(values = stats::setNames(unname(cols), unname(set_lab))) +
        ggplot2::scale_size_continuous(name = "n per group", range = c(1.2, 3.6),
                                       breaks = pretty(cal_df$n, n = 4)) +
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
        # **图例分区处理**（评审 3.8 的歧义只在**队列颜色**那一块：
        # 2×3 网格是 size 图例的，n=30..80 是连续刻度、横排读没有歧义；
        # 队列只有 2 项，歧义出在它和 size 图例混在一张网格里）。
        # 对策：colour 图例单列放**右侧**（队列名长，竖排不挤画布）；
        # size 图例保持底部横排 3 列。第一版把两个都改成底部单列，
        # 8 行图例把 80mm 高的画布压得只剩一条缝 —— 已回退。
        ggplot2::guides(colour = ggplot2::guide_legend(order = 1),
                        size = ggplot2::guide_legend(nrow = 1, order = 2)) +
        ggplot2::theme(legend.box = "vertical",
                       legend.position = "right")
      save_pdf(file.path(fig, "01-10-02-unit1-calibration.pdf"), print(p),
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
  # 区间是怎么算的，必须落盘 —— 否则读者无法判断它可不可比
  status$ci_method <- paste0(
    "pointwise 95% CI = AUC +/- qnorm(0.975) * SE，SE 取 ",
    "timeROC 返回对象的 inference$vect_sd_1（无竞争风险时 AUC = AUC_1）。",
    "未用 confint()：它的同时置信带要 2000 次 rnorm() 模拟，会引入随机性；",
    "点估计区间是闭式的，用闭式可让这一步保持确定性。")
  status$n_ci_finite <- if (is.null(roc_df)) 0L else
    sum(is.finite(roc_df$ci_low) & is.finite(roc_df$ci_high))
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
  # ---- 7. 差距清单补图（P1-3 risk-plot / P1-4 forest / P2 nomogram+DCA）----
  # 数据源：risk_df（已读 lasso_risk_scores.csv）+ lasso_status 的 EPV 模型系数。
  # 图名走 G4 组（risk plot 三联）与 01-10 的 03/04 单元。
  rp_ok <- FALSE
  rp_notes <- list()
  if (!is.null(risk_df) && nrow(risk_df) > 0L) {
    rp_ok <- tryCatch({
      # ---- G4: risk plot 三联（训练队列）----
      tr <- risk_df[risk_df$set == "training", , drop = FALSE]
      tr <- tr[order(tr$risk), , drop = FALSE]
      cutoff_rp <- stats::median(tr$risk)
      tr$grp <- ifelse(tr$risk > cutoff_rp, "high", "low")
      # 层1: risk score 排序散点
      p1 <- ggplot2::ggplot(tr, ggplot2::aes(x = seq_len(nrow(tr)), y = risk,
                                             colour = grp)) +
        ggplot2::geom_point(size = 0.8, alpha = 0.8) +
        ggplot2::scale_colour_manual(values = c(low = PAL$down, high = PAL$up),
                                     name = "risk group") +
        ggplot2::geom_hline(yintercept = cutoff_rp, linetype = "dashed",
                            colour = PAL$muted, linewidth = 0.4) +
        ggplot2::labs(title = "Risk score (sorted)",
                      x = NULL, y = "risk score") +
        theme_paper(9) + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                        axis.ticks.x = ggplot2::element_blank())
      save_pdf(file.path(fig, "01-10-05-unit1-risk-scores.pdf"), print(p1),
               width = W_DOUBLE, height = mm(55))
      # 层2: survival time 散点（y=time, 颜色=event）
      p2 <- ggplot2::ggplot(tr, ggplot2::aes(x = seq_len(nrow(tr)), y = time,
                                             colour = factor(event))) +
        ggplot2::geom_point(size = 0.8, alpha = 0.8) +
        ggplot2::scale_colour_manual(values = c("0" = PAL$down, "1" = PAL$up),
                                     labels = c("0" = "censored", "1" = "event"),
                                     name = "status") +
        ggplot2::labs(title = "Survival time (sorted by risk)",
                      x = "samples (sorted by risk)", y = sprintf("time (%s)",
                      cfg$survival$time_column)) +
        theme_paper(9)
      save_pdf(file.path(fig, "01-10-05-unit2-survival-time.pdf"), print(p2),
               width = W_DOUBLE, height = mm(55))
      # 层3: 签名基因热图（读 07 落盘的 signature_expr.csv —— 10 不碰表达矩阵本体）
      sig_file <- file.path(res, "signature_expr.csv")
      if (file.exists(sig_file)) {
        m_raw <- utils::read.csv(sig_file, row.names = 1)
        common_gsm <- intersect(colnames(m_raw), tr$gsm)
        if (ncol(m_raw) >= 2L && length(common_gsm) >= 2L) {
          m <- t(scale(t(m_raw[, common_gsm, drop = FALSE])))
          # 列按 risk 排序（与层1/层2 同序，三联对齐）
          m <- m[, order(match(colnames(m), tr$gsm)), drop = FALSE]
          ph <- pheatmap::pheatmap(m, cluster_rows = TRUE, cluster_cols = FALSE,
                                   scale = "none", border_color = NA,
                                   fontsize = 6, legend = TRUE, silent = TRUE,
                                   labels_col = rep("", ncol(m)))  # 样本名在图上无意义（GSM 编号），与三联的隐藏 x 轴一致
          save_pdf(file.path(fig, "01-10-05-unit3-signature-heatmap.pdf"),
                   grid::grid.draw(ph$gtable), width = W_DOUBLE, height = mm(80))
        } else {
          rp_notes <- c(rp_notes, "signature_expr.csv 列不足 —— 层3 热图跳过")
        }
      } else {
        rp_notes <- c(rp_notes, "signature_expr.csv 不存在 —— 层3 热图跳过")
      }
      TRUE
    }, error = function(e) {
      rp_notes <<- c(rp_notes, sprintf("risk plot 失败: %s", conditionMessage(e)))
      FALSE
    })
  } else {
    rp_notes <- c(rp_notes, "risk_df 为空 —— G4 组整体跳过")
  }
  # ---- 8. forest plot（差距清单 #10，SRC-2 Fig5A/B 惯例）------------------
  # 单因素 Cox（risk + 每个 EPV 签名基因逐个）HR(95%CI) 的左表右图合一。
  # 数据：risk_df（risk 的 HR）+ signature_expr.csv（逐基因 HR）。
  fr_ok <- FALSE
  if (!is.null(risk_df) && nrow(risk_df) > 0L) {
    fr_ok <- tryCatch({
      tr <- risk_df[risk_df$set == "training", , drop = FALSE]
      rows <- list()
      # risk 分数自身
      fit <- survival::coxph(survival::Surv(time, event) ~ risk, data = tr)
      s <- summary(fit)
      ci <- s$conf.int
      rows[[length(rows) + 1L]] <- data.frame(
        term = "risk score", hr = unname(ci[1, "exp(coef)"]),
        lo = unname(ci[1, "lower .95"]), hi = unname(ci[1, "upper .95"]),
        p = s$logtest["pvalue"], stringsAsFactors = FALSE)
      # 逐基因（若有 signature_expr.csv）
      sig_file <- file.path(res, "signature_expr.csv")
      if (file.exists(sig_file)) {
        m_raw <- utils::read.csv(sig_file, row.names = 1)
        for (g in colnames(m_raw)) {
          d2 <- data.frame(time = tr$time, event = tr$event,
                           z = as.numeric(m_raw[, g])[match(tr$gsm, colnames(m_raw))])
          d2 <- d2[is.finite(d2$z), , drop = FALSE]
          if (nrow(d2) < 10L) next
          fit_g <- survival::coxph(survival::Surv(time, event) ~ z, data = d2)
          sg <- summary(fit_g)
          cig <- sg$conf.int
          pg <- tryCatch(sg$logtest["pvalue"], error = function(e) NA_real_)
          rows[[length(rows) + 1L]] <- data.frame(
            term = g, hr = unname(cig[1, "exp(coef)"]),
            lo = unname(cig[1, "lower .95"]), hi = unname(cig[1, "upper .95"]),
            p = pg, stringsAsFactors = FALSE)
        }
      }
      fr <- do.call(rbind, rows)
      utils::write.csv(fr, file.path(res, "cox_univariate_hr.csv"), row.names = FALSE)
      fr$label <- sprintf("%s  HR=%.2f (%.2f-%.2f)", fr$term, fr$hr, fr$lo, fr$hi)
      fr <- fr[order(fr$hr), , drop = FALSE]
      fr$label <- factor(fr$label, levels = fr$label)
      p_f <- ggplot2::ggplot(fr, ggplot2::aes(x = hr, y = label)) +
        ggplot2::geom_point(size = 1.8, colour = PAL$primary) +
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi),
                                height = 0.2, colour = PAL$primary, linewidth = 0.4) +
        ggplot2::geom_vline(xintercept = 1, linetype = "dashed",
                            colour = PAL$muted, linewidth = 0.4) +
        ggplot2::scale_x_log10() +
        ggplot2::labs(title = "Univariate Cox per-gene forest (training)",
                      subtitle = wrap_subtitle(paste0(
                        "HR per +1 SD of expression (risk score row: per +1 risk). ",
                        "Dashed line = HR 1 (no effect); x axis log scale. ",
                        "Full table: cox_univariate_hr.csv"), W_DOUBLE),
                      x = "HR (log scale)", y = NULL) +
        theme_paper(9)
      save_pdf(file.path(fig, "01-10-06-unit1-forest-plot.pdf"), print(p_f),
               width = W_DOUBLE, height = mm(90))
      TRUE
    }, error = function(e) {
      log_warn(paste("forest plot 失败:", conditionMessage(e)))
      FALSE
    })
  }

  # ---- 9. nomogram + DCA（差距清单 #7/#8，三篇文献的标配组合）------------
  nm_ok <- FALSE
  nm_note <- ""
  if (!is.null(risk_df) && nrow(risk_df) > 0L) {
    nm_ok <- tryCatch({
      tr <- risk_df[risk_df$set == "training", , drop = FALSE]
      # nomogram 构造+绘制整体包独立 tryCatch —— rms 4.x 对单变量 cph 的
      # nomogram/survest 有 "x/y lengths differ" 兼容问题；失败记原因，
      # **DCA 不受牵连**（上一版它在同一 tryCatch 里被连带跳过）。
      nom_result <- tryCatch({
        # 只用 risk + config 里的事件时间列做 3 年生存预测（时间点取中位随访）
        ddat <- data.frame(time = tr$time, event = tr$event, risk = tr$risk)
        ddat <- ddat[is.finite(ddat$time) & is.finite(ddat$risk) & !is.na(ddat$event), , drop = FALSE]
        ddat$event <- as.integer(ddat$event)
        dd <- rms::datadist(ddat); options(datadist = dd)
        on.exit(options(datadist = NULL), add = TRUE)
        surv_h <- cfg$survival
        # 3 年点：如果时间单位是年取 3，天取 1095
        is_years <- grepl("year", surv_h$time_column, ignore.case = TRUE)
        t_pred <- if (is_years) 3 else 1095
        fit_n <- rms::cph(rms::Surv(time, event) ~ risk, data = ddat,
                          surv = TRUE, x = TRUE, y = TRUE, time.inc = t_pred)
        surv_prob <- rms::survest(fit_n, times = t_pred)
        if (is.null(surv_prob)) {
          nm_note <- "rms::survest 返回空 —— nomogram 跳过（surv=TRUE 但基线不可估）"
          log_warn(nm_note)
        } else {
          nom <- rms::nomogram(fit_n, fun = function(x) surv_prob,
                               funlabel = sprintf("%d-year survival probability", t_pred),
                               fun.at = c(0.9, 0.7, 0.5, 0.3, 0.1))
          # nomogram 的 base plot 在 rms 4.x 对单变量模型有 "x/y lengths differ"
          # 的已知问题 —— 画图失败不影响 nomogram 对象；包独立 tryCatch，
          # 失败时记原因，DCA 不受牵连。
          plot_ok <- tryCatch({
            pdf(file.path(fig, "01-10-03-unit1-nomogram.pdf"),
                width = 8.5, height = 5.5)
            plot(nom, xfrac = 0.35, cex.axis = 0.7, cex.var = 0.8)
            dev.off()
            png(file.path(fig, "01-10-03-unit1-nomogram.png"),
                width = 8.5, height = 5.5, units = "in", res = 300)
            plot(nom, xfrac = 0.35, cex.axis = 0.7, cex.var = 0.8)
            dev.off()
            TRUE
          }, error = function(e) {
            try(grDevices::dev.off(), silent = TRUE)
            nm_note <<- sprintf("nomogram 绘制失败 (rms 兼容性): %s",
                                conditionMessage(e))
            FALSE
          })
        }
        list(ok = plot_ok)
      }, error = function(e) {
        nm_note <<- sprintf("nomogram 失败 (rms 兼容性): %s", conditionMessage(e))
        list(ok = FALSE)
      })
      # ---- DCA：手写净获益（不引 dcurves 新包）----
      # 判据来自 vickers 2006：NB = TP/n - FP/n * (pt/(1-pt))，
      # T=|risk 高于阈值| 的人数，TP = 其中发生事件的。
      tr$grp <- tr$risk > stats::median(tr$risk)
      pt_seq <- seq(0.05, 0.95, by = 0.05)
      # 预测概率：用 cox 基线生存 + risk 的单调映射 —— 简化：直接用 risk 分数的秩/最大秩
      # （"阈值决策"只依赖排序，秩 = risk 即可）
      rr <- tr$risk
      ev <- tr$event
      n <- length(rr)
      nb_risk <- sapply(pt_seq, function(pt) {
        hi <- rr >= stats::quantile(rr, pt)
        TP <- sum(hi & ev == 1); FP <- sum(hi & ev == 0)
        TP / n - FP / n * (pt / (1 - pt))
      })
      nb_all <- sapply(pt_seq, function(pt) mean(ev == 1) - (1 - mean(ev == 1)) * (pt / (1 - pt)))
      dca_df <- rbind(
        data.frame(pt = pt_seq, nb = nb_risk, model = "risk score"),
        data.frame(pt = pt_seq, nb = nb_all, model = "treat all"))
      utils::write.csv(dca_df, file.path(res, "dca_net_benefit.csv"), row.names = FALSE)
      p_d <- ggplot2::ggplot(dca_df, ggplot2::aes(x = pt, y = nb, colour = model)) +
        ggplot2::geom_line(linewidth = 0.6) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                            colour = PAL$muted, linewidth = 0.4) +
        ggplot2::labs(title = "Decision curve analysis (training)",
                      subtitle = wrap_subtitle(paste0(
                        "Net benefit = TP/n - FP/n * (pt/(1-pt)). ",
                        "'treat all' = intervene on everyone. ",
                        "Threshold probability is on the RISK-SCORE quantile scale ",
                        "(not a calibrated probability) — curves show relative benefit ",
                        "of using the risk score vs treating everyone."), W_ONE_HALF),
                      x = "threshold probability (risk quantile)",
                      y = "net benefit", colour = NULL) +
        theme_paper(9) + ggplot2::theme(legend.position = "bottom")
      save_pdf(file.path(fig, "01-10-04-unit1-dca.pdf"), print(p_d),
               width = W_ONE_HALF, height = mm(72))
      TRUE
    }, error = function(e) {
      nm_note <<- sprintf("nomogram/DCA 失败: %s", conditionMessage(e))
      log_warn(nm_note)
      FALSE
    })
  }
  status$figures <- c(
    if (fig_ok) c("01-10-01-unit1-time-roc.png", "01-10-01-unit2-time-roc.png"),
    if (cal_ok) "01-10-02-unit1-calibration.png",
    if (rp_ok) c("01-10-05-unit1-risk-scores.png", "01-10-05-unit2-survival-time.png",
                  "01-10-05-unit3-signature-heatmap.png"),
    if (fr_ok) "01-10-06-unit1-forest-plot.png",
    if (nm_ok && isTRUE(nom_result$ok)) "01-10-03-unit1-nomogram.png",
    if (nm_ok) "01-10-04-unit1-dca.png")
  status$risk_plot_notes <- rp_notes
  status$nomogram_dca_note <- nm_note
  status$outputs <- c("time_roc.csv", if (!is.null(cal_df)) "calibration.csv",
                      if (length(rms_rows) > 0L) "rms_validate.csv",
                      "cox_univariate_hr.csv", "dca_net_benefit.csv")

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
