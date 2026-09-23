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
#'
#' **不再有默认配置。** 仓库现在有多个数据集各自的配置（`assets/config.<GSE>.yml`），
#' 静默默认到其中某一个正是"跑错数据集"的来源 —— 产物目录、分组、阈值全都不同，
#' 而日志开头那行 dataset 很容易被跳过。所以不给默认值，让调用方显式指定。
#'
#' @return 配置文件路径
parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  cfg <- NULL
  i <- 1L
  while (i <= length(argv)) {
    if (argv[i] %in% c("--config", "-c")) {
      if (i == length(argv)) stop("--config 需要一个路径参数")
      cfg <- argv[i + 1L]
      i <- i + 2L
    } else if (argv[i] %in% c("--help", "-h")) {
      cat("用法: Rscript <script>.R --config assets/config.<GSE>.yml\n")
      quit(save = "no", status = 0L)
    } else {
      # 允许直接传位置参数作为配置路径
      cfg <- argv[i]
      i <- i + 1L
    }
  }
  if (is.null(cfg)) {
    avail <- list.files("assets", pattern = "^config\\..*\\.ya?ml$")
    stop(sprintf(
      "没有指定配置文件。用法: Rscript <script>.R --config assets/config.<GSE>.yml\n  可用配置: %s",
      if (length(avail) == 0L) "（assets/ 下没有找到 config.*.yml）"
      else paste0("\n    ", paste(avail, collapse = "\n    "))))
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

  # ---- design_mode ---------------------------------------------------------
  # 两套**不同**的门禁，不是一个门禁加个阈值：
  #
  #   small_sample（默认）探索性小样本。n < 10。功效不足是预期内的，
  #               FDR 显著基因可能为 0，走 ranked_fallback 降级路径。
  #               WGCNA / LASSO 在这个尺度上做不了，门禁直接拒绝。
  #   cohort     队列级。n >= 15（WGCNA 的通行下限）。要求每个组
  #               至少有 COHORT_MIN_PER_GROUP 个样本，否则组间比较没有意义。
  #
  # **不要为了"让某个数据集跑通"把 small_sample 的上限调大。**
  # 那是把两种设计混成一个门禁，契约会变得含糊。要跑队列就显式写 cohort。
  cfg$design_mode <- cfg$design_mode %||% "small_sample"
  if (!cfg$design_mode %in% c("small_sample", "cohort")) {
    stop(sprintf("design_mode 只能是 small_sample 或 cohort，收到: %s", cfg$design_mode))
  }

  # ---- 输出目录按数据集分目录 ---------------------------------------------
  # results/<GSE>/ 与 data/<GSE>/。
  #
  # 原来所有产物平铺在 results/ 下，**跑第二个数据集会直接覆盖第一个**。
  # 放在这里派生而不是让每个脚本自己拼路径：全部 9 个脚本都读
  # cfg$output$results_dir，改这一处就全都跟着走，不会有漏改的。
  # 需要旧行为时在配置里显式写 output.results_dir 即可覆盖。
  cfg$output             <- utils::modifyList(
    list(results_dir = file.path("results", cfg$dataset_id),
         data_dir    = file.path("data", cfg$dataset_id)),
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

# ---- 运行清单（模块零规范层）-------------------------------------------------
#
# 规范来源：用户整合文档「模块零：语言与运行时规范」。
#   §0.2 跨语言接口 —— 只走 CSV；每次转换记录维度/metadata/丢失字段
#   §0.3 版本记录   —— sessionInfo()/installed.packages() 全量 + 关键工具单列
#   §0.3 随机种子   —— 所有随机过程固定种子并记录
#   §0.4 运行日志   —— 输入数据哈希、软件版本、关键参数、决策链、
#                      人工干预记录、跨语言转换记录
#
# 产物：results/<GSE>/run_manifest.json
#
# **为什么不塞进 state.json：** state.json 记的是"这一步跑没跑成"，每步重写；
# manifest 记的是"本轮是在什么条件下跑出来的"，是证据，写入后不该再变。
# 混在一起会让后者被前者覆盖。
#
# **诚实性要求（AGENTS.md 规则 24）：** 没做的分析、没装的工具、没确认的
# 复核节点，都要在 manifest 里留下痕迹，不能因为"不影响结论"就不写。

MANIFEST_NAME <- "run_manifest.json"

# 文档 §1 点名的 R 包。**没装的记 NA（JSON null），不省略键** ——
# 键消失和"值是 null"看起来完全不同，后者才说明"本该有但没装"。
KEY_PACKAGES <- c(
  "GEOquery", "limma", "WGCNA", "clusterProfiler", "GSVA", "glmnet",
  "survival", "survminer", "timeROC", "rms", "STRINGdb",
  # §1.5 转录因子调控（文档点名 TRRUST / ChEA3，本仓库未接入）
  "TRRUST", "ChEA3",
  # §1.7 / §1.8 虚拟扰动框架（文档标主语言 Python，本仓库无）
  "scTenifoldKnk", "PerturbNet", "RegVelo"
)

# ---- 文档点名、但本仓库用不了的工具（缺口登记）-----------------------------
#
# **`KEY_PACKAGES` 回答"装没装"，这张表回答"为什么"。** 两者不能互相替代：
# `"scTenifoldKnk": null` 说明"查过了，没装"，但**没说是"装不上"还是"不归本仓库管"**
# —— 而这两件事的后续动作完全不同（前者等上游，后者去 Part 2 找）。
#
# 姊妹项目 Python 侧（`scrna` / `spatial`）有同名的 `NAMED_TOOLS` /
# `probe_named_tools()` / `named_tools_note()`，三部分的清单因此可以并排读。
#
# `kind` 的取值在三个仓库里同义：
#   r_package      R/Bioconductor 包
#   not_on_cran    CRAN / Bioconductor 上都没有（只能走数据或 API）
#   web_service    是 web 服务，没有本地包
#   python_part    规范把这一节划给 Part 2（Python），本仓库只留交接
NAMED_TOOLS <- list(
  "TRRUST" = list(
    kind = "not_on_cran", section = "§1.5",
    reason = paste0(
      "CRAN / Bioconductor 上都没有（实测 CRAN 上 `TRRUST` / `trrust` 两个名字全无）。",
      "官方分发方式是 TSV 文件，所以**本仓库真的用上了** —— 走 base R ",
      "`download.file` 取官方 TSV（不引新依赖），缓存到 `data/<GSE>/trrust_cache/` ",
      "并记 sha256。**它是 dorothea 的独立交叉验证，不替换主来源**（硬性规则 27）。")
  ),
  "ChEA3" = list(
    kind = "web_service", section = "§1.5",
    reason = paste0(
      "Ma'ayan Lab 的 web 服务，**没有 CRAN / Bioconductor 包**",
      "（实测 CRAN 上 `ChEA3` / `chea3` 两个名字全无）。调用要 `httr` / `curl`，",
      "而 CI 的 R 包列表里没有 —— 加它等于为了一个交叉验证拉一条 HTTP 依赖链。",
      "**未接入，理由同时写进 `tf_status.json` 的 limitations。**")
  ),
  "scTenifoldKnk" = list(
    kind = "python_part", section = "§1.7",
    reason = paste0(
      "规范把 §1.7 虚拟敲除标为**保留框架、主语言 Python**，落地在 ",
      "`scrna-pipeline-skill/scripts/08_virtual_perturbation.py`。",
      "本仓库（Part 1）的职责是**产出候选靶基因表**并走 CSV 交接（§0.2），",
      "不重复实现一遍 —— 两个仓库各写一套会让\"哪份是结论\"变含糊。")
  ),
  "PerturbNet" = list(
    kind = "python_part", section = "§1.8",
    reason = paste0(
      "PyPI 上的 Python 包，且 §1.8 过表达标为**保留框架、主语言 Python**，",
      "落地在 Part 2 的 `08_virtual_perturbation.py`。",
      "本仓库只做 CSV 交接。")
  ),
  "RegVelo" = list(
    kind = "python_part", section = "§1.8",
    reason = paste0(
      "Python 包（RNA velocity 方向），§1.8 主语言是 Python，落地在 Part 2。",
      "**在 Part 2 那边它也没跑成** —— 需要 spliced/unspliced 层，本数据没有；",
      "理由记在 Part 2 的状态文件里，不在这里重复。")
  )
)

#' 把 NAMED_TOOLS 整理成可写进清单的登记表。
#'
#' **不尝试 load** —— 这一节的结论是"不是 R 包 / 不归本仓库管"，去
#' `requireNamespace` 只会反复报同一个 FALSE。真正的可用性判断在
#' `requireNamespace`，这里只做一次，用来区分"登记说用不了，但环境里其实有"
#' （那说明登记过期了，值得报出来）。
#'
#' `used` **不由这里决定** —— 它取决于本轮真的跑没跑（TRRUST 就是这样：
#' 不是 R 包，但真的用上了）。所以这里只给 `available/kind/section/reason`，
#' `used` 由调用方按已落盘的状态文件填。
probe_named_tools <- function(only = NULL) {
  out <- list()
  for (tool in names(NAMED_TOOLS)) {
    if (!is.null(only) && !(tool %in% only)) next
    meta <- NAMED_TOOLS[[tool]]
    avail <- isTRUE(requireNamespace(tool, quietly = TRUE))
    if (avail) {
      log_warn(sprintf(
        "工具 %s 登记为不可用，但环境里能 requireNamespace —— 登记需要更新", tool))
    }
    out[[tool]] <- list(available = avail, kind = meta$kind,
                        section = meta$section, reason = meta$reason)
  }
  out
}

#' 一句话说明文档点名的工具里哪些没用上、为什么。
named_tools_note <- function() {
  paste0(
    "文档 §1 点名的工具里，**§1.7 / §1.8 的虚拟扰动三个工具不归本仓库管** —— ",
    "规范把这两节标为主语言 Python，落地在 `scrna-pipeline-skill`；",
    "本仓库的职责是产出候选靶基因表并走 CSV 交接（§0.2）。",
    "§1.5 的 TRRUST **真的用上了**（官方 TSV，作为 dorothea 的独立交叉验证），",
    "ChEA3 未接入（web 服务，无本地包，调用需 httr/curl）。",
    "逐条理由见清单的 `named_tools` 字段。")
}

#' 把缺口登记写进清单。
#'
#' **`used` 由调用方传进来**（见 `probe_named_tools()` 的说明）——
#' 在这里自己判一遍就会和真正的执行结果分叉。
record_named_tools <- function(cfg, used = character(0)) {
  tools <- probe_named_tools()
  for (t in names(tools)) {
    tools[[t]]$used <- t %in% used
    # 既没说 used 也没写理由 = 没说清。理由一定非空，所以这里只补一个断言
    if (!tools[[t]]$used && !nzchar(tools[[t]]$reason %||% "")) {
      tools[[t]]$reason <- "**理由缺失**（登记不完整，必须补上）"
    }
  }
  m <- read_manifest(cfg)
  m$named_tools <- tools
  m$named_tools_note <- named_tools_note()
  m$dataset_id <- cfg$dataset_id
  m$updated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  write_manifest(cfg, m)
  invisible(tools)
}

manifest_path <- function(cfg) file.path(cfg$output$results_dir, MANIFEST_NAME)

#' 读清单。文件不存在返回空 list（不是 NULL，便于直接取字段）
read_manifest <- function(cfg) {
  p <- manifest_path(cfg)
  if (!file.exists(p)) return(list())
  tryCatch(jsonlite::fromJSON(p, simplifyVector = FALSE), error = function(e) list())
}

#' 清单专用写出：**NA 映射成 JSON null**，与 Python 侧的 None 对齐。
#' 共用的 write_json() 不带 na=，会把 NA 写成字符串 "NA" —— 那是两种
#' 不同的东西：null 是"没有这个值"，"NA" 是"值是字符串 NA"。
write_manifest <- function(cfg, m) {
  p <- manifest_path(cfg)
  tmp <- paste0(p, ".tmp")
  writeLines(jsonlite::toJSON(m, auto_unbox = TRUE, pretty = TRUE,
                              null = "null", na = "null"), tmp)
  if (file.exists(p)) unlink(p)
  file.rename(tmp, p)
  invisible(m)
}

#' 建立本轮清单骨架。**会清掉上一轮的内容** —— 清单描述的是本轮。
init_manifest <- function(cfg, language = "R") {
  m <- list(
    dataset_id = cfg$dataset_id,
    language = language,
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    seed = cfg$analysis$seed,
    versions = list(), key_versions = list(), inputs = list(),
    params = list(), decisions = list(), human_review = list(),
    cross_language = list()
  )
  write_manifest(cfg, m)
  invisible(m)
}

manifest_append <- function(cfg, key, entry) {
  m <- read_manifest(cfg)
  cur <- m[[key]]
  if (is.null(cur)) cur <- list()
  cur[[length(cur) + 1L]] <- entry
  m[[key]] <- cur
  m$dataset_id <- cfg$dataset_id
  m$updated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  write_manifest(cfg, m)
  invisible(m)
}

#' §0.4 文件哈希。
#'
#' 优先 digest::digest(algo = "sha256")；没装 digest 时退回 tools::md5sum。
#' **实际用的算法写进 hash_algo 字段** —— 用了哪种算法不能靠猜，
#' 换算法时前后两轮的哈希不可比。
file_hash <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (requireNamespace("digest", quietly = TRUE)) {
    return(list(algo = "sha256",
                value = digest::digest(path, algo = "sha256", file = TRUE)))
  }
  list(algo = "md5", value = unname(tools::md5sum(path)))
}

#' §0.4 输入数据哈希。文件不存在时**记 missing 而不是报错**
record_input <- function(cfg, path, label = NULL, required = TRUE) {
  h <- file_hash(path)
  entry <- list(
    label = label %||% basename(path), path = path, required = isTRUE(required),
    status = if (is.null(h)) "missing" else "present"
  )
  if (!is.null(h)) {
    entry$hash_algo <- h$algo
    entry$hash <- h$value
    entry$bytes <- unname(file.size(path))
  }
  manifest_append(cfg, "inputs", entry)
  invisible(entry)
}

#' §0.3 版本记录：全量已安装包 + 关键工具单独记版本。
#'
#' 脚本一律 `pkg::fun()` 写全名、不 attach，所以 `sessionInfo()$otherPkgs`
#' 是空的 —— **"没 attach"不等于"没用"**。用 installed.packages() 才是
#' 实际可用的全集。
capture_versions <- function(cfg, key_packages = KEY_PACKAGES) {
  full <- list()
  ip <- tryCatch(
    utils::installed.packages()[, c("Package", "Version"), drop = FALSE],
    error = function(e) NULL
  )
  if (!is.null(ip)) {
    # **必须 unname()。** `ip[i, "Version"]` 返回的是**带名字**的长度 1 字符
    # 向量（名字是 "Version"），而 jsonlite 对"有名字的原子向量"序列化成
    # 对象 —— 会写出 `{"dplyr": {"Version": "1.1.4"}}` 而不是
    # `{"dplyr": "1.1.4"}`。Python 侧写的是后者，两边 schema 就对不上了，
    # 而清单恰恰是要并排读的。这个错只有 CI 里看到 JSON 才发现。
    for (i in seq_len(nrow(ip))) {
      full[[ip[i, "Package"]]] <- unname(ip[i, "Version"])
    }
    full <- full[order(tolower(names(full)))]
  } else {
    log_warn("installed.packages() 失败 —— versions 会不完整")
  }

  key <- stats::setNames(vector("list", length(key_packages)), key_packages)
  for (p in key_packages) {
    v <- full[[p]]
    key[[p]] <- if (is.null(v)) NA_character_ else as.character(v)
  }

  m <- read_manifest(cfg)
  m$versions <- full
  m$key_versions <- key
  m$n_packages <- length(full)
  m$r_version <- paste(R.version$major, R.version$minor, sep = ".")
  m$platform <- R.version$platform
  m$dataset_id <- cfg$dataset_id
  m$updated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  write_manifest(cfg, m)

  missing <- key_packages[vapply(key_packages, function(p) is.na(key[[p]]), logical(1))]
  if (length(missing) > 0) {
    log_warn(sprintf("关键工具未安装（%d/%d）：%s",
                     length(missing), length(key_packages), paste(missing, collapse = ", ")))
  } else {
    log_info(sprintf("关键工具全部就位（%d 个），共记录 %d 个已安装包",
                     length(key_packages), length(full)))
  }
  invisible(key)
}

#' §0.4 关键参数完整记录（含随机种子）
record_params <- function(cfg, params) {
  m <- read_manifest(cfg)
  p <- m$params
  if (is.null(p)) p <- list()
  for (nm in names(params)) p[[nm]] <- params[[nm]]
  m$params <- p
  # **不能写 `m$seed <- cfg$analysis$seed`。** R 里 `x$k <- NULL` 是**删键**
  # 不是"设成空值" —— seed 配错/缺失时种子会从清单里静默消失，而
  # §0.3 要求所有随机过程都要记种子。用 list 赋值保住键，值可以是 null。
  m["seed"] <- list(cfg$analysis$seed)
  m$dataset_id <- cfg$dataset_id
  m$updated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  write_manifest(cfg, m)
  invisible(m)
}

#' §0.4 Agent 决策链：从原始问题到最终结论的每一步推理。
#' evidence 要写**支持这个选择的实际数字**，不是"因为这是通行做法"。
record_decision <- function(cfg, node, question, answer, evidence = "") {
  manifest_append(cfg, "decisions", list(
    node = node, question = question, answer = answer, evidence = evidence,
    at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  ))
}

#' §0.4 人工干预记录。
#'
#' status：pending（需确认，未确认）/ confirmed / overridden（人推翻了自动
#' 结果，note 写改成什么）/ not_needed（本数据集不涉及）。
#'
#' **默认 pending 而不是 confirmed。** 自动化流水线不能替人签字 ——
#' 把未确认的节点默认记成已确认，等于把复核节点变成摆设。
record_human_review <- function(cfg, node, required = TRUE,
                                status = "pending", note = "") {
  manifest_append(cfg, "human_review", list(
    node = node, required = isTRUE(required), status = status, note = note,
    at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  ))
}

#' §0.2 跨语言转换记录。
#'
#' 文档要求记录转换前后维度、metadata 字段数、丢失字段清单。桥接工具限定
#' zellkonverter / anndata2ri，**禁止 sceasy**（维护状态差、metadata 丢失
#' 风险高）。
#'
#' 本仓库与姊妹仓库之间只走 CSV，所以正常路径下 before/after 是行列数与
#' 列名集合；真正发生对象级转换时才填 tool。
record_cross_language <- function(cfg, src, dst, format,
                                  before = NULL, after = NULL,
                                  lost = character(0), tool = "", note = "") {
  manifest_append(cfg, "cross_language", list(
    src = src, dst = dst, format = format, tool = tool,
    before = before %||% list(), after = after %||% list(),
    lost_fields = sort(as.character(lost)), note = note,
    at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  ))
}

#' 给验收用的一行摘要
#'
#' **必需项缺失和可选项缺失分开报。** `clinical.csv` 是 required = FALSE ——
#' 本来就可以没有。把它算进"缺失"会让每个没有临床表的数据集都判失败，
#' 那是把"设计如此"当成"出错了"。两者都可见，但只有必需项判失败。
manifest_summary <- function(cfg) {
  m <- read_manifest(cfg)
  if (length(m) == 0) return(list(present = FALSE))
  pick <- function(k) m[[k]] %||% list()
  ins <- pick("inputs")
  statuses <- vapply(ins, function(i) i$status %||% "?", character(1))
  labels <- vapply(ins, function(i) i$label %||% "?", character(1))
  reqd <- vapply(ins, function(i) isTRUE(i$required), logical(1))
  is_missing <- statuses == "missing"
  hstat <- vapply(pick("human_review"), function(h) h$status %||% "?", character(1))
  hnode <- vapply(pick("human_review"), function(h) h$node %||% "?", character(1))
  list(
    present = TRUE,
    n_versions = length(pick("versions")),
    n_inputs = length(ins),
    inputs_missing = sort(labels[is_missing]),
    inputs_missing_required = sort(labels[is_missing & reqd]),
    n_decisions = length(pick("decisions")),
    human_review_pending = sort(hnode[hstat == "pending"]),
    n_cross_language = length(pick("cross_language"))
  )
}

# ---- 统一调色板 -------------------------------------------------------------
#
# **语义固定，所有图共用同一套，色值全部取自 SCI 发表常用色板。**
#
# | 角色 | 来源 |
# | --- | --- |
# | 上调 / tumor、下调 / normal | **ColorBrewer RdBu** 的两个端点 —— 基因组学论文里最通用的发散尺，方向色与发散尺同源 |
# | 显著性（连续） | **viridis**（Nature Methods 推荐；感知均匀、灰度单调、色盲安全） |
# | 模块 / 多组（分类） | **Okabe-Ito** 及其加深变体，穷举筛出的 4 色 |
#
# **色值是算出来的，不是挑出来的。** 判据取自 better-colors skill：
#   * 色相相差 15° 以内视为**同一个颜色** —— 承载不同含义的颜色必须拉开 15° 以上
#   * 颜色从不是唯一的语义载体（见各图的 shape / alpha / 位置编码）
#   * 报告任何对比度之前先测量，不要估
#
# 实测暴露过的四个问题（`_pal_sci.py` 可复算）：
#   1. 旧 p 值直方图的 `firebrick` 与 up 红**只差 0.4°** —— 按判据就是同一个颜色
#   2. 旧分类色板在 protanopia 下 up 红与棕距离 **0.002**（几乎重合）
#   3. magma 序列色中段距 up 红仅 **2.9°** —— 序列色中段就是"上调红"
#   4. **viridis 的暗端 `#365C8D` 距 down 蓝只有 3.1°** ——
#      所以 viridis 必须**截去暗端**，从 `#277F8E`（距 down 41.1°）起用
PAL <- list(
  up       = "#B2182B",   # ColorBrewer RdBu 红端     h=22.4  白底 6.6:1
  down     = "#2166AC",   # ColorBrewer RdBu 蓝端     h=252.4 白底 6.4:1
  ns       = "#BDBDBD",   # 不显著（中性灰，不承载色相语义）
  # **primary 必须存在**（与 scrna/spatial 的 Python 侧同名同值）：
  # 实测多处代码写 PAL$primary，而 R 侧原先没有这个键 —— NULL 被 c() 丢掉后
  # 颜色数比名字数少一个，报 `'names' attribute [N] must be the same length...`，
  # 整图/整步失败。值取 Okabe-Ito 蓝（与 Python 侧一致）。
  primary  = "#0072B2",
  mid      = "#F7F7F7",   # 发散色中点（RdBu 的中性色）
  ink      = "#1A1A1A",   # 文字 / 参考线（白底 17.4:1）
  muted    = "#666666",   # 副标题 / 次要说明
  edge     = "#8A8A8A",   # 网络边（中性灰，不承载语义）
  grid     = "#E8E8E8"
)

#' 发散色板：**ColorBrewer RdBu**（蓝 → 近白 → 红）
#'
#' 用于 Z-score 热图。**与方向色同源** —— 热图上的红就是火山图上"上调"的那个红，
#' 不需要读者重新学一遍。RdBu 是基因组学论文里最常用的发散尺。
pal_diverging <- function(n = 100) {
  grDevices::colorRampPalette(c(PAL$down, PAL$mid, PAL$up))(n)
}

#' 序列色板：**viridis**（截去与 down 蓝撞色相的暗端）
#'
#' 完整 viridis 的暗端三档色相是 318 / 292 / **255**，其中 `#365C8D` 距
#' down 蓝（252.4）只有 **3.1°** —— 按 15° 判据就是"下调蓝"，
#' 深色点会被读成"下调"。从 `#277F8E`（h=211.3，距 down 41.1°）起截断，
#' 保留 viridis 的感知均匀性与灰度单调性（亮度 0.552 → 0.918）。
pal_sequential <- function(n = 256) {
  grDevices::colorRampPalette(
    c("#277F8E", "#1FA187", "#4AC16D", "#A0DA39", "#FDE725"))(n)
}

#' 分类色板（模块 / 多组）：穷举出的色盲安全四色
#'
#' 这四个不是挑的，是搜出来的（`_pal_sci.py`）。**不是贪心** ——
#' 贪心会选出两个绿色（h=149.9 / 165.5，只差 15.6°）。这里在可行集里
#' 穷举所有 4 色组合，取**最小两两 OKLab 距离最大**的那组（实测 0.246）。
#' 约束三条：
#'   * 与 up / down / 序列色 的色相都拉开 15° 以上
#'   * 与它们、以及彼此之间，在 protanopia / deuteranopia / tritanopia 下距离 >= 0.05
#'   * 白底对比度 >= 2.0（实心圆点仍清晰可见的底线）
#'
#' 被剔除的典型：`#0072B2`（距 down 仅 13°）、`#1B7C8C` `#8E5A9E`（撞 down）、
#' `#C9B800`（撞序列色）、`#F0E442`（白底仅 1.32:1，小圆点看不见）、
#' 灰色（灰已被 ns / other 占用）。
#'
#' **超过 4 个模块时降级**：在四色之间插值，插出来的颜色**没有**经过上面的验证。
pal_categorical <- function(k) {
  base <- c("#2E6B3E", "#B5527A", "#56B4E9", "#E69F00")   # 绿 / 玫红 / 天蓝 / 橙
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

#' 把长副标题按画布宽度折行
#'
#' **ggplot 的副标题不会自动换行。** 超过画布宽度的部分被**静默裁掉** ——
#' 不是显示成省略号，是根本没有，所以从图上完全看不出"有字没显示"。
#' 实测 PPI 的副标题 455 字符、火山图 275 字符，都远超一行能容纳的长度。
#'
#' 每行字符数按画布宽度和字号估：副标题字号是 `base_size - 1.5`，
#' 无衬线体的平均字宽约 `0.52 em`，再留 15% 安全余量。
#'
#' @param x         副标题文本
#' @param fig_width 图的宽度（英寸），要和 `save_pdf()` 传的值一致
#' @param base_size 主题基准字号
wrap_subtitle <- function(x, fig_width = 8, base_size = 10) {
  size <- base_size - 1.5
  chars <- max(30L, floor(fig_width * 72 / (size * 0.52) * 0.85))
  paste(strwrap(x, width = chars), collapse = "\n")
}

#' 把多条说明各折各的行，再拼成一段多行图注
#'
#' `wrap_subtitle()` 用的是 `strwrap()`，它把输入当成**一个段落**，
#' 里面的 `\n` 会被当普通空白折叠掉 —— 所以直接传带换行的文本进去，
#' 结构化图注会被压成一坨。这里逐条折行再拼，行结构得以保留。
#'
#' @param lines     字符向量，每条是一行图注
#' @param fig_width 图宽（英寸），要和 `save_pdf()` 传的值一致
#' @param base_size 主题基准字号
wrap_caption <- function(lines, fig_width = 8, base_size = 10) {
  lines <- lines[!is.na(lines) & nzchar(lines)]
  if (length(lines) == 0L) return("")
  paste(vapply(lines, function(l) wrap_subtitle(l, fig_width, base_size),
               character(1)), collapse = "\n")
}

# ---- 行名标签放得下吗 ------------------------------------------------------
#
# **判据是算出来的，不是看出来的。** "标签挤不挤"如果靠肉眼判断，就会变成
# 一个没人能复现的印象：同一张图有人说挤有人说不挤，而且加一个基因、
# 改一次画布高度之后没人会想起来重新看一遍。
#
# 所以把它写成算术：一行标签需要 `fontsize + min_gap` 点的垂直空间，
# 画布能给的是 `height_in * 72 * panel_frac` 点（panel_frac 扣掉标题、
# 副标题、坐标轴和图例占掉的高度）。
#
# 实测：GSE42568 的 top50 热图 50 个基因、5pt 字号、5.75in 高，
# 每行只有约 3.6px 间隙（约 1.7pt），确实挤 —— 而 `show_rownames` 原来
# 写的是硬编码的 `length(genes) <= 60`，跟画布高度毫无关系。

#' 画布高度能容纳多少行标签
#'
#' @param height_in  图高（英寸），要和 `save_pdf()` / `pheatmap()` 传的值一致
#' @param fontsize   标签字号（pt）
#' @param panel_frac 绘图面板占图高的比例。0.75 是带标题、副标题、
#'   底部横排图例时的实测经验值；`pheatmap` 没有副标题、图例是右侧细色条，
#'   用 0.82 更准。
#' @param min_gap    相邻标签之间至少要留的空白（pt）。2.5pt 约合 10px（300dpi），
#'   是"一眼能分开两行"的下限；1.5pt 已经会糊成一片。
#' @return 整数上限
label_budget <- function(height_in, fontsize, panel_frac = 0.75, min_gap = 2.5) {
  avail <- height_in * 72 * panel_frac
  as.integer(floor(avail / (fontsize + min_gap)))
}

#' 这些标签放得下吗
#'
#' @inheritParams label_budget
#' @param n 标签行数
#' @return 逻辑值；`n` 为 0 或非有限值时返回 FALSE
fits_labels <- function(n, height_in, fontsize, panel_frac = 0.75, min_gap = 2.5) {
  if (!is.finite(n) || n <= 0L) return(FALSE)
  n <= label_budget(height_in, fontsize, panel_frac, min_gap)
}

#' 决定要不要显示行名，并把决定与依据写进日志
#'
#' 返回逻辑值，供 `show_rownames = ` 直接用。**日志里必须留下算式** ——
#' 否则"这张图为什么没有行名"又要靠猜。
#'
#' 传了 `figure` 就同时记进 `label_decisions.csv`（见
#' `write_label_decisions()`）。**日志会随 CI 日志一起滚掉，文件不会** ——
#' 想回答"这轮为什么把行名藏了"，翻文件比翻几万行日志快。
#'
#' @param what 标签指代的东西，用于日志（如 "热图基因"）
#' @param figure 图的名字（如 "01-04-01-unit1-top50-heatmap"）；`NULL` 表示只记日志不落盘
#' @return 逻辑值
decide_rownames <- function(n, height_in, fontsize, what = "行名",
                            panel_frac = 0.75, min_gap = 2.5, figure = NULL) {
  budget <- label_budget(height_in, fontsize, panel_frac, min_gap)
  ok <- fits_labels(n, height_in, fontsize, panel_frac, min_gap)
  log_info(sprintf(
    "%s标签: %d 行 / 画布可容纳 %d 行（高 %.2fin，字号 %gpt，行距余量 %gpt）-> %s",
    what, n, budget, height_in, fontsize, min_gap,
    if (ok) "显示行名" else "**不显示行名**（会糊成一片）"))
  if (!is.null(figure)) {
    record_label_decision(figure, what, n, height_in, fontsize,
                          panel_frac, min_gap, ok)
  }
  ok
}

# 本轮所有标签决策。放环境里而不是全局变量：`source()` 进来的脚本共享
# 同一个 globalenv，用 `<<-` 赋值会污染调用方；环境是显式的容器。
.label_decisions <- new.env(parent = emptyenv())
.label_decisions$rows <- list()

#' 记一条标签决策
#'
#' 落盘的理由和 `ppi_plot_layout.csv` 一样（规则 16）：**从图上反推
#' "这张为什么没行名"是猜**。记下 `n_labels` / `capacity` 和四个参数，
#' 决策就是可核对的数据而不是事后叙述。
#'
#' @inheritParams label_budget
#' @param figure 图名
#' @param what   标签种类
#' @param n      标签行数
#' @param shown  最终是否显示
record_label_decision <- function(figure, what, n, height_in, fontsize,
                                  panel_frac, min_gap, shown) {
  .label_decisions$rows[[length(.label_decisions$rows) + 1L]] <- data.frame(
    figure    = as.character(figure),
    label     = as.character(what),
    n_labels  = as.integer(n),
    capacity  = label_budget(height_in, fontsize, panel_frac, min_gap),
    height_in = as.numeric(height_in),
    fontsize  = as.numeric(fontsize),
    panel_frac = as.numeric(panel_frac),
    min_gap   = as.numeric(min_gap),
    shown     = isTRUE(shown),
    stringsAsFactors = FALSE)
  invisible(NULL)
}

#' 把本轮所有标签决策写成 CSV
#'
#' 即使一条都没有也写出带表头的空文件：**"文件不存在"和"没有需要标签的图"
#' 是两件事**，前者会让下游以为这步没跑。
write_label_decisions <- function(path) {
  rows <- .label_decisions$rows
  empty <- data.frame(figure = character(0), label = character(0),
                      n_labels = integer(0), capacity = integer(0),
                      height_in = numeric(0), fontsize = numeric(0),
                      panel_frac = numeric(0), min_gap = numeric(0),
                      shown = logical(0), stringsAsFactors = FALSE)
  df <- if (length(rows) == 0L) empty else do.call(rbind, rows)
  utils::write.csv(df, path, row.names = FALSE)
  log_info(sprintf("标签决策已落盘: %s（%d 条，其中 %d 条隐藏了行名）",
                   basename(path), nrow(df), sum(!df$shown)))
  invisible(df)
}

#' 组内协方差椭圆的坐标
#'
#' **不用 `stat_ellipse()`。** 实测它在每组 3 个样本时产出了空数据 —— 图上
#' 只有点、没有椭圆，而 ggplot 不会报错，坐标范围也没被撑大，所以从图上
#' 完全看不出"这一层没画"。自己算有三个好处：画得出来、坐标能落盘核对、
#' 半径系数是显式的而不是藏在 stat 内部。
#'
#' 半径按 `sqrt(qchisq(level, 2))`（`type = "norm"` 的口径）：
#' 把协方差当作**已知**。ggplot 默认的 `"t"` 用
#' `sqrt(2 * qf(level, 2, n-2))`，每组 3 个样本时是 6.16 倍标准差，
#' 椭圆会比数据范围大 6 倍、把点压成中心一小团，所以这里不用它。
#' 代价是**低估**了小样本下协方差本身的不确定性，这一点必须在图注里写明。
#'
#' @param x,y     同一组的点坐标
#' @param level   置信水平
#' @param n       椭圆上的点数
#' @return 数据框，列为 `x` / `y` / `radius_sd`
ellipse_points <- function(x, y, level = 0.95, n = 120L) {
  if (length(x) < 3L || stats::var(x) == 0 || stats::var(y) == 0) return(NULL)
  cov_m <- stats::cov(cbind(x, y))
  if (any(!is.finite(cov_m)) || det(cov_m) <= 0) return(NULL)
  eig <- eigen(cov_m, symmetric = TRUE)
  radius <- sqrt(stats::qchisq(level, 2))
  theta <- seq(0, 2 * pi, length.out = n)
  # 单位圆 -> 按特征值缩放 -> 旋到主轴方向
  pts <- cbind(cos(theta), sin(theta)) %*%
    diag(sqrt(pmax(eig$values, 0)), 2L) %*% t(eig$vectors)
  data.frame(
    x = mean(x) + radius * pts[, 1L],
    y = mean(y) + radius * pts[, 2L],
    radius_sd = radius
  )
}

#' 所有图共用的主题，保证字号、网格、留白一致
#'
#' **图例默认放底部、横向排列。** 右侧图例直接吃掉图的宽度 ——
#' 实测 PCA / 火山图 / dotplot / 热图都是右侧图例，横向占掉一大块，
#' 而底部横排在同样信息量下几乎不增加图幅。需要右侧的图单独覆盖。
theme_paper <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = PAL$grid, linewidth = 0.25),
      panel.border     = ggplot2::element_rect(colour = PAL$grid, linewidth = 0.4),
      plot.title       = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle    = ggplot2::element_text(colour = PAL$muted, size = base_size - 1.5),
      plot.margin      = ggplot2::margin(5, 6, 4, 5),
      # **图例一律框外右侧、纵向排列**（用户约定 v2，2026-09-23）。
      legend.position   = "right",                  # 框外右侧（不挤压主图）
      legend.justification = c(0.5, 0.5),
      # **同一个图例内部的键要竖着排**（不是横排）。
      # 实测反馈：`horizontal` 让一个图例的键横着铺开（volcano 的 direction 图例
      # 三个键横排、GSEA dotplot 的色标与 size 横排），既挤主图又难读。
      # 竖排之后图例块宽度 = 最宽的一个键，不再随键数增长。
      legend.direction  = "vertical",
      # **多个图例之间也要竖着摞**（`legend.box` 管图例之间，`legend.direction`
      # 管单个图例内部 —— 两个都要设成纵向）。
      # 实测 01-03-01-unit1-volcano-plot（colour + alpha 两个图例）：并排后总宽超出
      # 183 mm，两端被静默裁掉 —— 左端只剩 "…nificant"，右端计数整段消失。
      # ggplot 裁图例**不报警**，图上只表现为"这个图例好像短了一截"。
      legend.box        = "vertical",
      legend.title      = ggplot2::element_text(size = base_size - 1),
      legend.text       = ggplot2::element_text(size = base_size - 1.5),
  legend.key.size   = ggplot2::unit(0.7, "lines"),   # 右上角图例空间有限，键略缩
      legend.margin     = ggplot2::margin(1, 1, 1, 1),
      legend.box.spacing = ggplot2::unit(3, "pt"),
      legend.box.margin  = ggplot2::margin(0, 0, 0, 0)
    )
}

# ---- 画布尺寸（毫米）-------------------------------------------------------
#
# **为什么按毫米而不是英寸。** 期刊的栏宽是用毫米规定的，英寸是排版软件
# 内部单位。写英寸时"这图多宽"要靠换算才知道，写毫米时一眼就能对上投稿要求。
# 参考规范：K-Dense `scientific-visualization` skill。
#
# 三个标准宽度覆盖绝大多数情况；超宽图（类别数多）应**夹到 W_DOUBLE**，
# 而不是让宽度随类别数无限增长 —— 否则会画出装不进任何期刊一页的图。
# 姊妹项目 Python 侧（`scrna-pipeline-skill` / `spatial-pipeline-skill`）
# 用同一组数值，三部分文档的图幅因此可比。
W_SINGLE   <- 89 / 25.4    # 单栏
W_ONE_HALF <- 136 / 25.4   # 一栏半
W_DOUBLE   <- 183 / 25.4   # 双栏（通栏）

#' 毫米转英寸（`save_pdf()` 的宽高单位是英寸）
#' @param ... 一个或多个毫米数值
mm <- function(...) c(...) / 25.4

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
#' @param width,height 画布尺寸（英寸）。**用 `mm()` / `W_SINGLE` / `W_ONE_HALF` /
#'   `W_DOUBLE` 给值**，不要写裸英寸数字 —— 栏宽是按毫米规定的。
#' @param dpi   PNG 分辨率。**300 是投稿常规下限** —— 期刊对位图普遍要求
#'   300 dpi，线条图常要 600–1200。原来默认 150：实测 183 mm 宽的图只渲染成
#'   1080 px，放大一倍就糊，而"分辨率不够"在缩略图上完全看不出来。
#'   PDF 那一路是**矢量**的、与 dpi 无关，所以这个值只影响位图产物（PNG）。
#'   **在这里改而不是逐图传参**：分辨率是每一张图都有的属性，逐图传一定会漏
#'   （实测全仓库 12 个出图脚本、19+14 张图，没有一个调用点传过 dpi）。
save_pdf <- function(path, expr, width = W_DOUBLE, height = mm(64), dpi = 300) {
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

#' 把一组 `characteristics_ch1` 取值解析成 key -> value
#'
#' SOFT 里的形状是 `"key: value"`，一个样本有多行。
#' **同名 key 在一个样本里出现多次时用 " | " 连接** —— GEO 里确实有
#' （GSE20685 的 regimen 就有），直接取最后一个会静默丢数据。
#'
#' 没有 `"key:"` 前缀的自由文本单独归到 `unnamed`，**不伪装成某个字段**：
#' 把它按第一个词当 key 会凭空造出一个语义不明的列，下游按列名取性状时会中招。
#'
#' @param values 一个样本的全部 `Sample_characteristics_ch1` 取值
#' @return 命名 list
parse_characteristics <- function(values) {
  out <- list()
  for (v in values) {
    if (!nzchar(trimws(v))) next
    if (grepl(":", v, fixed = TRUE)) {
      k   <- trimws(sub(":.*$", "", v))
      val <- trimws(sub("^[^:]*:", "", v))
    } else {
      k <- "unnamed"; val <- trimws(v)
    }
    if (!nzchar(k)) next
    out[[k]] <- if (is.null(out[[k]])) val else paste(out[[k]], val, sep = " | ")
  }
  out
}

#' 从 SOFT 样本块建一张临床宽表
#'
#' step 00（主队列）与 step 07（外部验证队列）共用这一个实现 ——
#' 两处各写一份解析器，迟早会在"同名 key 怎么办"上分叉，
#' 而那种分叉不会报错，只会让两个队列的字段含义悄悄不一致。
#'
#' @param blocks  `split_soft_samples()` 的结果
#' @param gsm     样本 ID，顺序要与 `blocks` 一致
#' @return data.frame，第一列是 `gsm`
clinical_table <- function(blocks, gsm) {
  parsed <- lapply(blocks, function(b) parse_characteristics(soft_values(b, "Sample_characteristics_ch1")))
  keys <- unique(unlist(lapply(parsed, names)))
  out <- data.frame(gsm = gsm, stringsAsFactors = FALSE)
  for (k in keys) {
    out[[k]] <- vapply(parsed, function(p) {
      v <- p[[k]]
      if (is.null(v) || !nzchar(v)) NA_character_ else v
    }, character(1))
  }
  out
}

#' 临床字段画像：每个字段的覆盖度、类型、取值分布
#'
#' LASSO/Cox 的门禁要据此判断"有没有可用的终点、有多少个事件"。
#' **不在这里判定用哪个终点** —— 终点必须由 config 显式指定，
#' 自动配对 time/event 列在字段名不规整时会静默配错（实测 GSE20685 的
#' `event_death` 和 `follow_up_duration (years)` 名字里没有任何共同词）。
#'
#' @param clinical `clinical_table()` 的结果
#' @return 命名 list，每个字段一项
clinical_field_profile <- function(clinical) {
  out <- lapply(setdiff(colnames(clinical), "gsm"), function(k) {
    v <- clinical[[k]]
    present <- v[!is.na(v) & nzchar(v)]
    num <- suppressWarnings(as.numeric(present))
    numeric_ok <- length(present) > 0L && !anyNA(num)
    info <- list(n_present = length(present),
                 n_missing = sum(is.na(v) | !nzchar(v)),
                 n_distinct = length(unique(present)))
    if (numeric_ok) {
      info$kind <- "numeric"
      info$min <- min(num); info$max <- max(num); info$median <- stats::median(num)
    } else {
      info$kind <- "categorical"
      tab <- sort(table(present), decreasing = TRUE)
      # 取值太多就只留前 12 个，避免 JSON 里塞进 327 个样本标题之类的东西
      info$top_values <- as.list(head(as.integer(tab), 12L))
      names(info$top_values) <- names(head(tab, 12L))
      info$n_levels <- length(tab)
    }
    info
  })
  names(out) <- setdiff(colnames(clinical), "gsm")
  out
}

# ---- 基因/统计小工具 -------------------------------------------------------

#' 判断表达矩阵是 log2 尺度还是线性强度尺度
#'
#' **为什么必须查这件事。** 本流水线从 GEO **series matrix** 出发，
#' 而那是提交者放上去的东西 —— 有的已经 log2（Affymetrix 的 RMA/MAS5 输出、
#' 大多数 Illumina 输出），有的还是线性荧光强度。下游 limma 假定数据是
#' log 尺度：把线性数据直接喂进去，**差异倍数会变成"强度比的对数"而不是
#' 表达比的对数**，火山图、阈值、GSEA 排序全部偏掉，而**结果看起来完全正常**
#' —— 没有报错、没有警告，只有错的数字。泛化到任意数据集后这是最容易踩的坑。
#'
#' 判据（按可靠性排序，任一命中即定论）：
#'   1. **出现负值** -> 一定是 log 尺度（线性强度不可能为负）
#'   2. **99 分位 > 100** -> 一定是线性（log2 的芯片强度极少超过 ~20）
#'   3. 否则看中位数：> 50 判线性，<= 50 判 log2
#'
#' 注意 3 是启发式。所以返回值里带 `reason`，调用方把它落盘 ——
#' 判断依据可核对，而不是"程序说是就是"。
#'
#' @param expr 表达矩阵（探针 x 样本）
#' @return list(scale = "log2"|"linear"|"unknown", reason, q01, q50, q99, neg_frac)
detect_expr_scale <- function(expr) {
  v <- as.numeric(as.matrix(expr))
  v <- v[is.finite(v)]
  if (length(v) == 0L) {
    return(list(scale = "unknown", reason = "矩阵里没有有限值",
                q01 = NA_real_, q50 = NA_real_, q99 = NA_real_, neg_frac = NA_real_))
  }
  q <- stats::quantile(v, c(0.01, 0.5, 0.99), na.rm = TRUE, names = FALSE)
  neg_frac <- mean(v < 0)
  base <- list(q01 = q[1L], q50 = q[2L], q99 = q[3L], neg_frac = neg_frac)

  if (neg_frac > 0.001) {
    return(c(list(scale = "log2",
                  reason = sprintf("有 %.2f%% 的值 < 0，线性荧光强度不可能为负", 100 * neg_frac)),
             base))
  }
  if (q[3L] > 100) {
    return(c(list(scale = "linear",
                  reason = sprintf("99 分位 = %.1f > 100（log2 芯片强度极少超过 ~20）", q[3L])),
             base))
  }
  if (q[2L] > 50) {
    return(c(list(scale = "linear",
                  reason = sprintf("中位数 = %.1f > 50", q[2L])),
             base))
  }
  c(list(scale = "log2",
         reason = sprintf("无负值，99 分位 = %.1f <= 100，中位数 = %.1f <= 50", q[3L], q[2L])),
    base)
}

#' 从临床/特征表里识别候选批次字段
#'
#' **为什么不能"没找到批次就直接往下跑"。** 用户文档那条警告是对的：
#' 批次与分组共线时做 ComBat 会把真实信号一起抹掉。但反过来的错误更常见 ——
#' **根本不查批次**，于是把"分组本身就是一个批次"当成生物学差异报出去。
#' 所以这里必须留下"查过、结论是什么"的记录，而不是静默跳过。
#'
#' 候选来自列名（GEO 的 characteristics 解析后列名形如 `batch`、`plate`、
#' `scan_date`）。**常量列被排除** —— 只有一个取值的列不是批次，是元数据噪声。
#'
#' @return 列名（character(1)）或 NULL
detect_batch_field <- function(clinical) {
  if (is.null(clinical) || ncol(clinical) == 0L) return(NULL)
  cand <- grep("batch|plate|chip|slide|scan|sentrix|hyb|process|run|date|center|site",
               colnames(clinical), ignore.case = TRUE, value = TRUE)
  if (length(cand) == 0L) return(NULL)
  cand <- cand[vapply(cand, function(cn) {
    v <- clinical[[cn]]
    length(unique(v[!is.na(v)])) > 1L
  }, logical(1))]
  if (length(cand) == 0L) return(NULL)
  cand[1L]
}

#' 评估批次与分组的混杂程度
#'
#' **判据是"每个组是否跨了多个批次"，不是相关系数。**
#' 完全共线（每组恰好落在一个批次里，且各组批次不重叠）时，
#' 批次效应与分组效应在数学上不可分离 —— 此时做批次校正等于
#' 把要检验的效应本身减掉，得到"校正后无差异"，而流程全绿。
#'
#' Cramér's V 作为辅助指标一并给出（V 越接近 1 关联越强），
#' 但**结论以"每组批次数"为准**：小样本下 V 的估计本身不稳。
#'
#' @return list(n_batch, n_group, batches_per_group, groups_per_batch,
#'              cramers_v, confounded, verdict, reason)
assess_batch_confounding <- function(batch, group) {
  keep <- !is.na(batch) & !is.na(group)
  batch <- as.character(batch[keep]); group <- as.character(group[keep])
  tb <- table(batch, group)
  n_batch <- nrow(tb); n_group <- ncol(tb)
  batches_per_group <- apply(tb, 2L, function(cc) sum(cc > 0L))
  groups_per_batch  <- apply(tb, 1L, function(rr) sum(rr > 0L))

  # Cramér's V：手算，不用 chisq.test —— 小样本/稀疏列联表下它会警告或报错，
  # 而那种报错会被当成"检查失败"，掩盖真正要看的结构。
  n <- sum(tb)
  v <- NA_real_
  if (n > 0L && n_batch > 1L && n_group > 1L) {
    e <- outer(rowSums(tb), colSums(tb)) / n
    ok <- e > 0
    chi2 <- sum((tb[ok] - e[ok])^2 / e[ok])
    v <- sqrt(chi2 / (n * min(n_batch - 1L, n_group - 1L)))
  }

  # 完全共线：没有任何一个组跨了 2 个以上批次
  confounded <- n_batch > 1L && n_group > 1L && all(batches_per_group <= 1L)
  list(
    n_batch = n_batch, n_group = n_group,
    batches_per_group = as.list(batches_per_group),
    groups_per_batch = as.list(groups_per_batch),
    cramers_v = if (is.na(v)) NA_real_ else round(v, 4),
    confounded = confounded,
    verdict = if (n_batch < 2L) "no_batch_variation"
              else if (confounded) "fully_confounded"
              else "partially_crossed",
    reason = if (n_batch < 2L)
      "只有一个批次（或没有批次信息），无从校正"
    else if (confounded)
      sprintf("完全共线：%d 个组各自只落在一个批次里 —— 批次与分组不可分离，校正会抹掉真实效应",
              n_group)
    else
      sprintf("可分离：至少一个组跨了 %d 个批次", max(batches_per_group))
  )
}

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
