# ============================================================================
# 00_validate_inputs.R — 前置校验与分组
# ============================================================================
# spec 的 validate_inputs 步骤 + fetch_geo 的 validate 部分。
#
# 这是整个流水线的硬门禁：数据集必须是人源、芯片、乳腺癌，且满足
# config 里 design_mode 对应的样本量判据：
#   small_sample  总样本 < 10，每组 >= 3
#   cohort        总样本 >= 15，每组 >= 10
# 任一不满足直接 stop()，不进入任何分析。
#
# 只用 base R 抓 GEO SOFT 元数据，因此这一步在装任何 Bioconductor 包之前就能跑完。
# 输出：data/group.csv, data/meta.csv, data/platform.txt,
#       data/clinical.csv, data/clinical_fields.json, data/geo_metadata.json
# ============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
})

# ---- bootstrap: 定位并加载 lib/common.R（Rscript 与 source 两种方式都适用）--
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

# 两套门禁，由 config 的 design_mode 选择。**不是同一个门禁换个阈值** ——
# 两种设计的失败模式不同，所以判据不同：
#
#   small_sample  探索性小样本。上限 10 是硬性的：n >= 10 时全基因组 BH 校正
#                 通常能出结果，就不再需要 ranked_fallback 那套降级路径和
#                 相应的措辞约束，应该走 cohort 模式。
#   cohort        队列级。下限 15 来自 WGCNA 的通行要求（低于 15 相关矩阵
#                 不稳定，模块划分不可复现）。每组下限 10 是 limma 给出
#                 可用统计量的实际要求。
MIN_PER_GROUP        <- 3L    # small_sample：每组最少样本数
SMALL_SAMPLE_MAX     <- 10L   # small_sample：总样本数上限（不含）
COHORT_MIN_SAMPLES   <- 15L   # cohort：总样本数下限（WGCNA 通行下限）
COHORT_MIN_PER_GROUP <- 10L   # cohort：每组最少样本数

#' 把 SOFT 行按 ^SAMPLE 切成每样本一块
split_soft_samples <- function(lines) {
  starts <- grep("^\\^SAMPLE = ", lines)
  if (length(starts) == 0L) return(list())
  ends <- c(starts[-1L] - 1L, length(lines))
  lapply(seq_along(starts), function(i) lines[starts[i]:ends[i]])
}

#' 判定一个样本属于哪个组；命中 0 个或多个组都返回相应诊断
classify_sample <- function(values, group_values) {
  haystack <- tolower(paste(values, collapse = " | "))
  hit <- names(group_values)[vapply(group_values, function(pats) {
    any(vapply(tolower(pats), function(p) grepl(p, haystack, fixed = TRUE), logical(1)))
  }, logical(1))]
  hit
}

run_00_validate_inputs <- function(cfg) {
  log_info("=== 步骤 00：输入校验与分组 ===")
  ensure_dirs(cfg)

  gse <- cfg$dataset_id
  log_info(sprintf("数据集: %s  分组字段: %s  对比: %s vs %s",
                   gse, cfg$group_field, cfg$contrast[1], cfg$contrast[2]))

  # ---- 1. series 级元数据 --------------------------------------------------
  series_lines <- fetch_geo_soft(gse, targ = "self")
  series_type  <- soft_value(series_lines, "Series_type")
  platform_id  <- soft_value(series_lines, "Series_platform_id")
  title        <- soft_value(series_lines, "Series_title")
  sample_taxid <- soft_values(series_lines, "Series_sample_taxid")
  design       <- soft_value(series_lines, "Series_overall_design", "")

  if (is.na(series_type)) stop("GEO 未返回 Series_type，无法校验数据类型")
  if (!is.na(cfg$platform_id) && nzchar(cfg$platform_id) &&
      !is.na(platform_id) && !identical(platform_id, cfg$platform_id)) {
    stop(sprintf("配置的 platform_id=%s 与 GEO 实际平台 %s 不一致",
                 cfg$platform_id, platform_id))
  }

  # ---- 2. 硬门禁：数据类型必须为芯片 --------------------------------------
  if (grepl("sequencing", series_type, ignore.case = TRUE)) {
    stop(sprintf(
      "数据类型不合规：%s 的类型是「%s」，属于测序而非基因芯片。\n  spec 要求 GES/GEO 基因芯片，本流水线不适用于 RNA-seq。",
      gse, series_type))
  }
  if (!grepl("array", series_type, ignore.case = TRUE)) {
    stop(sprintf("数据类型不合规：%s 的类型是「%s」，不是 expression profiling by array",
                 gse, series_type))
  }
  log_info(sprintf("数据类型 OK: %s", series_type))

  # ---- 3. 样本级元数据与分组 ----------------------------------------------
  gsm_lines <- fetch_geo_soft(gse, targ = "gsm")
  blocks <- split_soft_samples(gsm_lines)
  if (length(blocks) == 0L) stop(sprintf("未能从 %s 解析出任何样本", gse))

  meta <- do.call(rbind, lapply(blocks, function(b) {
    data.frame(
      gsm          = soft_value(b, "Sample_geo_accession"),
      title        = soft_value(b, "Sample_title"),
      source_name  = soft_value(b, "Sample_source_name_ch1", ""),
      organism     = soft_value(b, "Sample_organism_ch1", ""),
      taxid        = soft_value(b, "Sample_taxid_ch1", ""),
      group_field  = paste(soft_values(b, paste0("Sample_", cfg$group_field)), collapse = " | "),
      stringsAsFactors = FALSE
    )
  }))

  # ---- 4. 硬门禁：物种必须为人源 ------------------------------------------
  organisms <- unique(meta$organism[nzchar(meta$organism)])
  if (length(organisms) == 0L) {
    # 回退到 series 级 taxid
    if (!any(grepl("9606", sample_taxid))) {
      stop(sprintf("物种不合规：无法确认 %s 为人源（taxid 9606）", gse))
    }
    organisms <- "Homo sapiens"
  }
  if (!all(grepl("Homo sapiens", organisms, fixed = TRUE))) {
    stop(sprintf("物种不合规：%s 包含非人源样本: %s", gse, paste(organisms, collapse = ", ")))
  }
  log_info(sprintf("物种 OK: %s (%d 个样本)", paste(organisms, collapse = ", "), nrow(meta)))

  # ---- 5. 硬门禁：样本量（按 design_mode 分别判定）-------------------------
  n <- nrow(meta)
  if (identical(cfg$design_mode, "small_sample")) {
    if (n >= SMALL_SAMPLE_MAX) {
      stop(sprintf(
        "样本量不合规：design_mode=small_sample 要求 < %d，%s 共 %d 个样本。\n  这个规模应该用 design_mode: cohort（下限 %d）。\n  不要为了让它通过而调大上限 —— 两种设计的降级路径和措辞约束不同。",
        SMALL_SAMPLE_MAX, gse, n, COHORT_MIN_SAMPLES))
    }
    log_info(sprintf("样本量 OK: %d (< %d, design_mode=small_sample)", n, SMALL_SAMPLE_MAX))
  } else {
    if (n < COHORT_MIN_SAMPLES) {
      stop(sprintf(
        "样本量不合规：design_mode=cohort 要求 >= %d，%s 只有 %d 个样本。\n  低于此规模 WGCNA 的相关矩阵不稳定、模块划分不可复现。\n  这个规模应该用 design_mode: small_sample。",
        COHORT_MIN_SAMPLES, gse, n))
    }
    log_info(sprintf("样本量 OK: %d (>= %d, design_mode=cohort)", n, COHORT_MIN_SAMPLES))
  }

  # ---- 6. 分组 ------------------------------------------------------------
  hits <- lapply(strsplit(meta$group_field, " \\| "), classify_sample, cfg$group_values)
  meta$group <- vapply(hits, function(h) if (length(h) == 1L) h else NA_character_, character(1))

  ambiguous <- which(vapply(hits, length, integer(1)) > 1L)
  if (length(ambiguous) > 0L) {
    stop(sprintf(
      "分组歧义：以下样本同时命中多个组，请收紧 config 的 group_values 模式：\n  %s",
      paste(sprintf("%s (%s) -> %s", meta$gsm[ambiguous], meta$title[ambiguous],
                    vapply(hits[ambiguous], paste, character(1), collapse = "+")),
            collapse = "\n  ")))
  }
  unmatched <- which(is.na(meta$group))
  if (length(unmatched) > 0L) {
    stop(sprintf(
      "分组失败：以下样本未命中任何组，请检查 config 的 group_field / group_values：\n  %s\n  样本 %s 的实际取值: %s",
      paste(sprintf("%s (%s)", meta$gsm[unmatched], meta$title[unmatched]), collapse = "\n  "),
      cfg$group_field,
      paste(unique(meta$group_field[unmatched]), collapse = " ;; ")))
  }

  counts <- table(factor(meta$group, levels = names(cfg$group_values)))
  log_info(sprintf("分组结果: %s", paste(sprintf("%s=%d", names(counts), as.integer(counts)),
                                          collapse = ", ")))

  # ---- 7. 硬门禁：每组样本数下限（按 design_mode）--------------------------
  min_per_group <- if (identical(cfg$design_mode, "cohort")) COHORT_MIN_PER_GROUP else MIN_PER_GROUP
  thin <- names(counts)[as.integer(counts) < min_per_group]
  if (length(thin) > 0L) {
    stop(sprintf("分组样本数不足：组 %s 的样本数 < %d（design_mode=%s），limma 无法给出有意义的统计量",
                 paste(thin, collapse = ", "), min_per_group, cfg$design_mode))
  }
  if (!all(cfg$contrast %in% names(cfg$group_values))) {
    stop("contrast 引用了未定义的组")
  }

  # ---- 8. 配对关系（可选）-------------------------------------------------
  #
  # 配对关系在 config 里**显式声明**，不从样本标题里猜后缀 ——
  # "TNBC tissue 1" / "matched normal breast tissues 1" 这种后缀看着能对上，
  # 但换个数据集就失效，而且猜错了不会报错，只会让配对分析静默变错。
  patient <- rep(NA_character_, n)
  if (isTRUE(cfg$paired)) {
    pairs <- cfg$pairs
    if (is.null(pairs) || length(pairs) == 0L) {
      stop("config 里 paired: true 但没给 pairs；配对关系必须显式声明")
    }
    seen <- character(0)
    for (k in seq_along(pairs)) {
      pr <- unlist(pairs[[k]])
      if (length(pr) != 2L) stop(sprintf("pairs[[%d]] 必须恰好两个 GSM", k))
      if (!all(pr %in% meta$gsm)) {
        stop(sprintf("pairs[[%d]] 含不属于本数据集的样本: %s", k,
                     paste(setdiff(pr, meta$gsm), collapse = ", ")))
      }
      if (any(pr %in% seen)) stop(sprintf("样本在 pairs 里重复出现: %s",
                                          paste(intersect(pr, seen), collapse = ", ")))
      seen <- c(seen, pr)
      grp <- meta$group[match(pr, meta$gsm)]
      if (!setequal(grp, cfg$contrast)) {
        stop(sprintf("pairs[[%d]] (%s) 必须一例来自 %s、一例来自 %s，实际为 %s",
                     k, paste(pr, collapse = "+"), cfg$contrast[1L], cfg$contrast[2L],
                     paste(grp, collapse = "+")))
      }
      patient[match(pr, meta$gsm)] <- sprintf("pair%02d", k)
    }
    if (anyNA(patient)) {
      stop(sprintf("paired: true 但以下样本没有配对: %s",
                   paste(meta$gsm[is.na(patient)], collapse = ", ")))
    }
    log_info(sprintf("配对校验通过: %d 对（每对一例 %s、一例 %s）",
                     length(pairs), cfg$contrast[1L], cfg$contrast[2L]))
  }

  # ---- 8b. 临床字段（全量，不只是分组字段）--------------------------------
  #
  # **这一步原来一个临床字段都没落盘。** `group_field` 只保留了 config 指定的
  # 那一个，其余 characteristics（生存时间、生存事件、年龄、分级、ER 状态）
  # 在 step 00 之后就消失了 —— 于是 step 07 做 LASSO/Cox 时没有任何终点可用。
  #
  # 到 step 07 再回头去下 SOFT 会把"数据源"拆成两处（同一份数据两个入口），
  # 所以在这里一次抽全，之后所有步骤都从 clinical.csv 读。
  # 解析器在 common.R（`clinical_table()`），因为 step 07 的外部验证队列
  # 也要用同一份 —— 两处各写一份迟早会在"同名 key 怎么办"上分叉。
  clinical <- clinical_table(blocks, meta$gsm)
  utils::write.csv(clinical, file.path(cfg$output$data_dir, "clinical.csv"), row.names = FALSE)

  field_profile <- clinical_field_profile(clinical)
  write_json(file.path(cfg$output$data_dir, "clinical_fields.json"), field_profile)
  log_info(sprintf("临床字段: %d 个（%s）", length(field_profile),
                   paste(utils::head(names(field_profile), 6), collapse = ", ")))

  # ---- 9. 落盘 ------------------------------------------------------------
  meta_out <- meta[, c("gsm", "title", "source_name", "organism", "taxid", "group", "group_field")]
  group_out <- data.frame(gsm = meta_out$gsm, group = meta_out$group, stringsAsFactors = FALSE)
  if (isTRUE(cfg$paired)) group_out$patient <- patient
  utils::write.csv(meta_out, file.path(cfg$output$data_dir, "meta.csv"), row.names = FALSE)
  utils::write.csv(group_out, file.path(cfg$output$data_dir, "group.csv"), row.names = FALSE)
  writeLines(c(sprintf("platform_id\t%s", platform_id),
               sprintf("platform_title\t%s", soft_value(series_lines, "Series_platform_title", "NA")),
               sprintf("series_type\t%s", series_type),
               sprintf("n_samples\t%d", n)),
             file.path(cfg$output$data_dir, "platform.txt"))

  write_json(file.path(cfg$output$data_dir, "geo_metadata.json"), list(
    dataset_id = gse,
    title = title,
    series_type = series_type,
    platform_id = platform_id,
    n_samples = n,
    design_mode = cfg$design_mode,
    min_per_group = min_per_group,
    organisms = organisms,
    group_counts = as.list(as.integer(counts)),
    group_names = names(counts),
    overall_design = design,
    contrast = cfg$contrast,
    paired = cfg$paired,
    validated = TRUE,
    validated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  ))

  log_info(sprintf("已写出 %s/group.csv, meta.csv, platform.txt, clinical.csv, clinical_fields.json, geo_metadata.json",
                   cfg$output$data_dir))
  invisible(meta_out)
}

if (!GEO_ORCHESTRATED()) {
  cfg <- load_config()
  run_00_validate_inputs(cfg)
}
