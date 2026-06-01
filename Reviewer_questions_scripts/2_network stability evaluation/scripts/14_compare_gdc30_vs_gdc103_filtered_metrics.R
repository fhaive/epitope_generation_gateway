#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "/data/fsluma/pipelines/Epitope_Generation_Gateway/TCGA_melanoma/sample_specific_networks"
WORK <- file.path(ROOT, "reviewer_C_external_benchmark")

DIR30 <- file.path(WORK, "gdc_primary_original30_filtered_metrics")
DIR103 <- file.path(WORK, "tcga_primary_benchmark_filtered_metrics")

OUT <- file.path(WORK, "gdc30_vs_gdc103_filtered_metric_comparison")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TOP_PROPS <- c(0.01, 0.05, 0.10, 0.20)

METRICS <- list(
  degree = list(
    subdir = "Network_Metrics_Degree",
    suffix = "_degree.tsv",
    value_col = "degree"
  ),
  strength = list(
    subdir = "Network_Metrics_Strength",
    suffix = "_strength.tsv",
    value_col = "strength"
  ),
  betweenness = list(
    subdir = "Network_Metrics_Betweenness",
    suffix = "_betweenness.tsv",
    value_col = "betweenness"
  ),
  lci = list(
    subdir = "Network_Metrics_LargestComponentImpact",
    suffix = "_impact.tsv",
    value_col = "largest_component_impact"
  )
)

message2 <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ...)
  cat("\n")
  flush.console()
}

read_metric <- function(base, sample, metric_name) {
  cfg <- METRICS[[metric_name]]
  f <- file.path(base, cfg$subdir, paste0(sample, cfg$suffix))

  if (!file.exists(f)) {
    stop("Missing metric file: ", f)
  }

  dt <- fread(f)

  if (!("gene" %in% names(dt))) {
    stop("Missing gene column in: ", f)
  }

  if (!(cfg$value_col %in% names(dt))) {
    stop("Missing value column ", cfg$value_col, " in: ", f, ". Columns: ", paste(names(dt), collapse = ", "))
  }

  dt[, .(
    gene = as.character(gene),
    value = as.numeric(get(cfg$value_col))
  )]
}

top_genes <- function(dt, prop) {
  x <- copy(dt)
  x <- x[!is.na(value)]

  if (nrow(x) == 0) {
    return(character(0))
  }

  k <- max(1, ceiling(prop * nrow(x)))
  setorder(x, -value, gene)

  x$gene[seq_len(min(k, nrow(x)))]
}

get_samples <- function() {
  files <- list.files(
    file.path(DIR30, "Network_Metrics_Degree"),
    pattern = "_degree\\.tsv$",
    full.names = FALSE
  )

  sub("_degree\\.tsv$", "", files)
}

compare_one <- function(sample, metric_name) {
  x30 <- read_metric(DIR30, sample, metric_name)
  x103 <- read_metric(DIR103, sample, metric_name)

  m <- merge(
    x30,
    x103,
    by = "gene",
    suffixes = c("_gdc30", "_gdc103")
  )

  m <- m[!is.na(value_gdc30) & !is.na(value_gdc103)]

  if (nrow(m) < 3) {
    stop("Too few common genes for ", sample, " ", metric_name, ": ", nrow(m))
  }

  cor_dt <- data.table(
    comparison = "GDC-primary n30 vs GDC-primary n103",
    sample = sample,
    metric = metric_name,
    n_common_genes = nrow(m),
    spearman_rank_correlation = suppressWarnings(cor(m$value_gdc30, m$value_gdc103, method = "spearman")),
    kendall_rank_correlation = suppressWarnings(cor(m$value_gdc30, m$value_gdc103, method = "kendall"))
  )

  top_list <- list()

  for (p in TOP_PROPS) {
    top30 <- top_genes(x30, p)
    top103 <- top_genes(x103, p)

    inter <- length(intersect(top30, top103))
    union <- length(union(top30, top103))
    min_top <- min(length(top30), length(top103))

    top_list[[as.character(p)]] <- data.table(
      comparison = "GDC-primary n30 vs GDC-primary n103",
      sample = sample,
      metric = metric_name,
      top_prop = p,
      gdc30_top_n = length(top30),
      gdc103_top_n = length(top103),
      intersection_n = inter,
      union_n = union,
      overlap_fraction_of_gdc30_top = ifelse(length(top30) > 0, inter / length(top30), NA_real_),
      overlap_fraction_min_top = ifelse(min_top > 0, inter / min_top, NA_real_),
      jaccard = ifelse(union > 0, inter / union, NA_real_)
    )
  }

  list(
    cor = cor_dt,
    top = rbindlist(top_list, use.names = TRUE, fill = TRUE)
  )
}

samples <- sort(get_samples())

message2("Samples found in GDC n30 metric folder: ", length(samples))
print(samples)

cor_results <- list()
top_results <- list()
status <- list()

for (metric_name in names(METRICS)) {
  message2("Metric: ", metric_name)

  for (sample in samples) {
    task_start <- Sys.time()

    res <- tryCatch({
      comp <- compare_one(sample, metric_name)

      list(
        cor = comp$cor,
        top = comp$top,
        status = data.table(
          sample = sample,
          metric = metric_name,
          status = "success",
          elapsed_seconds = as.numeric(difftime(Sys.time(), task_start, units = "secs")),
          message = NA_character_
        )
      )
    }, error = function(e) {
      list(
        cor = data.table(),
        top = data.table(),
        status = data.table(
          sample = sample,
          metric = metric_name,
          status = "failed",
          elapsed_seconds = as.numeric(difftime(Sys.time(), task_start, units = "secs")),
          message = conditionMessage(e)
        )
      )
    })

    if (nrow(res$cor) > 0) {
      cor_results[[length(cor_results) + 1]] <- res$cor
    }

    if (nrow(res$top) > 0) {
      top_results[[length(top_results) + 1]] <- res$top
    }

    status[[length(status) + 1]] <- res$status
  }
}

cor_dt <- rbindlist(cor_results, use.names = TRUE, fill = TRUE)
top_dt <- rbindlist(top_results, use.names = TRUE, fill = TRUE)
status_dt <- rbindlist(status, use.names = TRUE, fill = TRUE)

fwrite(status_dt, file.path(OUT, "gdc30_vs_gdc103_filtered_metric_status.tsv"), sep = "\t")

cat("\nStatus summary:\n")
print(status_dt[, .N, by = .(metric, status)])

if (nrow(status_dt[status == "failed"]) > 0) {
  cat("\nFailed rows:\n")
  print(status_dt[status == "failed"])
}

if (nrow(cor_dt) == 0 || nrow(top_dt) == 0) {
  stop("No successful comparisons")
}

fwrite(cor_dt, file.path(OUT, "gdc30_vs_gdc103_filtered_metric_rank_correlations.tsv"), sep = "\t")
fwrite(top_dt, file.path(OUT, "gdc30_vs_gdc103_filtered_metric_topk_overlap.tsv"), sep = "\t")

rank_summary <- cor_dt[, .(
  n_samples = .N,
  median_spearman = median(spearman_rank_correlation, na.rm = TRUE),
  mean_spearman = mean(spearman_rank_correlation, na.rm = TRUE),
  min_spearman = min(spearman_rank_correlation, na.rm = TRUE),
  q025_spearman = quantile(spearman_rank_correlation, 0.025, na.rm = TRUE),
  q25_spearman = quantile(spearman_rank_correlation, 0.25, na.rm = TRUE),
  q75_spearman = quantile(spearman_rank_correlation, 0.75, na.rm = TRUE),
  q975_spearman = quantile(spearman_rank_correlation, 0.975, na.rm = TRUE),
  max_spearman = max(spearman_rank_correlation, na.rm = TRUE),
  median_kendall = median(kendall_rank_correlation, na.rm = TRUE),
  min_common_genes = min(n_common_genes, na.rm = TRUE),
  max_common_genes = max(n_common_genes, na.rm = TRUE)
), by = .(comparison, metric)]

top_summary <- top_dt[, .(
  n_samples = .N,
  median_overlap_fraction_of_gdc30_top = median(overlap_fraction_of_gdc30_top, na.rm = TRUE),
  mean_overlap_fraction_of_gdc30_top = mean(overlap_fraction_of_gdc30_top, na.rm = TRUE),
  min_overlap_fraction_of_gdc30_top = min(overlap_fraction_of_gdc30_top, na.rm = TRUE),
  q025_overlap_fraction_of_gdc30_top = quantile(overlap_fraction_of_gdc30_top, 0.025, na.rm = TRUE),
  q25_overlap_fraction_of_gdc30_top = quantile(overlap_fraction_of_gdc30_top, 0.25, na.rm = TRUE),
  q75_overlap_fraction_of_gdc30_top = quantile(overlap_fraction_of_gdc30_top, 0.75, na.rm = TRUE),
  q975_overlap_fraction_of_gdc30_top = quantile(overlap_fraction_of_gdc30_top, 0.975, na.rm = TRUE),
  max_overlap_fraction_of_gdc30_top = max(overlap_fraction_of_gdc30_top, na.rm = TRUE),
  median_overlap_fraction_min_top = median(overlap_fraction_min_top, na.rm = TRUE),
  median_jaccard = median(jaccard, na.rm = TRUE),
  min_jaccard = min(jaccard, na.rm = TRUE)
), by = .(comparison, metric, top_prop)][order(metric, top_prop)]

fwrite(rank_summary, file.path(OUT, "gdc30_vs_gdc103_filtered_metric_rank_summary.tsv"), sep = "\t")
fwrite(top_summary, file.path(OUT, "gdc30_vs_gdc103_filtered_metric_topk_summary.tsv"), sep = "\t")

display <- merge(
  rank_summary,
  top_summary[top_prop == 0.01, .(
    comparison,
    metric,
    top1_median_overlap = median_overlap_fraction_of_gdc30_top,
    top1_q025_overlap = q025_overlap_fraction_of_gdc30_top,
    top1_min_overlap = min_overlap_fraction_of_gdc30_top,
    top1_median_jaccard = median_jaccard
  )],
  by = c("comparison", "metric"),
  all.x = TRUE
)

metric_order <- c("degree", "strength", "lci", "betweenness")
display[, metric_order := match(metric, metric_order)]

display <- display[order(metric_order), .(
  Comparison = comparison,
  Metric = toupper(metric),
  `N samples` = n_samples,
  `Median Spearman rho` = round(median_spearman, 3),
  `2.5% Spearman rho` = round(q025_spearman, 3),
  `Minimum Spearman rho` = round(min_spearman, 3),
  `Median Kendall tau` = round(median_kendall, 3),
  `Top 1% median overlap` = round(100 * top1_median_overlap, 1),
  `Top 1% 2.5% overlap` = round(100 * top1_q025_overlap, 1),
  `Top 1% minimum overlap` = round(100 * top1_min_overlap, 1),
  `Top 1% median Jaccard` = round(top1_median_jaccard, 3),
  `Common genes min-max` = paste0(min_common_genes, "-", max_common_genes)
)]

fwrite(display, file.path(OUT, "gdc30_vs_gdc103_filtered_metric_supervisor_summary.tsv"), sep = "\t")

cat("\nFiltered metric benchmark summary:\n\n")
print(display)

cat("\nTop-k summary:\n\n")
print(top_summary)

message2("Done. Output: ", OUT)
