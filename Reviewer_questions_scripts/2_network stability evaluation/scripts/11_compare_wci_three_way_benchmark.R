#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "/data/fsluma/pipelines/Epitope_Generation_Gateway/TCGA_melanoma/sample_specific_networks"
WORK <- file.path(ROOT, "reviewer_C_external_benchmark")

DIR_ORIGINAL_N30 <- file.path(ROOT, "Network_Metrics_Full", "WCI")
DIR_GDC_N30 <- file.path(WORK, "gdc_primary_original30_lioness_wci", "Network_Metrics_Full", "WCI")
DIR_GDC_N103 <- file.path(WORK, "tcga_primary_benchmark_lioness_wci", "Network_Metrics_Full", "WCI")

OUT <- file.path(WORK, "wci_three_way_benchmark_comparison")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TOP_PROPS <- c(0.01, 0.05, 0.10, 0.20)

TARGETS <- c(
  "Sample_TCGA-BF-A1PV",
  "Sample_TCGA-BF-A3DJ",
  "Sample_TCGA-BF-A3DL",
  "Sample_TCGA-BF-A3DM",
  "Sample_TCGA-BF-A3DN",
  "Sample_TCGA-BF-A5EO",
  "Sample_TCGA-BF-A5EP",
  "Sample_TCGA-BF-A5EQ",
  "Sample_TCGA-BF-A5ER",
  "Sample_TCGA-BF-A5ES",
  "Sample_TCGA-BF-AAOU",
  "Sample_TCGA-BF-AAOX",
  "Sample_TCGA-BF-AAP1",
  "Sample_TCGA-BF-AAP2",
  "Sample_TCGA-BF-AAP4",
  "Sample_TCGA-BF-AAP6",
  "Sample_TCGA-BF-AAP7",
  "Sample_TCGA-BF-AAP8",
  "Sample_TCGA-D3-A5GT",
  "Sample_TCGA-D9-A3Z4",
  "Sample_TCGA-D9-A4Z2",
  "Sample_TCGA-D9-A4Z3",
  "Sample_TCGA-EB-A1NK",
  "Sample_TCGA-EB-A3HV",
  "Sample_TCGA-EB-A3XB",
  "Sample_TCGA-EB-A3XC",
  "Sample_TCGA-EB-A3XD",
  "Sample_TCGA-EB-A3XE",
  "Sample_TCGA-EB-A3XF",
  "Sample_TCGA-EB-A41A"
)

COMPARISONS <- list(
  list(
    label = "Original EGG n30 vs GDC-primary n30",
    left_name = "Original EGG n30",
    right_name = "GDC-primary n30",
    left_dir = DIR_ORIGINAL_N30,
    right_dir = DIR_GDC_N30
  ),
  list(
    label = "GDC-primary n30 vs GDC-primary n103",
    left_name = "GDC-primary n30",
    right_name = "GDC-primary n103",
    left_dir = DIR_GDC_N30,
    right_dir = DIR_GDC_N103
  ),
  list(
    label = "Original EGG n30 vs GDC-primary n103",
    left_name = "Original EGG n30",
    right_name = "GDC-primary n103",
    left_dir = DIR_ORIGINAL_N30,
    right_dir = DIR_GDC_N103
  )
)

message2 <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ...)
  cat("\n")
  flush.console()
}

read_wci <- function(f) {
  if (!file.exists(f)) stop("Missing WCI file: ", f)

  dt <- fread(f)

  if (!all(c("gene", "wci") %in% names(dt))) {
    stop("Expected columns gene,wci in: ", f, " columns=", paste(names(dt), collapse = ","))
  }

  dt[, .(gene = as.character(gene), value = as.numeric(wci))]
}

top_genes <- function(dt, prop) {
  x <- copy(dt)
  x <- x[!is.na(value)]
  k <- max(1, ceiling(prop * nrow(x)))
  setorder(x, -value, gene)
  x$gene[seq_len(min(k, nrow(x)))]
}

compare_one <- function(sample_id, cmp) {
  left_file <- file.path(cmp$left_dir, paste0(sample_id, "_wci.tsv"))
  right_file <- file.path(cmp$right_dir, paste0(sample_id, "_wci.tsv"))

  left <- read_wci(left_file)
  right <- read_wci(right_file)

  m <- merge(left, right, by = "gene", suffixes = c("_left", "_right"))
  m <- m[!is.na(value_left) & !is.na(value_right)]

  if (nrow(m) < 3) {
    stop("Too few common genes for ", sample_id, " in ", cmp$label, ": ", nrow(m))
  }

  cor_dt <- data.table(
    comparison = cmp$label,
    left = cmp$left_name,
    right = cmp$right_name,
    sample = sample_id,
    metric = "wci",
    n_common_genes = nrow(m),
    spearman_rank_correlation = suppressWarnings(cor(m$value_left, m$value_right, method = "spearman")),
    kendall_rank_correlation = suppressWarnings(cor(m$value_left, m$value_right, method = "kendall"))
  )

  top_list <- list()

  for (p in TOP_PROPS) {
    left_top <- top_genes(left, p)
    right_top <- top_genes(right, p)

    inter <- length(intersect(left_top, right_top))
    union <- length(union(left_top, right_top))
    min_top <- min(length(left_top), length(right_top))

    top_list[[as.character(p)]] <- data.table(
      comparison = cmp$label,
      left = cmp$left_name,
      right = cmp$right_name,
      sample = sample_id,
      metric = "wci",
      top_prop = p,
      left_top_n = length(left_top),
      right_top_n = length(right_top),
      intersection_n = inter,
      union_n = union,
      overlap_fraction_of_left_top = ifelse(length(left_top) > 0, inter / length(left_top), NA_real_),
      overlap_fraction_min_top = ifelse(min_top > 0, inter / min_top, NA_real_),
      jaccard = ifelse(union > 0, inter / union, NA_real_)
    )
  }

  list(
    cor = cor_dt,
    top = rbindlist(top_list, use.names = TRUE, fill = TRUE)
  )
}

cor_results <- list()
top_results <- list()
status <- list()

for (cmp in COMPARISONS) {
  message2("Comparison: ", cmp$label)

  for (sample_id in TARGETS) {
    task_start <- Sys.time()

    res <- tryCatch({
      out <- compare_one(sample_id, cmp)

      list(
        cor = out$cor,
        top = out$top,
        status = data.table(
          comparison = cmp$label,
          sample = sample_id,
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
          comparison = cmp$label,
          sample = sample_id,
          status = "failed",
          elapsed_seconds = as.numeric(difftime(Sys.time(), task_start, units = "secs")),
          message = conditionMessage(e)
        )
      )
    })

    if (nrow(res$cor) > 0) cor_results[[length(cor_results) + 1]] <- res$cor
    if (nrow(res$top) > 0) top_results[[length(top_results) + 1]] <- res$top
    status[[length(status) + 1]] <- res$status
  }
}

cor_dt <- rbindlist(cor_results, use.names = TRUE, fill = TRUE)
top_dt <- rbindlist(top_results, use.names = TRUE, fill = TRUE)
status_dt <- rbindlist(status, use.names = TRUE, fill = TRUE)

fwrite(cor_dt, file.path(OUT, "wci_three_way_rank_correlations.tsv"), sep = "\t")
fwrite(top_dt, file.path(OUT, "wci_three_way_topk_overlap.tsv"), sep = "\t")
fwrite(status_dt, file.path(OUT, "wci_three_way_status.tsv"), sep = "\t")

cat("\nStatus summary:\n")
print(status_dt[, .N, by = .(comparison, status)])

if (nrow(cor_dt) == 0) {
  stop("No successful comparisons")
}

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
), by = .(comparison, left, right, metric)]

top_summary <- top_dt[, .(
  n_samples = .N,
  median_overlap_fraction_of_left_top = median(overlap_fraction_of_left_top, na.rm = TRUE),
  q025_overlap_fraction_of_left_top = quantile(overlap_fraction_of_left_top, 0.025, na.rm = TRUE),
  min_overlap_fraction_of_left_top = min(overlap_fraction_of_left_top, na.rm = TRUE),
  median_overlap_fraction_min_top = median(overlap_fraction_min_top, na.rm = TRUE),
  median_jaccard = median(jaccard, na.rm = TRUE),
  min_jaccard = min(jaccard, na.rm = TRUE)
), by = .(comparison, left, right, metric, top_prop)][order(comparison, top_prop)]

fwrite(rank_summary, file.path(OUT, "wci_three_way_rank_summary.tsv"), sep = "\t")
fwrite(top_summary, file.path(OUT, "wci_three_way_topk_summary.tsv"), sep = "\t")

display <- merge(
  rank_summary,
  top_summary[top_prop == 0.01, .(
    comparison,
    left,
    right,
    metric,
    top1_median_overlap = median_overlap_fraction_of_left_top,
    top1_q025_overlap = q025_overlap_fraction_of_left_top,
    top1_min_overlap = min_overlap_fraction_of_left_top,
    top1_median_jaccard = median_jaccard
  )],
  by = c("comparison", "left", "right", "metric"),
  all.x = TRUE
)

display <- display[, .(
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

fwrite(display, file.path(OUT, "wci_three_way_supervisor_summary.tsv"), sep = "\t")

cat("\nThree-way WCI benchmark summary:\n\n")
print(display)

cat("\nTop-k summary:\n\n")
print(top_summary)

message2("Done. Output: ", OUT)
