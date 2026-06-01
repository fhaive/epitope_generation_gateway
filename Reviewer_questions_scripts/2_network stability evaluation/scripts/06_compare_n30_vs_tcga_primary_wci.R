#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "/data/fsluma/pipelines/Epitope_Generation_Gateway/TCGA_melanoma/sample_specific_networks"
WORK <- file.path(ROOT, "reviewer_C_external_benchmark")

ORIG_WCI_DIR <- file.path(ROOT, "Network_Metrics_Full", "WCI")
BENCH_WCI_DIR <- file.path(WORK, "tcga_primary_benchmark_lioness_wci", "Network_Metrics_Full", "WCI")

OUT <- file.path(WORK, "tcga_primary_benchmark_wci_comparison")
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

message2 <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ...)
  cat("\n")
  flush.console()
}

read_wci <- function(f) {
  if (!file.exists(f)) {
    stop("Missing WCI file: ", f)
  }

  dt <- fread(f)

  if (nrow(dt) == 0) {
    stop("Empty WCI file: ", f)
  }

  gene_candidates <- c("gene", "Gene", "GENE", "ensembl_gene_id", "gene_id")
  gene_col <- gene_candidates[gene_candidates %in% names(dt)][1]

  if (is.na(gene_col)) {
    gene_col <- names(dt)[1]
  }

  value_candidates <- c(
    "wci", "WCI", "value", "Value",
    "weighted_connectivity_index",
    "Weighted_Connectivity_Index",
    "weighted_connectivity"
  )

  value_col <- value_candidates[value_candidates %in% names(dt)][1]

  if (is.na(value_col)) {
    stop(
      "Could not find WCI/value column in: ", f,
      "\nColumns are: ", paste(names(dt), collapse = ", ")
    )
  }

  out <- dt[, .(
    gene = as.character(get(gene_col)),
    value = as.numeric(get(value_col))
  )]

  out <- out[!is.na(gene) & nzchar(gene)]
  out
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

compare_one_sample <- function(sample_id) {
  orig_file <- file.path(ORIG_WCI_DIR, paste0(sample_id, "_wci.tsv"))
  bench_file <- file.path(BENCH_WCI_DIR, paste0(sample_id, "_wci.tsv"))

  orig <- read_wci(orig_file)
  bench <- read_wci(bench_file)

  m <- merge(
    orig,
    bench,
    by = "gene",
    suffixes = c("_n30", "_tcga_primary")
  )

  m <- m[!is.na(value_n30) & !is.na(value_tcga_primary)]

  if (nrow(m) < 3) {
    stop("Too few common genes after merging WCI files for ", sample_id, ": ", nrow(m))
  }

  spearman <- suppressWarnings(cor(m$value_n30, m$value_tcga_primary, method = "spearman"))
  kendall <- suppressWarnings(cor(m$value_n30, m$value_tcga_primary, method = "kendall"))

  cor_dt <- data.table(
    sample = sample_id,
    metric = "wci",
    benchmark = "TCGA-SKCM Primary Tumor",
    n_common_genes = nrow(m),
    spearman_rank_correlation = spearman,
    kendall_rank_correlation = kendall
  )

  topk_list <- list()

  for (p in TOP_PROPS) {
    orig_top <- top_genes(orig, p)
    bench_top <- top_genes(bench, p)

    inter <- length(intersect(orig_top, bench_top))
    union <- length(union(orig_top, bench_top))
    min_top <- min(length(orig_top), length(bench_top))

    topk_list[[as.character(p)]] <- data.table(
      sample = sample_id,
      metric = "wci",
      benchmark = "TCGA-SKCM Primary Tumor",
      top_prop = p,
      n30_top_n = length(orig_top),
      benchmark_top_n = length(bench_top),
      intersection_n = inter,
      union_n = union,
      overlap_fraction_of_n30_top = ifelse(length(orig_top) > 0, inter / length(orig_top), NA_real_),
      overlap_fraction_min_top = ifelse(min_top > 0, inter / min_top, NA_real_),
      jaccard = ifelse(union > 0, inter / union, NA_real_)
    )
  }

  list(
    cor = cor_dt,
    topk = rbindlist(topk_list, use.names = TRUE, fill = TRUE)
  )
}

message2("Original WCI dir: ", ORIG_WCI_DIR)
message2("Benchmark WCI dir: ", BENCH_WCI_DIR)

message2("Original WCI files found: ", length(list.files(ORIG_WCI_DIR, pattern = "_wci\\.tsv$")))
message2("Benchmark WCI files found: ", length(list.files(BENCH_WCI_DIR, pattern = "_wci\\.tsv$")))

cor_results <- list()
topk_results <- list()
status <- list()

for (i in seq_along(TARGETS)) {
  sample_id <- TARGETS[i]
  message2("Comparing WCI target ", i, "/", length(TARGETS), ": ", sample_id)

  task_start <- Sys.time()

  res <- tryCatch({
    comp <- compare_one_sample(sample_id)

    list(
      cor = comp$cor,
      topk = comp$topk,
      status = data.table(
        sample = sample_id,
        status = "success",
        elapsed_seconds = as.numeric(difftime(Sys.time(), task_start, units = "secs")),
        message = NA_character_
      )
    )
  }, error = function(e) {
    list(
      cor = data.table(),
      topk = data.table(),
      status = data.table(
        sample = sample_id,
        status = "failed",
        elapsed_seconds = as.numeric(difftime(Sys.time(), task_start, units = "secs")),
        message = conditionMessage(e)
      )
    )
  })

  if (nrow(res$cor) > 0) {
    cor_results[[length(cor_results) + 1]] <- res$cor
  }

  if (nrow(res$topk) > 0) {
    topk_results[[length(topk_results) + 1]] <- res$topk
  }

  status[[length(status) + 1]] <- res$status
}

cor_dt <- rbindlist(cor_results, use.names = TRUE, fill = TRUE)
topk_dt <- rbindlist(topk_results, use.names = TRUE, fill = TRUE)
status_dt <- rbindlist(status, use.names = TRUE, fill = TRUE)

fwrite(status_dt, file.path(OUT, "wci_comparison_status.tsv"), sep = "\t")

cat("\nStatus summary:\n")
print(status_dt[, .N, by = status])

if (nrow(status_dt[status == "failed"]) > 0) {
  cat("\nFirst failed rows:\n")
  print(head(status_dt[status == "failed"], 10))
}

if (nrow(cor_dt) == 0 || nrow(topk_dt) == 0) {
  stop(
    "No successful WCI comparisons. See: ",
    file.path(OUT, "wci_comparison_status.tsv")
  )
}

fwrite(cor_dt, file.path(OUT, "wci_rank_correlations_n30_vs_tcga_primary.tsv"), sep = "\t")
fwrite(topk_dt, file.path(OUT, "wci_topk_overlap_n30_vs_tcga_primary.tsv"), sep = "\t")

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
), by = .(benchmark, metric)]

topk_summary <- topk_dt[, .(
  n_samples = .N,
  median_overlap_fraction_of_n30_top = median(overlap_fraction_of_n30_top, na.rm = TRUE),
  mean_overlap_fraction_of_n30_top = mean(overlap_fraction_of_n30_top, na.rm = TRUE),
  min_overlap_fraction_of_n30_top = min(overlap_fraction_of_n30_top, na.rm = TRUE),
  q025_overlap_fraction_of_n30_top = quantile(overlap_fraction_of_n30_top, 0.025, na.rm = TRUE),
  q25_overlap_fraction_of_n30_top = quantile(overlap_fraction_of_n30_top, 0.25, na.rm = TRUE),
  q75_overlap_fraction_of_n30_top = quantile(overlap_fraction_of_n30_top, 0.75, na.rm = TRUE),
  q975_overlap_fraction_of_n30_top = quantile(overlap_fraction_of_n30_top, 0.975, na.rm = TRUE),
  max_overlap_fraction_of_n30_top = max(overlap_fraction_of_n30_top, na.rm = TRUE),
  median_jaccard = median(jaccard, na.rm = TRUE),
  min_jaccard = min(jaccard, na.rm = TRUE)
), by = .(benchmark, metric, top_prop)][order(top_prop)]

fwrite(rank_summary, file.path(OUT, "wci_rank_summary_n30_vs_tcga_primary.tsv"), sep = "\t")
fwrite(topk_summary, file.path(OUT, "wci_topk_summary_n30_vs_tcga_primary.tsv"), sep = "\t")

display <- merge(
  rank_summary,
  topk_summary[top_prop == 0.01, .(
    benchmark,
    metric,
    top1_median_overlap = median_overlap_fraction_of_n30_top,
    top1_q025_overlap = q025_overlap_fraction_of_n30_top,
    top1_min_overlap = min_overlap_fraction_of_n30_top,
    top1_median_jaccard = median_jaccard
  )],
  by = c("benchmark", "metric"),
  all.x = TRUE
)

display_clean <- display[, .(
  Benchmark = benchmark,
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

fwrite(display_clean, file.path(OUT, "wci_benchmark_supervisor_summary.tsv"), sep = "\t")

cat("\nWCI n=30 vs TCGA-primary benchmark summary:\n\n")
print(display_clean)

cat("\nTop-k summary:\n\n")
print(topk_summary)

message2("Done. Output: ", OUT)
