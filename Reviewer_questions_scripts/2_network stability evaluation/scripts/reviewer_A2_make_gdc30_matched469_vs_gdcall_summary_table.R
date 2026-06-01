#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

ROOT <- "/data/fsluma/pipelines/Epitope_Generation_Gateway"

SRC <- file.path(
  ROOT,
  "TCGA_melanoma/sample_specific_networks/reviewer_C_external_benchmark/gdc30_matched469_vs_gdcall_all_metrics_comparison"
)

OUT <- file.path(
  ROOT,
  "TCGA_melanoma/rescue_final_analysis/reviewer_A2_network_stability_external_gdc469_benchmark"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

supervisor_file <- file.path(SRC, "gdc30_matched469_vs_gdcall_all_metrics_supervisor_summary.tsv")
topk_file <- file.path(SRC, "gdc30_matched469_vs_gdcall_all_metrics_topk_summary.tsv")

if (!file.exists(supervisor_file)) {
  stop("Missing supervisor summary: ", supervisor_file)
}

if (!file.exists(topk_file)) {
  stop("Missing top-k summary: ", topk_file)
}

sup <- fread(supervisor_file)
topk <- fread(topk_file)

# Standardize metric names.
sup[, metric_key := tolower(Metric)]
topk[, metric_key := tolower(metric)]

metric_order <- c("wci", "degree", "strength", "lci", "betweenness")

metric_label <- c(
  wci = "WCI",
  degree = "Degree",
  strength = "Strength",
  lci = "LCI",
  betweenness = "Betweenness"
)

network_level <- c(
  wci = "Full LIONESS network",
  degree = "PPI/membrane-filtered network",
  strength = "PPI/membrane-filtered network",
  lci = "PPI/membrane-filtered network",
  betweenness = "PPI/membrane-filtered network"
)

# Extract top 20% summary from top-k table.
top20 <- topk[abs(as.numeric(top_prop) - 0.20) < 1e-9]

top20 <- top20[, .(
  metric_key,
  top20_median_overlap = as.numeric(median_overlap_fraction),
  top20_q025_overlap = as.numeric(q025_overlap_fraction),
  top20_min_overlap = as.numeric(min_overlap_fraction),
  top20_median_jaccard = as.numeric(median_jaccard),
  top20_median_intersection_n = as.numeric(median_intersection_n),
  top20_median_top_n = as.numeric(median_top_n)
)]

tab <- merge(sup, top20, by = "metric_key", all.x = TRUE)

tab[, metric_order := match(metric_key, metric_order)]
setorder(tab, metric_order)

# Numeric table keeps proportions numeric for downstream use.
numeric <- tab[, .(
  metric = metric_label[metric_key],
  benchmark_comparison = Comparison,
  network_level = network_level[metric_key],
  n_samples = as.integer(`N samples`),
  median_spearman = as.numeric(`Median Spearman rho`),
  q025_spearman = as.numeric(`2.5% Spearman rho`),
  min_spearman = as.numeric(`Minimum Spearman rho`),
  median_kendall = as.numeric(`Median Kendall tau`),
  top1_median_overlap = as.numeric(`Top 1% median overlap`) / 100,
  top1_q025_overlap = as.numeric(`Top 1% 2.5% overlap`) / 100,
  top1_min_overlap = as.numeric(`Top 1% minimum overlap`) / 100,
  top1_median_jaccard = as.numeric(`Top 1% median Jaccard`),
  top20_median_overlap,
  top20_q025_overlap,
  top20_min_overlap,
  top20_median_jaccard,
  n30_genes_min_max = `n30 genes min-max`,
  n469_genes_min_max = `n469 genes min-max`,
  common_genes_min_max = `Common genes min-max`
)]

pct <- function(x) {
  ifelse(is.na(x), NA_character_, paste0(round(100 * x, 1), "%"))
}

display <- numeric[, .(
  Metric = metric,
  `Benchmark comparison` = benchmark_comparison,
  `Network level` = network_level,
  N = n_samples,
  `Median Spearman rho` = round(median_spearman, 3),
  `2.5% Spearman rho` = round(q025_spearman, 3),
  `Minimum Spearman rho` = round(min_spearman, 3),
  `Median Kendall tau` = round(median_kendall, 3),
  `Top 1% median overlap` = pct(top1_median_overlap),
  `Top 1% 2.5% overlap` = pct(top1_q025_overlap),
  `Top 1% minimum overlap` = pct(top1_min_overlap),
  `Top 1% median Jaccard` = round(top1_median_jaccard, 3),
  `Top 20% median overlap` = pct(top20_median_overlap),
  `Top 20% 2.5% overlap` = pct(top20_q025_overlap),
  `Top 20% minimum overlap` = pct(top20_min_overlap),
  `Top 20% median Jaccard` = round(top20_median_jaccard, 3),
  `n30 genes min-max` = n30_genes_min_max,
  `n469 genes min-max` = n469_genes_min_max,
  `Common genes min-max` = common_genes_min_max
)]

fwrite(
  numeric,
  file.path(OUT, "Table_S_A2_external_benchmark_gdc30_matched469_vs_gdcall_n469_summary_numeric.tsv"),
  sep = "\t"
)

fwrite(
  display,
  file.path(OUT, "Table_S_A2_external_benchmark_gdc30_matched469_vs_gdcall_n469_summary.tsv"),
  sep = "\t"
)

# Copy source summaries for provenance.
fwrite(
  sup,
  file.path(OUT, "source_gdc30_matched469_vs_gdcall_supervisor_summary.tsv"),
  sep = "\t"
)

fwrite(
  topk,
  file.path(OUT, "source_gdc30_matched469_vs_gdcall_topk_summary.tsv"),
  sep = "\t"
)

cat("\nWrote:\n")
cat(file.path(OUT, "Table_S_A2_external_benchmark_gdc30_matched469_vs_gdcall_n469_summary.tsv"), "\n\n")

cat("Display table:\n")
print(display)
