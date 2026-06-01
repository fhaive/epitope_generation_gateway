#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "/data/fsluma/pipelines/Epitope_Generation_Gateway/TCGA_melanoma/sample_specific_networks"

summary_file <- file.path(
  BASE,
  "reviewer_A_combined_supervisor_summary",
  "combined_supervisor_summary.tsv"
)

if (!file.exists(summary_file)) {
  stop("Missing summary file: ", summary_file)
}

dt <- fread(summary_file)

cat("\n============================================================\n")
cat("Reviewer A / Network stability: main jackknife overview\n")
cat("============================================================\n\n")

cat("Design:\n")
cat("  30 target patients x 29 leave-one-background-sample-out perturbations = 870 comparisons\n")
cat("  Each jackknife metric ranking is compared against the original full-cohort ranking\n")
cat("  for the same target patient.\n\n")

cat("Main result table:\n\n")
print(dt)

cat("\nCompact interpretation:\n\n")

for (i in seq_len(nrow(dt))) {
  metric <- dt$Metric[i]
  med <- dt$`Median Spearman rho`[i]
  q025 <- dt$`2.5% Spearman rho`[i]
  minv <- dt$`Minimum Spearman rho`[i]
  top1 <- dt$`Top 1% median overlap`[i]
  top1_q025 <- dt$`Top 1% 2.5% overlap`[i]
  top1_min <- dt$`Top 1% minimum overlap`[i]

  cat(sprintf(
    "  %s: median Spearman rho = %.3f; 2.5%% = %.3f; minimum = %.3f; top-1%% median overlap = %.1f%%; top-1%% 2.5%% overlap = %.1f%%; top-1%% minimum overlap = %.1f%%.\n",
    metric, med, q025, minv, top1, top1_q025, top1_min
  ))
}

cat("\nReviewer-facing takeaway:\n")
cat("  WCI, degree, strength, and LCI show very high jackknife stability.\n")
cat("  Betweenness is more sensitive, as expected for a shortest-path-based metric,\n")
cat("  but remains concordant with the original rankings.\n\n")

cat("Most important numbers to cite:\n")
cat("  WCI median Spearman rho = 0.997; minimum = 0.990\n")
cat("  Degree median Spearman rho = 0.998; minimum = 0.994\n")
cat("  Strength median Spearman rho = 0.997; minimum = 0.989\n")
cat("  LCI median Spearman rho = 0.996; minimum = 0.987\n")
cat("  Betweenness median Spearman rho = 0.933; minimum = 0.862\n")
cat("  Median top-1% overlap across metrics = 91.8% to 95.0%\n\n")

cat("Files:\n")
cat("  ", summary_file, "\n")
cat("  ", file.path(BASE, "reviewer_A_section1_wci_merged/global_rank_summary.tsv"), "\n")
cat("  ", file.path(BASE, "reviewer_A_section1_wci_merged/global_topk_summary.tsv"), "\n")
cat("  ", file.path(BASE, "reviewer_A_2C_metrics_merged/global_filtered_metric_rank_summary.tsv"), "\n")
cat("  ", file.path(BASE, "reviewer_A_2C_metrics_merged/global_filtered_metric_topk_summary.tsv"), "\n")
cat("\n")
