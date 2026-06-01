#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "/data/fsluma/pipelines/Epitope_Generation_Gateway/TCGA_melanoma/sample_specific_networks"

WCI_RANK <- file.path(BASE, "reviewer_A_section1_wci_merged", "global_rank_summary.tsv")
WCI_TOPK <- file.path(BASE, "reviewer_A_section1_wci_merged", "global_topk_summary.tsv")

FILT_RANK <- file.path(BASE, "reviewer_A_2C_metrics_merged", "global_filtered_metric_rank_summary.tsv")
FILT_TOPK <- file.path(BASE, "reviewer_A_2C_metrics_merged", "global_filtered_metric_topk_summary.tsv")

OUT <- file.path(BASE, "reviewer_A_combined_supervisor_summary")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

wci_rank <- fread(WCI_RANK)
wci_topk <- fread(WCI_TOPK)

filt_rank <- fread(FILT_RANK)
filt_topk <- fread(FILT_TOPK)

wci_rank[, network_level := "Full LIONESS network"]
filt_rank[, network_level := "PPI/membrane-filtered network"]

rank_all <- rbindlist(list(wci_rank, filt_rank), use.names = TRUE, fill = TRUE)

wci_top1 <- wci_topk[top_prop == 0.01]
filt_top1 <- filt_topk[top_prop == 0.01]

top1_all <- rbindlist(list(wci_top1, filt_top1), use.names = TRUE, fill = TRUE)

combined <- merge(
  rank_all,
  top1_all[, .(
    metric,
    top1_median_overlap = median_overlap_fraction,
    top1_mean_overlap = mean_overlap_fraction,
    top1_min_overlap = min_overlap_fraction,
    top1_q025_overlap = q025_overlap_fraction,
    top1_q25_overlap = q25_overlap_fraction,
    top1_q75_overlap = q75_overlap_fraction,
    top1_q975_overlap = q975_overlap_fraction,
    top1_median_jaccard = median_jaccard,
    top1_min_jaccard = min_jaccard
  )],
  by = "metric",
  all.x = TRUE
)

# Nice ordering
metric_order <- c("wci", "degree", "strength", "lci", "betweenness")
combined[, metric := factor(metric, levels = metric_order)]
setorder(combined, metric)
combined[, metric := as.character(metric)]

# Clean display table
display <- combined[, .(
  Metric = toupper(metric),
  `Network level` = network_level,
  `N jackknife comparisons` = n_replicates,
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

fwrite(combined, file.path(OUT, "combined_full_numeric_summary.tsv"), sep = "\t")
fwrite(display, file.path(OUT, "combined_supervisor_summary.tsv"), sep = "\t")

# Markdown table for easy copy/paste
md_file <- file.path(OUT, "combined_supervisor_summary.md")

con <- file(md_file, open = "wt")

writeLines("| Metric | Network level | N | Median Spearman ρ | 2.5% Spearman ρ | Minimum Spearman ρ | Median Kendall τ | Top 1% median overlap | Top 1% 2.5% overlap | Top 1% minimum overlap | Top 1% median Jaccard | Common genes |", con)
writeLines("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|", con)

for (i in seq_len(nrow(display))) {
  row <- display[i]
  writeLines(paste0(
    "| ", row$Metric,
    " | ", row$`Network level`,
    " | ", row$`N jackknife comparisons`,
    " | ", row$`Median Spearman rho`,
    " | ", row$`2.5% Spearman rho`,
    " | ", row$`Minimum Spearman rho`,
    " | ", row$`Median Kendall tau`,
    " | ", row$`Top 1% median overlap`, "%",
    " | ", row$`Top 1% 2.5% overlap`, "%",
    " | ", row$`Top 1% minimum overlap`, "%",
    " | ", row$`Top 1% median Jaccard`,
    " | ", row$`Common genes min-max`,
    " |"
  ), con)
}

close(con)

cat("\nSupervisor-friendly table:\n\n")
print(display)

cat("\nMarkdown table:\n\n")
cat(readLines(md_file), sep = "\n")

cat("\n\nFiles written:\n")
cat(file.path(OUT, "combined_supervisor_summary.tsv"), "\n")
cat(file.path(OUT, "combined_supervisor_summary.md"), "\n")
cat(file.path(OUT, "combined_full_numeric_summary.tsv"), "\n")
