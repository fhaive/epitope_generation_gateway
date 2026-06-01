#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

ROOT <- "/data/fsluma/pipelines/Epitope_Generation_Gateway"

OUT <- file.path(
  ROOT,
  "TCGA_melanoma/rescue_final_analysis/reviewer_C5_membrane_impact/weighted_borda_ppi_only_ablation"
)

validation <- fread(file.path(OUT, "weighted_borda_validation_rank_correlations.tsv"))
ablation <- fread(file.path(OUT, "weighted_borda_ppi_only_ablation_rank_correlations.tsv"))
topk <- fread(file.path(OUT, "weighted_borda_topk_overlap.tsv"))
promoted <- fread(file.path(OUT, "membrane_promoted_examples_per_patient_summary.tsv"))

validation <- validation[, .(
  sample,
  n_epitopes,
  validation_spearman_final_vs_recomputed = spearman,
  validation_kendall_final_vs_recomputed = kendall
)]

ablation <- ablation[, .(
  sample,
  ablation_spearman_final_vs_ppi_only = spearman,
  ablation_kendall_final_vs_ppi_only = kendall
)]

top10 <- topk[
  comparison == "ablation_final_membrane_vs_ppi_only" & top_n == 10,
  .(
    sample,
    top10_overlap_n = top_overlap_n,
    top10_overlap_fraction = top_overlap_fraction,
    top10_jaccard = top_jaccard,
    original_top10_membrane_n = a_top_membrane_n,
    ppi_only_top10_membrane_n = b_top_membrane_n
  )
]

top20 <- topk[
  comparison == "ablation_final_membrane_vs_ppi_only" & top_n == 20,
  .(
    sample,
    top20_overlap_n = top_overlap_n,
    top20_overlap_fraction = top_overlap_fraction,
    top20_jaccard = top_jaccard,
    original_top20_membrane_n = a_top_membrane_n,
    ppi_only_top20_membrane_n = b_top_membrane_n
  )
]

promoted <- promoted[, .(
  sample,
  total_epitopes,
  membrane_epitopes,
  membrane_promoted_n,
  membrane_promoted_gain_ge_5_n,
  membrane_promoted_gain_ge_10_n,
  membrane_promoted_into_top10_n = membrane_in_original_top10_not_ppi_only_n,
  membrane_promoted_into_top20_n = membrane_in_original_top20_not_ppi_only_n
)]

tab <- Reduce(
  function(x, y) merge(x, y, by = "sample", all = TRUE),
  list(validation, ablation, top10, top20, promoted)
)

setorder(tab, sample)

display <- copy(tab)

num_cols <- names(display)[vapply(display, is.numeric, logical(1))]
for (cc in num_cols) {
  if (grepl("fraction|spearman|kendall|jaccard", cc)) {
    display[, (cc) := round(get(cc), 3)]
  }
}

# Add percentage columns for easier reading.
display[, top10_overlap_percent := round(100 * top10_overlap_fraction, 1)]
display[, top20_overlap_percent := round(100 * top20_overlap_fraction, 1)]

# Reorder columns for supplement readability.
ordered_cols <- c(
  "sample",
  "n_epitopes",
  "total_epitopes",
  "membrane_epitopes",

  "validation_spearman_final_vs_recomputed",
  "validation_kendall_final_vs_recomputed",

  "ablation_spearman_final_vs_ppi_only",
  "ablation_kendall_final_vs_ppi_only",

  "top10_overlap_n",
  "top10_overlap_percent",
  "top10_jaccard",
  "original_top10_membrane_n",
  "ppi_only_top10_membrane_n",

  "top20_overlap_n",
  "top20_overlap_percent",
  "top20_jaccard",
  "original_top20_membrane_n",
  "ppi_only_top20_membrane_n",

  "membrane_promoted_n",
  "membrane_promoted_gain_ge_5_n",
  "membrane_promoted_gain_ge_10_n",
  "membrane_promoted_into_top10_n",
  "membrane_promoted_into_top20_n"
)

ordered_cols <- ordered_cols[ordered_cols %in% names(display)]
display <- display[, ..ordered_cols]

out_file <- file.path(OUT, "Table_S_C5_membrane_expansion_ablation_summary.tsv")
fwrite(display, out_file, sep = "\t")

# Also create a compact median/summary table for quick checking only.
summary <- data.table(
  n_patients = nrow(display),
  median_validation_spearman = median(tab$validation_spearman_final_vs_recomputed, na.rm = TRUE),
  min_validation_spearman = min(tab$validation_spearman_final_vs_recomputed, na.rm = TRUE),
  median_ablation_spearman = median(tab$ablation_spearman_final_vs_ppi_only, na.rm = TRUE),
  min_ablation_spearman = min(tab$ablation_spearman_final_vs_ppi_only, na.rm = TRUE),
  median_ablation_kendall = median(tab$ablation_kendall_final_vs_ppi_only, na.rm = TRUE),
  median_top10_overlap_n = median(tab$top10_overlap_n, na.rm = TRUE),
  median_top10_overlap_percent = median(100 * tab$top10_overlap_fraction, na.rm = TRUE),
  median_top20_overlap_n = median(tab$top20_overlap_n, na.rm = TRUE),
  median_top20_overlap_percent = median(100 * tab$top20_overlap_fraction, na.rm = TRUE),
  total_membrane_epitopes = sum(tab$membrane_epitopes, na.rm = TRUE),
  total_membrane_promoted_n = sum(tab$membrane_promoted_n, na.rm = TRUE),
  total_membrane_promoted_into_top10_n = sum(tab$membrane_promoted_into_top10_n, na.rm = TRUE),
  total_membrane_promoted_into_top20_n = sum(tab$membrane_promoted_into_top20_n, na.rm = TRUE)
)

fwrite(
  summary,
  file.path(OUT, "Table_S_C5_membrane_expansion_ablation_summary_MEDIANS_FOR_RESPONSE.tsv"),
  sep = "\t"
)

cat("\nWrote main supplementary table:\n")
cat(out_file, "\n\n")

cat("Summary for response:\n")
print(summary)

cat("\nFirst rows of supplementary table:\n")
print(display[1:min(.N, 10)])
