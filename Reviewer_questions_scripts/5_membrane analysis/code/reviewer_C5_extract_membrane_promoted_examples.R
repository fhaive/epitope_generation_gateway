#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "/data/fsluma/pipelines/Epitope_Generation_Gateway"

IN_DIR <- file.path(
  ROOT,
  "TCGA_melanoma/rescue_final_analysis/reviewer_C5_membrane_impact/weighted_borda_ppi_only_ablation/ablation_final_epitopes"
)

OUT <- file.path(
  ROOT,
  "TCGA_melanoma/rescue_final_analysis/reviewer_C5_membrane_impact/weighted_borda_ppi_only_ablation"
)

files <- list.files(
  IN_DIR,
  pattern = "_weighted_ppi_only_ablation_epitopes\\.csv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No ablation final epitope files found in: ", IN_DIR)
}

sample_from_file <- function(f) {
  x <- basename(f)
  x <- sub("_weighted_ppi_only_ablation_epitopes\\.csv$", "", x)
  x
}

all <- rbindlist(lapply(files, function(f) {
  dt <- fread(f)
  dt[, sample := sample_from_file(f)]
  dt
}), use.names = TRUE, fill = TRUE)

needed <- c(
  "sample",
  "Gene.Name",
  "HLA.Allele",
  "Peptide.Length",
  "MT.Epitope.Seq",
  "WT.Epitope.Seq",
  "Borda_Rank",
  "ppi_only_weighted_borda_rank",
  "rank_gain_with_membrane",
  "is_membrane_gene"
)

missing <- setdiff(needed, names(all))
if (length(missing) > 0) {
  stop("Missing columns: ", paste(missing, collapse = ", "))
}

all[, Borda_Rank := as.numeric(Borda_Rank)]
all[, ppi_only_weighted_borda_rank := as.numeric(ppi_only_weighted_borda_rank)]
all[, rank_gain_with_membrane := as.numeric(rank_gain_with_membrane)]

# Positive rank_gain_with_membrane means:
# ppi_only_rank - final_rank > 0, so the membrane-expanded ranking is better.
membrane_promoted <- all[
  is_membrane_gene == TRUE &
    !is.na(rank_gain_with_membrane) &
    rank_gain_with_membrane > 0
]

# Strongest individual membrane-derived upward shifts, regardless of top-k boundary.
membrane_promoted_examples <- membrane_promoted[
  order(-rank_gain_with_membrane, Borda_Rank)
]

# Cases where a membrane epitope enters top-k only after membrane expansion.
topk_list <- list()

for (k in c(10, 20, 50, 100)) {
  x <- all[
    is_membrane_gene == TRUE &
      Borda_Rank <= k &
      ppi_only_weighted_borda_rank > k
  ]

  if (nrow(x) > 0) {
    x[, top_n := k]
    topk_list[[as.character(k)]] <- x
  }
}

membrane_promoted_into_topk <- if (length(topk_list) > 0) {
  rbindlist(topk_list, use.names = TRUE, fill = TRUE)
} else {
  data.table()
}

# Cases where membrane expansion improves a membrane epitope by at least selected thresholds.
threshold_summary <- rbindlist(lapply(c(1, 2, 5, 10, 20, 50), function(th) {
  x <- all[
    is_membrane_gene == TRUE &
      !is.na(rank_gain_with_membrane) &
      rank_gain_with_membrane >= th
  ]

  data.table(
    min_rank_gain = th,
    n_epitopes = nrow(x),
    n_patients = uniqueN(x$sample),
    median_final_rank = ifelse(nrow(x) > 0, median(x$Borda_Rank, na.rm = TRUE), NA_real_),
    best_final_rank = ifelse(nrow(x) > 0, min(x$Borda_Rank, na.rm = TRUE), NA_real_),
    median_rank_gain = ifelse(nrow(x) > 0, median(x$rank_gain_with_membrane, na.rm = TRUE), NA_real_),
    max_rank_gain = ifelse(nrow(x) > 0, max(x$rank_gain_with_membrane, na.rm = TRUE), NA_real_)
  )
}), use.names = TRUE)

# Per-patient counts.
per_patient_summary <- all[, .(
  total_epitopes = .N,
  membrane_epitopes = sum(is_membrane_gene == TRUE, na.rm = TRUE),
  membrane_promoted_n = sum(is_membrane_gene == TRUE & rank_gain_with_membrane > 0, na.rm = TRUE),
  membrane_promoted_gain_ge_5_n = sum(is_membrane_gene == TRUE & rank_gain_with_membrane >= 5, na.rm = TRUE),
  membrane_promoted_gain_ge_10_n = sum(is_membrane_gene == TRUE & rank_gain_with_membrane >= 10, na.rm = TRUE),
  membrane_in_original_top10_not_ppi_only_n = sum(is_membrane_gene == TRUE & Borda_Rank <= 10 & ppi_only_weighted_borda_rank > 10, na.rm = TRUE),
  membrane_in_original_top20_not_ppi_only_n = sum(is_membrane_gene == TRUE & Borda_Rank <= 20 & ppi_only_weighted_borda_rank > 20, na.rm = TRUE),
  membrane_in_original_top50_not_ppi_only_n = sum(is_membrane_gene == TRUE & Borda_Rank <= 50 & ppi_only_weighted_borda_rank > 50, na.rm = TRUE)
), by = sample][order(sample)]

overall_summary <- data.table(
  n_patients = uniqueN(all$sample),
  total_epitopes = nrow(all),
  total_membrane_epitopes = sum(all$is_membrane_gene == TRUE, na.rm = TRUE),
  membrane_promoted_n = nrow(membrane_promoted),
  patients_with_any_membrane_promoted = uniqueN(membrane_promoted$sample),
  membrane_promoted_into_top10_n = sum(per_patient_summary$membrane_in_original_top10_not_ppi_only_n),
  patients_with_membrane_promoted_into_top10 = uniqueN(per_patient_summary[membrane_in_original_top10_not_ppi_only_n > 0, sample]),
  membrane_promoted_into_top20_n = sum(per_patient_summary$membrane_in_original_top20_not_ppi_only_n),
  patients_with_membrane_promoted_into_top20 = uniqueN(per_patient_summary[membrane_in_original_top20_not_ppi_only_n > 0, sample]),
  membrane_promoted_into_top50_n = sum(per_patient_summary$membrane_in_original_top50_not_ppi_only_n),
  patients_with_membrane_promoted_into_top50 = uniqueN(per_patient_summary[membrane_in_original_top50_not_ppi_only_n > 0, sample])
)

# Compact example columns for reviewer/supplement.
example_cols <- c(
  "sample",
  "Gene.Name",
  "HLA.Allele",
  "Peptide.Length",
  "MT.Epitope.Seq",
  "WT.Epitope.Seq",
  "Median.MT.IC50.Score",
  "Tumor DNA VAF",
  "Borda_Rank",
  "ppi_only_weighted_borda_rank",
  "rank_gain_with_membrane",
  "Borda_Score",
  "ppi_only_weighted_borda_score",
  "Net_Degree",
  "ppi_only_degree",
  "Net_Strength",
  "ppi_only_strength",
  "Net_Impact",
  "ppi_only_impact",
  "Net_Betweenness",
  "ppi_only_betweenness",
  "Net_WCI",
  "is_membrane_gene"
)

example_cols <- example_cols[example_cols %in% names(all)]

fwrite(
  overall_summary,
  file.path(OUT, "membrane_promoted_examples_overall_summary.tsv"),
  sep = "\t"
)

fwrite(
  threshold_summary,
  file.path(OUT, "membrane_promoted_rank_gain_threshold_summary.tsv"),
  sep = "\t"
)

fwrite(
  per_patient_summary,
  file.path(OUT, "membrane_promoted_examples_per_patient_summary.tsv"),
  sep = "\t"
)

fwrite(
  membrane_promoted_examples[, ..example_cols],
  file.path(OUT, "membrane_promoted_epitopes_all.tsv"),
  sep = "\t"
)

fwrite(
  membrane_promoted_examples[1:min(.N, 50), ..example_cols],
  file.path(OUT, "membrane_promoted_epitopes_top50_examples.tsv"),
  sep = "\t"
)

if (nrow(membrane_promoted_into_topk) > 0) {
  topk_cols <- c("top_n", example_cols)
  topk_cols <- topk_cols[topk_cols %in% names(membrane_promoted_into_topk)]

  fwrite(
    membrane_promoted_into_topk[order(top_n, Borda_Rank), ..topk_cols],
    file.path(OUT, "membrane_promoted_into_topk_epitopes.tsv"),
    sep = "\t"
  )
} else {
  fwrite(
    data.table(),
    file.path(OUT, "membrane_promoted_into_topk_epitopes.tsv"),
    sep = "\t"
  )
}

cat("\nOverall summary:\n")
print(overall_summary)

cat("\nRank-gain threshold summary:\n")
print(threshold_summary)

cat("\nPer-patient summary, patients with promoted membrane epitopes in top 20:\n")
print(per_patient_summary[membrane_in_original_top20_not_ppi_only_n > 0])

cat("\nTop promoted membrane examples:\n")
print(membrane_promoted_examples[1:min(.N, 20), ..example_cols])
