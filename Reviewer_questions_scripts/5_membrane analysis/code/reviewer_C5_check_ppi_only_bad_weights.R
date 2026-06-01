#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})

ROOT <- "/data/fsluma/pipelines/Epitope_Generation_Gateway"
IN_DIR <- file.path(
  ROOT,
  "TCGA_melanoma/rescue_final_analysis/reviewer_C5_membrane_impact/ppi_only_filtered_networks/filtered_networks_rds"
)

files <- list.files(IN_DIR, pattern = "_ppi_only_filtered\\.rds$", full.names = TRUE)

res <- rbindlist(lapply(files, function(f) {
  g <- readRDS(f)
  w <- E(g)$weight

  if (is.null(w)) {
    w <- rep(1, ecount(g))
  }

  data.table(
    file = basename(f),
    vertices = vcount(g),
    edges = ecount(g),
    weight_na = sum(is.na(w)),
    weight_nan = sum(is.nan(w)),
    weight_inf = sum(is.infinite(w)),
    weight_nonfinite = sum(!is.finite(w)),
    weight_zero = sum(is.finite(w) & abs(w) == 0),
    weight_finite = sum(is.finite(w))
  )
}), use.names = TRUE, fill = TRUE)

print(res)
cat("\nSummary:\n")
print(res[, .(
  graphs = .N,
  total_edges = sum(edges),
  total_nonfinite = sum(weight_nonfinite),
  max_nonfinite_per_graph = max(weight_nonfinite),
  total_zero = sum(weight_zero)
)])

fwrite(
  res,
  file.path(
    ROOT,
    "TCGA_melanoma/rescue_final_analysis/reviewer_C5_membrane_impact/ppi_only_filtered_networks/ppi_only_bad_weight_qc.tsv"
  ),
  sep = "\t"
)
