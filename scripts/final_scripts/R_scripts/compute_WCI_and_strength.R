#!/usr/bin/env Rscript

library(Matrix)
library(dplyr)
library(tibble)
library(readr)

# --- Command-line arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript compute_WCI_and_strength.R input_rds output_strength_tsv output_wci_tsv")
}
input_rds <- args[1]
output_strength <- args[2]
output_wci <- args[3]

# --- Load network ---
cat("Input file:", input_rds, "\n")
S <- readRDS(input_rds)
genes <- rownames(S)

# --- Compute metrics ---
deg   <- rowSums(S != 0)
strength <- rowSums(abs(S))
total_weight <- sum(abs(S))

# --- Write strength file ---
df_strength <- tibble(
  gene = genes,
  degree = deg,
  strength = strength
)
write_tsv(df_strength, output_strength)
cat("Wrote strength file:", output_strength, "\n")

# --- Weighted Connectivity Impact (WCI) ---
if (total_weight > 0) {
  wci <- strength / total_weight
} else {
  wci <- rep(NA_real_, length(strength))
}
df_wci <- tibble(
  gene = genes,
  wci = wci
)
write_tsv(df_wci, output_wci)
cat("Wrote WCI file:", output_wci, "\n")
cat("Done.\n")
