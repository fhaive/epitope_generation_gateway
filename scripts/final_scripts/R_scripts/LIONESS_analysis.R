library(Matrix)
library(progress)
library(tidyverse)

# Get command-line arguments: input vst_mat.csv and output directory
args <- commandArgs(trailingOnly = TRUE)
input_file_rel <- args[1]  # Relative path to vst_mat.csv
output_dir_rel <- args[2]  # Relative path to output directory (Sample_Specific_Networks)

# Derive project root by getting script's directory and navigating up
full_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", full_args[grep("^--file=", full_args)])
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "../../.."))

# Debug: Print project root and relative paths
cat("Project root:", project_root, "\n")
cat("Relative input file:", input_file_rel, "\n")
cat("Relative output directory:", output_dir_rel, "\n")

# Resolve relative paths to absolute paths
input_file <- normalizePath(file.path(project_root, input_file_rel))
output_dir <- normalizePath(file.path(project_root, output_dir_rel), mustWork = FALSE)

# Debug: Print resolved absolute paths
cat("Absolute input file:", input_file, "\n")
cat("Absolute output directory:", output_dir, "\n")

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read input expression matrix
expr_mat <- read.csv(input_file, row.names = 1, check.names = FALSE)

# Network function for Pearson correlation
netFun <- function(m) cor(t(m), method = "pearson")
make_sparse <- function(M) as(Matrix(M, sparse = TRUE), "dgCMatrix")

# LIONESS computation
N <- ncol(expr_mat)
G_alpha <- netFun(expr_mat)

# Progress bar
pb <- progress_bar$new(total = N,
                       format = "  [:bar] :current/:total ETA: :eta")

for (q in seq_len(N)) {
  samp <- colnames(expr_mat)[q]
  pb$tick()

  # Leave-one-out global minus sample q
  G_minus <- netFun(expr_mat[, -q, drop = FALSE])

  # LIONESS formula
  E_q <- N * (G_alpha - G_minus) + G_minus

  # Convert to sparse and save
  S_q <- make_sparse(E_q)

  # Save as RDS
  saveRDS(S_q, file = file.path(output_dir, paste0(samp, ".rds")))

  # Clean up before next iteration
  rm(S_q, E_q, G_minus)
  gc()
}
