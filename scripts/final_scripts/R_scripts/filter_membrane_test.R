# Filter_Membrane_Only.R
library(Matrix)
library(yaml)

args <- commandArgs(trailingOnly = TRUE)
input_rds     <- args[1]
output_rds    <- args[2]
output_mm     <- args[3]
output_genes  <- args[4]
config_file   <- args[5]
# humannet_file <- args[6]   # Not used for membrane filtering only
mem_genes_file <- args[7]

cat("Input RDS:", input_rds, "\n")
cat("Output RDS:", output_rds, "\n")
cat("Matrix Market:", output_mm, "\n")
cat("Gene list:", output_genes, "\n")
cat("Config File:", config_file, "\n")
cat("Membrane Genes File:", mem_genes_file, "\n")

# --- Config ---
config   <- read_yaml(config_file)
MEM_PROP <- config$mem_prop   # e.g. 0.05 for top 5%

# --- Load data ---
S <- readRDS(input_rds)
mem_genes <- as.character(readRDS(mem_genes_file))
G <- rownames(S)
n <- length(G)

mem_idx <- which(G %in% mem_genes)
if (length(mem_idx) == 0) stop("No membrane genes found in network!")

# All pairs (i, j) with i < j and either i or j is a membrane gene
all_pairs <- expand.grid(i = mem_idx, j = seq_len(n))
all_pairs <- all_pairs[all_pairs$i < all_pairs$j, , drop = FALSE]
mem_edges <- unique(as.matrix(all_pairs))

# Get edge weights
if (nrow(mem_edges) > 0) {
  mem_weights <- abs(S[mem_edges])
  cutoff <- quantile(mem_weights, 1 - MEM_PROP, na.rm = TRUE)
  keep_mem <- mem_edges[mem_weights >= cutoff, , drop = FALSE]
} else {
  keep_mem <- matrix(nrow = 0, ncol = 2)
}

# Build filtered matrix
if (nrow(keep_mem) == 0) {
  filtered_mat <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                              dims = dim(S), dimnames = dimnames(S))
} else {
  mem_i <- keep_mem[, 1]
  mem_j <- keep_mem[, 2]
  mem_x <- S[as.matrix(keep_mem)]
  filtered_mat <- sparseMatrix(
    i = c(mem_i, mem_j),
    j = c(mem_j, mem_i),
    x = c(mem_x, mem_x),
    dims = dim(S),
    dimnames = dimnames(S)
  )
}

# Save outputs
saveRDS(filtered_mat, output_rds)
Matrix::writeMM(filtered_mat, output_mm)
writeLines(rownames(filtered_mat), output_genes)
cat("✓ Done\n")
