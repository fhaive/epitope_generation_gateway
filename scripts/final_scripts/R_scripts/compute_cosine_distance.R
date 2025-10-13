#!/usr/bin/env Rscript

suppressMessages({
  library(Matrix)
  library(igraph)
  library(proxyC)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript compute_cosine_distance.R output_prefix label input_files...")
}

out_prefix <- args[1]
label <- args[2]
input_files <- args[3:length(args)]

# First, collect all unique gene names across all samples
all_genes <- character()
sample_strengths <- list()

cat(sprintf("Found %d sample files\n", length(input_files)))

for (file in input_files) {
  net <- readRDS(file)
  if (inherits(net, "igraph")) {
    adj <- as_adjacency_matrix(net, attr="weight", sparse=TRUE)
    genes <- rownames(adj)
    str_i <- rowSums(abs(adj))
  } else if (inherits(net, "Matrix") || inherits(net, "matrix")) {
    adj <- as(net, "dgCMatrix")
    genes <- rownames(adj)
    str_i <- rowSums(abs(adj))
  } else {
    stop("Input files must be igraph or matrix objects.")
  }
  all_genes <- union(all_genes, genes)
  sample_strengths[[file]] <- list(str=str_i, genes=genes)
  rm(net, adj, str_i, genes)
  gc()
}

# Ensure all genes present as rows, all samples as columns
all_genes <- sort(unique(all_genes))
n_genes <- length(all_genes)
n_samples <- length(input_files)
strength_mat <- matrix(0, nrow=n_genes, ncol=n_samples,
                       dimnames=list(all_genes, basename(input_files)))

# Fill by gene names
for (i in seq_along(input_files)) {
  dat <- sample_strengths[[input_files[i]]]
  str_vec <- dat$str
  genes_vec <- dat$genes
  # Only fill existing genes for this sample
  strength_mat[genes_vec, i] <- str_vec
}

cat("Computing cosine similarity and distance...\n")
sim_mat <- proxyC::simil(strength_mat, method = "cosine", margin = 2L)
dist_mat <- 1 - as.matrix(sim_mat)

cat("Saving distance matrix and plots...\n")
# Save the matrix as RDS and CSV
saveRDS(dist_mat, file = paste0(out_prefix, "_dist_matrix.rds"))
write.csv(dist_mat, file = paste0(out_prefix, "_dist_matrix.csv"), row.names = TRUE)

# Save heatmaps
pheatmap(dist_mat,
         main = sprintf("Cosine distance (%s)", label),
         cluster_rows = TRUE, cluster_cols = TRUE,
         fontsize = 8,
         filename = paste0(out_prefix, "_dist_heatmap.png"))

pheatmap(sim_mat,
         main = sprintf("Cosine similarity (%s)", label),
         cluster_rows = TRUE, cluster_cols = TRUE,
         fontsize = 8,
         filename = paste0(out_prefix, "_sim_heatmap.png"))

cat("Done!\n")
