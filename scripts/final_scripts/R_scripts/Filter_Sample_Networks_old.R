library(yaml)
library(stringr)
library(Matrix)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
input_rds <- args[1]
output_rds <- args[2]
output_mm <- args[3]
output_genes <- args[4]
config_file <- args[5]
humannet_file <- args[6]
mem_genes_file <- args[7]

cat("Input RDS:", input_rds, "\n")
cat("Output RDS:", output_rds, "\n")
cat("Matrix Market:", output_mm, "\n")
cat("Gene list:", output_genes, "\n")
cat("Config File:", config_file, "\n")
cat("HumanNet File:", humannet_file, "\n")
cat("Membrane Genes File:", mem_genes_file, "\n")

# Load config
config <- read_yaml(config_file)
top_prop <- config$top_prop
mem_prop <- config$mem_prop

# Load HumanNet (PPI) edge list
humannet <- readRDS(humannet_file)
ppi_edges <- humannet %>%
  transmute(g1 = as.character(ensembl1),
            g2 = as.character(ensembl2)) %>%
  unique()

# Load membrane gene list
mem_genes <- as.character(readRDS(mem_genes_file))

# Read dense matrix
samp <- str_remove(basename(input_rds), "\\.rds$")
cat("Processing sample:", samp, "\n")
S_full <- readRDS(input_rds)
gene_names <- rownames(S_full)
n <- nrow(S_full)
filtered_mat <- matrix(0, nrow = n, ncol = n, dimnames = dimnames(S_full))

# Make a lookup: gene_name -> integer index
gene_idx <- setNames(seq_along(gene_names), gene_names)

# Preprocess: For each node, get its PPI neighbors
cat("Building PPI neighbor lists...\n")
ppi_neighbors <- vector("list", n)
names(ppi_neighbors) <- gene_names

for (edge in seq_len(nrow(ppi_edges))) {
  g1 <- ppi_edges$g1[edge]
  g2 <- ppi_edges$g2[edge]
  if (!is.na(gene_idx[g1]) && !is.na(gene_idx[g2])) {
    # Add each as neighbor to the other
    ppi_neighbors[[g1]] <- c(ppi_neighbors[[g1]], g2)
    ppi_neighbors[[g2]] <- c(ppi_neighbors[[g2]], g1)
  }
}

cat("Filtering top % PPI edges per node...\n")
for (i in seq_len(n)) {
  node1 <- gene_names[i]
  nbrs <- unique(ppi_neighbors[[node1]])
  nbr_idx <- gene_idx[nbrs]
  nbr_idx <- nbr_idx[!is.na(nbr_idx)]
  if (length(nbr_idx) > 0) {
    weights <- abs(S_full[i, nbr_idx])
    cutoff <- quantile(weights, 1 - top_prop, na.rm = TRUE)
    keep <- nbr_idx[weights >= cutoff]
    for (j in keep) {
      filtered_mat[i, j] <- S_full[i, j]
      filtered_mat[j, i] <- S_full[j, i] # ensure symmetry
    }
  }
  if (i %% 1000 == 0) cat("Processed", i, "nodes for PPI\n")
}

# Membrane gene global edge filtering (same as before)
cat("Filtering top mem_prop membrane-involved edges...\n")
mem_idx <- which(gene_names %in% mem_genes)
mem_mask <- outer(mem_idx, 1:n, function(i, j) TRUE) | outer(1:n, mem_idx, function(i, j) TRUE)
mem_edges <- which(mem_mask, arr.ind = TRUE)
mem_edges <- mem_edges[mem_edges[, 1] <= mem_edges[, 2], , drop = FALSE]

if(nrow(mem_edges) > 0){
  mem_weights <- abs(S_full[mem_edges])
  global_cutoff <- quantile(mem_weights, 1 - mem_prop, na.rm = TRUE)
  keep_mem <- mem_edges[mem_weights >= global_cutoff, , drop = FALSE]
  for (k in seq_len(nrow(keep_mem))) {
    i <- keep_mem[k, 1]
    j <- keep_mem[k, 2]
    filtered_mat[i, j] <- S_full[i, j]
    filtered_mat[j, i] <- S_full[j, i]
  }
}

saveRDS(filtered_mat, file = output_rds)
cat("✓ Filtered network written to", output_rds, "\n")

filtered_sparse <- as(filtered_mat, "dgCMatrix")
writeMM(filtered_sparse, file = output_mm)
writeLines(rownames(filtered_sparse), file = output_genes)
cat("✓ Matrix Market and gene list written\n")

# Cleanup
rm(S_full, filtered_mat, filtered_sparse, ppi_edges, ppi_neighbors, gene_names, gene_idx, mem_genes, mem_idx, mem_mask, mem_edges, mem_weights, keep_mem)
gc()
