library(yaml)
library(stringr)
library(Matrix)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
input_rds     <- args[1]
output_rds    <- args[2]
output_mm     <- args[3]
output_genes  <- args[4]
config_file   <- args[5]
humannet_file <- args[6]
mem_genes_file <- args[7]

cat("Input RDS:", input_rds, "\n")
cat("Output RDS:", output_rds, "\n")
cat("Matrix Market:", output_mm, "\n")
cat("Gene list:", output_genes, "\n")
cat("Config File:", config_file, "\n")
cat("HumanNet File:", humannet_file, "\n")
cat("Membrane Genes File:", mem_genes_file, "\n")

# --- Config ---
config   <- read_yaml(config_file)
TOP_PROP <- config$top_prop   # e.g. 0.20
MEM_PROP <- config$mem_prop   # e.g. 0.05

# --- Reference sets ---
humannet <- readRDS(humannet_file)
ppi_keys <- humannet %>%
  transmute(g1 = pmin(ensembl1, ensembl2),
            g2 = pmax(ensembl1, ensembl2),
            edge_key = paste(g1, g2, sep = "_")) %>%
  pull(edge_key) %>% unique()
mem_genes <- as.character(readRDS(mem_genes_file))

# --- Read network ---
samp <- str_remove(basename(input_rds), "\\.rds$")
cat("Processing sample:", samp, "\n")
S <- readRDS(input_rds)      # dense or sparse
G <- rownames(S); n <- length(G)

# --- PPI edge collection (summary-based, robust) ---
cat("Collecting PPI edge weights...\n")
sm    <- summary(S)
regs  <- rownames(S)[sm$i]
tars  <- colnames(S)[sm$j]
wts   <- sm$x

edge_key <- ifelse(regs < tars,
                   paste(regs, tars, sep = "_"),
                   paste(tars, regs, sep = "_"))

ppi_idx  <- which(edge_key %in% ppi_keys)
regs_ppi <- regs[ppi_idx]
tars_ppi <- tars[ppi_idx]
wts_ppi  <- wts[ppi_idx]

# --- Top-X% of PPI edges ---
if (length(wts_ppi) > 0) {
  cutoff20 <- quantile(abs(wts_ppi), 1 - TOP_PROP)
  top_idx  <- which(abs(wts_ppi) >= cutoff20)
} else {
  top_idx <- integer(0)
}

# --- Membrane edge filtering: top MEM_PROP globally ---
cat("Filtering top membrane edges globally...\n")
mem_idx <- which(G %in% mem_genes)
mem_edges <- which(outer(mem_idx, seq_len(n), Vectorize(function(i, j) TRUE)) |
                   outer(seq_len(n), mem_idx, Vectorize(function(i, j) TRUE)), arr.ind = TRUE)
mem_edges <- mem_edges[mem_edges[,1] < mem_edges[,2], , drop = FALSE]
if (nrow(mem_edges) > 0) {
  mem_weights <- abs(S[mem_edges])
  global_cutoff <- quantile(mem_weights, 1 - MEM_PROP, na.rm = TRUE)
  keep_mem <- mem_edges[mem_weights >= global_cutoff, , drop = FALSE]
} else {
  keep_mem <- matrix(nrow = 0, ncol = 2)
}

# --- Build filtered matrix ---
keep_idx <- unique(top_idx)
if (length(keep_idx) == 0 && nrow(keep_mem) == 0) {
  filtered_mat <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                              dims = dim(S), dimnames = dimnames(S))
} else {
  # From PPI top edges
  ii <- sm$i[ppi_idx][keep_idx]
  jj <- sm$j[ppi_idx][keep_idx]
  xx <- wts_ppi[keep_idx]

  # From membrane edges
  mem_i <- keep_mem[, 1]
  mem_j <- keep_mem[, 2]
  mem_x <- S[as.matrix(keep_mem)]

  # Build sparse matrix (PPI + membrane, both directions)
  filtered_mat <- sparseMatrix(
    i = c(ii, jj, mem_i, mem_j),
    j = c(jj, ii, mem_j, mem_i),
    x = c(xx, xx, mem_x, mem_x),
    dims = dim(S),
    dimnames = dimnames(S)
  )
}

# --- Output ---
saveRDS(filtered_mat, output_rds)
writeMM(filtered_mat, output_mm)
writeLines(rownames(filtered_mat), output_genes)
cat("✓ Done for", samp, "\n")

# --- Cleanup ---
rm(S, filtered_mat, sm, regs, tars, wts, edge_key,
   ppi_idx, regs_ppi, tars_ppi, wts_ppi, top_idx, ii, jj, xx,
   mem_i, mem_j, mem_x, mem_edges, mem_weights, keep_mem)
gc()
