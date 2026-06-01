#!/usr/bin/env Rscript

suppressMessages({
  library(igraph)
  library(Matrix)
  library(dplyr)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
  stop("Usage: Filter_PPI_Sample_Networks.R <input_rds> <output_rds> <output_mm> <output_genes> <config_file> <humannet_file>")
}

input_rds    <- args[1]
output_rds   <- args[2]
output_mm    <- args[3]
output_genes <- args[4]
config_file  <- args[5]
humannet_file <- args[6]

cat("Input RDS:", input_rds, "\n")
cat("Output PPI-filtered RDS:", output_rds, "\n")
cat("Matrix Market:", output_mm, "\n")
cat("Gene list:", output_genes, "\n")
cat("Config File:", config_file, "\n")
cat("HumanNet File:", humannet_file, "\n")

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_mm), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_genes), recursive = TRUE, showWarnings = FALSE)

test_network <- readRDS(input_rds)
HumanNet <- readRDS(humannet_file)
config <- yaml::read_yaml(config_file)

top_prop <- 0.20
if (!is.null(config$ppi_filter) && !is.null(config$ppi_filter$top_prop)) {
  top_prop <- as.numeric(config$ppi_filter$top_prop)
} else if (!is.null(config$top_prop)) {
  top_prop <- as.numeric(config$top_prop)
}

if (is.na(top_prop) || top_prop <= 0 || top_prop > 1) {
  stop("ppi_filter.top_prop must be > 0 and <= 1")
}

cat(sprintf("Using PPI per-node top_prop = %.4f\n", top_prop))

# -------- 1. Build weighted PPI graph from HumanNet for this sample --------

ppi_edges <- HumanNet %>%
  select(from = ensembl1, to = ensembl2) %>%
  filter(from != to)

g_ppi <- igraph::graph_from_data_frame(ppi_edges, directed = FALSE)

edge_df <- igraph::as_data_frame(g_ppi, what = "edges")

idx_from <- match(edge_df$from, rownames(test_network))
idx_to   <- match(edge_df$to, colnames(test_network))

edge_df$weight <- NA_real_

valid1 <- !is.na(idx_from) & !is.na(idx_to)
edge_df$weight[valid1] <- test_network[cbind(idx_from[valid1], idx_to[valid1])]

missing <- is.na(edge_df$weight)

if (any(missing)) {
  idx_from2 <- match(edge_df$to[missing], rownames(test_network))
  idx_to2   <- match(edge_df$from[missing], colnames(test_network))
  valid2 <- !is.na(idx_from2) & !is.na(idx_to2)
  edge_df$weight[which(missing)[valid2]] <- test_network[cbind(idx_from2[valid2], idx_to2[valid2])]
}

E(g_ppi)$weight <- edge_df$weight
g_ppi_weighted <- delete_edges(g_ppi, which(is.na(E(g_ppi)$weight)))

cat(sprintf("Weighted PPI graph: %d nodes, %d edges\n", vcount(g_ppi_weighted), ecount(g_ppi_weighted)))

# -------- 2. Per-node top-proportion PPI edge filtering --------

start_time <- Sys.time()
top_eids <- c()
n_nodes <- vcount(g_ppi_weighted)

for (i in seq_len(n_nodes)) {
  v <- V(g_ppi_weighted)[i]
  inc_e <- incident(g_ppi_weighted, v, mode = "all")
  wts <- abs(E(g_ppi_weighted)[inc_e]$weight)

  if (length(wts) > 0) {
    n_keep <- max(1, ceiling(length(wts) * top_prop))
    top_idx <- order(wts, decreasing = TRUE)[seq_len(n_keep)]
    top_eids <- c(top_eids, inc_e[top_idx])
  }

  if (i %% 1000 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat(sprintf("Processed %d of %d nodes in %.1f seconds\n", i, n_nodes, elapsed))
    flush.console()
  }
}

top_eids <- unique(top_eids)
g_ppi_top <- subgraph.edges(g_ppi_weighted, eids = top_eids, delete.vertices = FALSE)

cat(sprintf("PPI-only filtered graph: %d nodes, %d edges\n", vcount(g_ppi_top), ecount(g_ppi_top)))

# -------- 3. Save outputs --------

saveRDS(g_ppi_top, file = output_rds)

adjmat <- as_adjacency_matrix(g_ppi_top, attr = "weight", sparse = TRUE, type = "both")
writeMM(adjmat, output_mm)
writeLines(rownames(adjmat), output_genes)

cat("✓ Done PPI-only filtering\n")
