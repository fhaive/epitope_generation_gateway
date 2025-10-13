#!/usr/bin/env Rscript

# Libraries
suppressMessages({
  library(igraph)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(yaml)
  library(stringr)
})

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

# Load data
mem_genes <- readRDS(mem_genes_file)
HumanNet <- readRDS(humannet_file)
test_network <- readRDS(input_rds)
config <- read_yaml(config_file)

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

# -------- 2. Per-node top-20% PPI edges filtering --------
top_prop <- 0.20
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

# -------- 3. Add membrane-touched edges: top N by abs(weight), N = 5% of unfiltered PPI edges --------
mem_row_idx <- which(rownames(test_network) %in% mem_genes)
mem_col_idx <- which(colnames(test_network) %in% mem_genes)
mat_idx_mem_row <- expand.grid(i = mem_row_idx, j = seq_len(nrow(test_network)))
mat_idx_mem_col <- expand.grid(i = seq_len(nrow(test_network)), j = mem_col_idx)
mem_edges_idx <- rbind(mat_idx_mem_row, mat_idx_mem_col)
mem_edges_idx <- mem_edges_idx[mem_edges_idx$i != mem_edges_idx$j, ]
a <- pmin(mem_edges_idx$i, mem_edges_idx$j)
b <- pmax(mem_edges_idx$i, mem_edges_idx$j)
edges_sorted <- cbind(a, b)
edges_unique <- unique(edges_sorted)
colnames(edges_unique) <- c("i", "j")
mem_edges_idx <- edges_unique

all_genes <- rownames(test_network)
mem_edge_df <- data.frame(
  from = all_genes[mem_edges_idx[, 1]],
  to = all_genes[mem_edges_idx[, 2]],
  weight = test_network[mem_edges_idx]
)
mem_edge_df <- mem_edge_df[mem_edge_df$from != mem_edge_df$to, ]

ppi_top_df <- igraph::as_data_frame(g_ppi_top, what = "edges")
ppi_top_edge_keys <- paste(pmin(ppi_top_df$from, ppi_top_df$to), pmax(ppi_top_df$from, ppi_top_df$to), sep = "_")
mem_edge_keys <- paste(pmin(mem_edge_df$from, mem_edge_df$to), pmax(mem_edge_df$from, mem_edge_df$to), sep = "_")
new_mem_edges <- mem_edge_df[!mem_edge_keys %in% ppi_top_edge_keys, ]

ppi_weighted_df <- igraph::as_data_frame(g_ppi_weighted, what = "edges")
n_to_add <- ceiling(0.05 * nrow(ppi_weighted_df))
cat(sprintf("Will add %d membrane-touched edges (5%% of g_ppi_weighted edges = %d)\n", n_to_add, nrow(ppi_weighted_df)))
top_N_mem_edges <- new_mem_edges[order(abs(new_mem_edges$weight), decreasing = TRUE), ][1:n_to_add, ]

all_vertex_names <- unique(c(
  V(g_ppi_top)$name,
  mem_genes,
  top_N_mem_edges$from,
  top_N_mem_edges$to
))

all_edges <- rbind(
  ppi_top_df[, c("from", "to", "weight")],
  top_N_mem_edges[, c("from", "to", "weight")]
)
all_edge_keys <- paste(pmin(all_edges$from, all_edges$to), pmax(all_edges$from, all_edges$to), sep = "_")
all_edges <- all_edges[!duplicated(all_edge_keys), ]

final_graph <- igraph::graph_from_data_frame(all_edges, directed = FALSE, vertices = all_vertex_names)

cat(sprintf("Final graph: %d nodes, %d edges\n", vcount(final_graph), ecount(final_graph)))

# -------- 4. Save outputs (RDS, Matrix Market, gene list) --------
saveRDS(final_graph, file = output_rds)

# Convert to adjacency matrix for Matrix Market
adjmat <- as_adjacency_matrix(final_graph, attr = "weight", sparse = TRUE, type = "both")
writeMM(adjmat, output_mm)
writeLines(rownames(adjmat), output_genes)

cat("✓ Done!\n")
