library(DESeq2)
library(tidyverse)
library(data.table)
library(purrr)
library(igraph)


# Get command-line arguments: output file followed by count files
args <- commandArgs(trailingOnly = TRUE)
output_file_rel <- args[1]  # Relative output file path
count_files_rel <- args[-1]  # Relative input count file paths

# Check if count files are provided
if (length(count_files_rel) == 0) {
  stop("No count files provided.")
}

# Get project root from current working directory (set by Snakemake)
project_root <- getwd()

# Debug: Print project root and relative paths
cat("Project root:", project_root, "\n")
cat("Relative output file:", output_file_rel, "\n")
cat("Relative count files:", count_files_rel, "\n")

# Resolve relative paths to absolute paths
output_file <- normalizePath(file.path(project_root, output_file_rel), mustWork = FALSE)
count_files <- sapply(count_files_rel, function(f) normalizePath(file.path(project_root, f)))

# Debug: Print resolved absolute paths
cat("Absolute output file:", output_file, "\n")
cat("Absolute count files:", count_files, "\n")

# Ensure output directory exists
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Extract sample names from file paths
sample_names <- sapply(count_files, function(f) {
  basename <- basename(f)
  sub("_CancerRNA.*$", "", basename)
})

# Read and merge count files
counts_list <- lapply(count_files, function(f) {
  if (!file.exists(f)) stop("File does not exist: ", f)
  dt <- fread(f, skip = 1)
  sample_id <- sub("_CancerRNA.*$", "", basename(f))
  setnames(dt, old = tail(names(dt), 1), new = sample_id)
  dt[, .SD, .SDcols = c("Geneid", sample_id)]
})

merged_counts <- purrr::reduce(counts_list, full_join, by = "Geneid")

# Strip Ensembl version tags and merge duplicates
count_mat <- as.matrix(merged_counts[,-1])
rownames(count_mat) <- merged_counts$Geneid
storage.mode(count_mat) <- "integer"

gene_ids <- sub("\\..*$", "", rownames(count_mat))  # Remove ".10"
count_mat <- rowsum(count_mat, group = gene_ids)   # Sum duplicates

# DESeq2 normalization (VST)
colData <- data.frame(condition = rep("all", ncol(count_mat)),
                      row.names = colnames(count_mat))
if (!all(colnames(count_mat) %in% sample_names)) {
  stop("Sample names mismatch.")
}

dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = colData,
                              design = ~ 1)

# Dynamic threshold: floor(total_samples / 5)
total_samples <- ncol(count_mat)
threshold <- floor(total_samples / 5)
dds75 <- dds[rowSums(counts(dds) >= 15) >= threshold, ]  # Keep genes with counts
dds_vst <- vst(dds75)
vst_mat <- assay(dds_vst)

# Save output
write.csv(vst_mat, output_file, row.names = TRUE)
