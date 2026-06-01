#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
})

COUNT_FILE <- "tcga_skcm_expression_matrix/tcga_skcm_star_counts_unstranded_primary_tumor_one_sample_per_patient.tsv"

OUT_DIR <- "tcga_skcm_expression_matrix"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_VST <- file.path(OUT_DIR, "tcga_skcm_vst_mat.csv")
OUT_VST_TSV <- file.path(OUT_DIR, "tcga_skcm_vst_mat.tsv")
OUT_FILTERED_COUNTS <- file.path(OUT_DIR, "tcga_skcm_counts_filtered_for_vst.tsv")
OUT_QC <- file.path(OUT_DIR, "tcga_skcm_vst_qc_summary.tsv")

message2 <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ...)
  cat("\n")
  flush.console()
}

message2("Reading count matrix: ", COUNT_FILE)

counts_dt <- fread(COUNT_FILE)

gene_col <- names(counts_dt)[1]
genes <- counts_dt[[gene_col]]

counts <- as.matrix(counts_dt[, -1, with = FALSE])
rownames(counts) <- genes
storage.mode(counts) <- "integer"

message2("Raw matrix: ", nrow(counts), " genes x ", ncol(counts), " samples")

# Basic low-count filtering for DESeq2 stability.
# Keep genes with at least 10 counts in at least 10 samples.
keep <- rowSums(counts >= 10) >= 10
counts_filt <- counts[keep, , drop = FALSE]

message2("Filtered matrix for VST: ", nrow(counts_filt), " genes x ", ncol(counts_filt), " samples")

coldata <- data.frame(
  sample = colnames(counts_filt),
  row.names = colnames(counts_filt)
)

dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData = coldata,
  design = ~ 1
)

message2("Estimating size factors")
dds <- estimateSizeFactors(dds)

message2("Running vst")
vst_obj <- vst(dds, blind = TRUE)
vst_mat <- assay(vst_obj)

# Write as CSV matching original style: first column rownames when read.csv(row.names=1)
vst_out <- data.frame(gene_id = rownames(vst_mat), vst_mat, check.names = FALSE)
fwrite(vst_out, OUT_VST_TSV, sep = "\t")

# CSV with first column as row names-like field
fwrite(vst_out, OUT_VST)

counts_out <- data.frame(gene_id = rownames(counts_filt), counts_filt, check.names = FALSE)
fwrite(counts_out, OUT_FILTERED_COUNTS, sep = "\t")

qc <- data.table(
  raw_genes = nrow(counts),
  filtered_genes = nrow(counts_filt),
  samples = ncol(counts_filt),
  min_library_size = min(colSums(counts)),
  median_library_size = median(colSums(counts)),
  max_library_size = max(colSums(counts))
)

fwrite(qc, OUT_QC, sep = "\t")

message2("Done")
message2("VST CSV: ", OUT_VST)
message2("VST TSV: ", OUT_VST_TSV)
message2("Filtered counts: ", OUT_FILTERED_COUNTS)
message2("QC: ", OUT_QC)

print(qc)
