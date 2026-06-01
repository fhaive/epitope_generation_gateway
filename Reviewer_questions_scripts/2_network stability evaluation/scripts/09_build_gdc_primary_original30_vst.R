#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
})

COUNT_FILE <- "tcga_skcm_expression_matrix/tcga_skcm_star_counts_unstranded_primary_tumor_one_sample_per_patient.tsv"

OUT_DIR <- "gdc_primary_original30_expression_matrix"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_COUNTS <- file.path(OUT_DIR, "gdc_primary_original30_counts.tsv")
OUT_FILTERED_COUNTS <- file.path(OUT_DIR, "gdc_primary_original30_counts_filtered_for_vst.tsv")
OUT_VST_TSV <- file.path(OUT_DIR, "gdc_primary_original30_vst_mat.tsv")
OUT_VST_CSV <- file.path(OUT_DIR, "gdc_primary_original30_vst_mat.csv")
OUT_QC <- file.path(OUT_DIR, "gdc_primary_original30_vst_qc_summary.tsv")

TARGETS <- c(
  "Sample_TCGA-BF-A1PV",
  "Sample_TCGA-BF-A3DJ",
  "Sample_TCGA-BF-A3DL",
  "Sample_TCGA-BF-A3DM",
  "Sample_TCGA-BF-A3DN",
  "Sample_TCGA-BF-A5EO",
  "Sample_TCGA-BF-A5EP",
  "Sample_TCGA-BF-A5EQ",
  "Sample_TCGA-BF-A5ER",
  "Sample_TCGA-BF-A5ES",
  "Sample_TCGA-BF-AAOU",
  "Sample_TCGA-BF-AAOX",
  "Sample_TCGA-BF-AAP1",
  "Sample_TCGA-BF-AAP2",
  "Sample_TCGA-BF-AAP4",
  "Sample_TCGA-BF-AAP6",
  "Sample_TCGA-BF-AAP7",
  "Sample_TCGA-BF-AAP8",
  "Sample_TCGA-D3-A5GT",
  "Sample_TCGA-D9-A3Z4",
  "Sample_TCGA-D9-A4Z2",
  "Sample_TCGA-D9-A4Z3",
  "Sample_TCGA-EB-A1NK",
  "Sample_TCGA-EB-A3HV",
  "Sample_TCGA-EB-A3XB",
  "Sample_TCGA-EB-A3XC",
  "Sample_TCGA-EB-A3XD",
  "Sample_TCGA-EB-A3XE",
  "Sample_TCGA-EB-A3XF",
  "Sample_TCGA-EB-A41A"
)

message2 <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ...)
  cat("\n")
  flush.console()
}

message2("Reading primary-tumor count matrix: ", COUNT_FILE)

dt <- fread(COUNT_FILE)
gene_col <- names(dt)[1]

missing_targets <- setdiff(TARGETS, names(dt))
if (length(missing_targets) > 0) {
  stop("Missing target columns: ", paste(missing_targets, collapse = ", "))
}

sub_dt <- dt[, c(gene_col, TARGETS), with = FALSE]
fwrite(sub_dt, OUT_COUNTS, sep = "\t")

genes <- sub_dt[[gene_col]]
counts <- as.matrix(sub_dt[, -1, with = FALSE])
rownames(counts) <- genes
storage.mode(counts) <- "integer"

message2("Original-30 raw count matrix: ", nrow(counts), " genes x ", ncol(counts), " samples")

# Same simple low-count rule used for the n=103 public-GDC benchmark.
keep <- rowSums(counts >= 10) >= 10
counts_filt <- counts[keep, , drop = FALSE]

message2("Filtered for VST: ", nrow(counts_filt), " genes x ", ncol(counts_filt), " samples")

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

vst_out <- data.frame(gene_id = rownames(vst_mat), vst_mat, check.names = FALSE)
fwrite(vst_out, OUT_VST_TSV, sep = "\t")
fwrite(vst_out, OUT_VST_CSV)

counts_out <- data.frame(gene_id = rownames(counts_filt), counts_filt, check.names = FALSE)
fwrite(counts_out, OUT_FILTERED_COUNTS, sep = "\t")

qc <- data.table(
  matrix = "GDC-primary original 30 only",
  raw_genes = nrow(counts),
  filtered_genes = nrow(counts_filt),
  samples = ncol(counts_filt),
  min_library_size = min(colSums(counts)),
  median_library_size = median(colSums(counts)),
  max_library_size = max(colSums(counts))
)

fwrite(qc, OUT_QC, sep = "\t")

message2("Done")
print(qc)

cat("\nFiles:\n")
cat(OUT_COUNTS, "\n")
cat(OUT_FILTERED_COUNTS, "\n")
cat(OUT_VST_TSV, "\n")
cat(OUT_QC, "\n")
