#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

VST_FILE <- "tcga_skcm_expression_matrix/tcga_skcm_vst_mat.tsv"

original_30 <- c(
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

dt <- fread(VST_FILE, nrows = 1)
cols <- names(dt)
sample_cols <- setdiff(cols, "gene_id")

present <- original_30 %in% sample_cols

res <- data.table(
  sample = original_30,
  present = present
)

fwrite(res, "tcga_skcm_expression_matrix/original_30_presence_in_full_tcga_vst.tsv", sep = "\t")

cat("Original 30 present:", sum(present), "/ 30\n")
if (any(!present)) {
  cat("Missing:\n")
  print(res[present == FALSE])
}

cat("Full TCGA-SKCM VST sample count:", length(sample_cols), "\n")
