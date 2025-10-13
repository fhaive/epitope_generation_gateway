#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(optparse)
})

opt_list <- list(
  make_option("--merged_in", type="character"),
  make_option("--depmap",    type="character"),
  make_option("--intogen",   type="character"),
  make_option("--out_csv",   type="character")
)
opt <- parse_args(OptionParser(option_list = opt_list))

stopifnot(file.exists(opt$merged_in),
          file.exists(opt$depmap),
          file.exists(opt$intogen))

# Read inputs ---------------------------------------------------------------
merged <- as.data.frame(fread(opt$merged_in))
dep    <- as.data.frame(fread(opt$depmap))
intg   <- as.data.frame(fread(opt$intogen))

# Clean & standardize keys --------------------------------------------------
# merged has gene column "Gene.Name" (dot-style)
merged <- merged %>%
  mutate(.SYMBOL_KEY = toupper(trimws(Gene.Name)))

# depmap: keep useful cols and rename Chronos median -> requested name
# (Filter out odd rows like "V1" or blank genes if present)
dep <- dep %>%
  filter(!is.na(Gene), Gene != "", Gene != "V1") %>%
  transmute(
    .SYMBOL_KEY = toupper(trimws(Gene)),
    Depmap_survivability_score = suppressWarnings(as.numeric(PanCancer_Median_Chronos))
  )

# intogen: keep and rename requested fields
intg <- intg %>%
  transmute(
    .SYMBOL_KEY          = toupper(trimws(SYMBOL)),
    Intogen_Driver_role  = ROLE,
    Intogen_Cohort       = COHORT,
    Intogen_CancerType   = CANCER_TYPE
  )

# Join (left) ---------------------------------------------------------------
out <- merged %>%
  left_join(dep,  by = ".SYMBOL_KEY") %>%
  left_join(intg, by = ".SYMBOL_KEY") %>%
  select(-.SYMBOL_KEY)

# Write ---------------------------------------------------------------------
fwrite(out, opt$out_csv)
cat("[annotate_epitopes] Wrote:", opt$out_csv, "\n")
