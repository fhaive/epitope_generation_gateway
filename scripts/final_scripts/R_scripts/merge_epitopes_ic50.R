#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(yaml)
  library(optparse)
  library(purrr)
})

# ---------- CLI ----------
opt_list <- list(
  make_option("--somatic_mhci",  type="character", default=NA),
  make_option("--somatic_mhcii", type="character", default=NA),
  make_option("--fusion_mhci",   type="character", default=NA),
  make_option("--fusion_mhcii",  type="character", default=NA),
  make_option("--splice_mhci",   type="character", default=NA),
  make_option("--splice_mhcii",  type="character", default=NA),
  make_option("--columns_yaml",  type="character", default=NA),
  make_option("--out_csv",       type="character", default="epitopes_merged_ic50.csv")
)
opt <- parse_args(OptionParser(option_list = opt_list))

# ---------- helpers ----------
exists_nonempty <- function(p) {
  !is.na(p) && nzchar(p) && file.exists(p) && isTRUE(file.info(p)$size > 0)
}

read_tab <- function(path) {
  if (!exists_nonempty(path)) return(data.frame())
  as.data.frame(fread(path, sep = "auto", header = TRUE, data.table = FALSE, showProgress = FALSE))
}

pull_or_na <- function(df, col) {
  if (!is.na(col) && nzchar(col) && col %in% names(df)) df[[col]] else rep(NA, nrow(df))
}

# Load config (exact names only)
if (!exists_nonempty(opt$columns_yaml)) stop("Missing or unreadable --columns_yaml")
cfg <- yaml::read_yaml(opt$columns_yaml)

# Get ordered list of column names: column1, column2, ...
get_cols <- function(source, hla_class) {
  blk <- cfg$required_columns[[source]][[hla_class]]
  if (is.null(blk)) return(character(0))
  idx <- order(as.integer(gsub("^column", "", names(blk))))
  unname(as.character(blk[idx]))
}

# Return extras beyond the first N required columns for a source/hla
get_extras <- function(source, hla_class, required_n) {
  blk <- cfg$required_columns[[source]][[hla_class]]
  if (is.null(blk)) return(character(0))
  nums <- as.integer(gsub("^column", "", names(blk)))
  keep <- nums > required_n
  if (!any(keep)) return(character(0))
  ord <- order(nums[keep])
  unname(as.character(blk[keep][ord]))
}

standard_cols <- c(
  "Gene.Name","HLA.Allele","Peptide.Length","MT.Epitope.Seq","WT.Epitope.Seq",
  "Median.MT.IC50.Score","Mutation.Source","HLA.Class","Fusion_Genes"
)

append_extras <- function(out, df, extras) {
  if (!length(extras)) return(out)
  extra_keep <- setdiff(intersect(extras, names(df)), standard_cols)
  if (!length(extra_keep)) return(out)
  bind_cols(out, df[, extra_keep, drop = FALSE])
}

# ------- per-source builders (STRICT: use only configured names) -------
build_somatic <- function(df, cols, extras, hlabel) {
  if (!nrow(df)) return(tibble())
  need <- 6
  cols <- c(cols, rep(NA_character_, max(0, need - length(cols))))[seq_len(need)]
  out <- tibble(
    Gene.Name              = pull_or_na(df, cols[1]),
    HLA.Allele             = pull_or_na(df, cols[2]),
    Peptide.Length         = suppressWarnings(as.integer(pull_or_na(df, cols[3]))),
    `MT.Epitope.Seq`       = pull_or_na(df, cols[4]),
    `WT.Epitope.Seq`       = pull_or_na(df, cols[5]),
    `Median.MT.IC50.Score` = suppressWarnings(as.numeric(pull_or_na(df, cols[6]))),
    Mutation.Source        = "Somatic",
    HLA.Class              = hlabel,
    Fusion_Genes           = NA_character_
  )
  append_extras(out, df, extras)
}

build_fusion <- function(df, cols, extras, hlabel) {
  # now expects 5 required: gene, hla, peptide.length, epitope.seq, median.ic50
  if (!nrow(df)) return(tibble())
  need <- 5
  cols <- c(cols, rep(NA_character_, max(0, need - length(cols))))[seq_len(need)]

  mt_seq  <- pull_or_na(df, cols[4])
  pep_len <- pull_or_na(df, cols[3])
  if (all(is.na(pep_len)) && !all(is.na(mt_seq))) pep_len <- nchar(mt_seq)

  out <- tibble(
    Gene.Name              = pull_or_na(df, cols[1]),
    HLA.Allele             = pull_or_na(df, cols[2]),
    Peptide.Length         = suppressWarnings(as.integer(pep_len)),
    `MT.Epitope.Seq`       = mt_seq,
    `WT.Epitope.Seq`       = NA_character_,
    `Median.MT.IC50.Score` = suppressWarnings(as.numeric(pull_or_na(df, cols[5]))),
    Mutation.Source        = "Fusion",
    HLA.Class              = hlabel,
    Fusion_Genes           = pull_or_na(df, cols[1])
  )

  out <- append_extras(out, df, extras)

  out %>%
    mutate(Gene.Name = str_split(Gene.Name, "_")) %>%
    unnest(Gene.Name)
}

build_splicing <- function(df, cols, extras, hlabel) {
  if (!nrow(df)) return(tibble())
  need <- 5
  cols <- c(cols, rep(NA_character_, max(0, need - length(cols))))[seq_len(need)]

  mt_seq  <- pull_or_na(df, cols[4])
  pep_len <- pull_or_na(df, cols[3])
  if (all(is.na(pep_len)) && !all(is.na(mt_seq))) pep_len <- nchar(mt_seq)

  out <- tibble(
    Gene.Name              = pull_or_na(df, cols[1]),
    HLA.Allele             = pull_or_na(df, cols[2]),
    Peptide.Length         = suppressWarnings(as.integer(pep_len)),
    `MT.Epitope.Seq`       = mt_seq,
    `WT.Epitope.Seq`       = NA_character_,
    `Median.MT.IC50.Score` = suppressWarnings(as.numeric(pull_or_na(df, cols[5]))),
    Mutation.Source        = "Splice Variant",
    HLA.Class              = hlabel,
    Fusion_Genes           = NA_character_
  )

  append_extras(out, df, extras)
}

# ---------- read inputs ----------
som_mhci  <- read_tab(opt$somatic_mhci)
som_mhcii <- read_tab(opt$somatic_mhcii)
fus_mhci  <- read_tab(opt$fusion_mhci)
fus_mhcii <- read_tab(opt$fusion_mhcii)
spl_mhci  <- read_tab(opt$splice_mhci)
spl_mhcii <- read_tab(opt$splice_mhcii)

# ---------- build frames using ONLY configured names ----------
som_mhci_cols   <- get_cols("somatic",  "MHCI")
som_mhcii_cols  <- get_cols("somatic",  "MHCII")
fus_mhci_cols   <- get_cols("fusion",   "MHCI")
fus_mhcii_cols  <- get_cols("fusion",   "MHCII")
spl_mhci_cols   <- get_cols("splicing", "MHCI")
spl_mhcii_cols  <- get_cols("splicing", "MHCII")

som_mhci_extra  <- get_extras("somatic",  "MHCI",  6)
som_mhcii_extra <- get_extras("somatic",  "MHCII", 6)
fus_mhci_extra  <- get_extras("fusion",   "MHCI",  5)  # updated required_n
fus_mhcii_extra <- get_extras("fusion",   "MHCII", 5)  # updated required_n
spl_mhci_extra  <- get_extras("splicing", "MHCI",  5)
spl_mhcii_extra <- get_extras("splicing", "MHCII", 5)

lst <- list(
  build_somatic (som_mhci,  som_mhci_cols,  som_mhci_extra,  "MHC Class I"),
  build_somatic (som_mhcii, som_mhcii_cols, som_mhcii_extra, "MHC Class II"),
  build_fusion  (fus_mhci,  fus_mhci_cols,  fus_mhci_extra,  "MHC Class I"),
  build_fusion  (fus_mhcii, fus_mhcii_cols, fus_mhcii_extra, "MHC Class II"),
  build_splicing(spl_mhci,  spl_mhci_cols,  spl_mhci_extra,  "MHC Class I"),
  build_splicing(spl_mhcii, spl_mhcii_cols, spl_mhcii_extra, "MHC Class II")
)

res <- lst %>% discard(~ is.null(.x) || !nrow(.x))

if (!length(res)) {
  empty <- setNames(data.frame(matrix(nrow = 0, ncol = length(standard_cols))), standard_cols)
  fwrite(empty, opt$out_csv)
  quit(save="no", status=0)
}

merged <- bind_rows(res) %>%
  distinct() %>%
  arrange(HLA.Class, Gene.Name)

extra_cols_all <- setdiff(names(merged), standard_cols)
merged <- merged[, c(standard_cols, extra_cols_all), drop = FALSE]

fwrite(merged, opt$out_csv)
cat("[merge_epitopes_ic50] Wrote:", opt$out_csv, "\n")
