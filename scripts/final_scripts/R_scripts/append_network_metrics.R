#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(AnnotationDbi)
})

# -------- args --------
opt <- list(
  make_option("--annot_in",   type="character", help="Annotated epitope CSV (with Gene.Name)"),
  make_option("--betweenness",type="character", default=NA, help="Betweenness TSV"),
  make_option("--degree",     type="character", default=NA, help="Degree TSV"),
  make_option("--impact",     type="character", default=NA, help="Largest component impact TSV"),
  make_option("--strength",   type="character", default=NA, help="Strength TSV (filtered network)"),
  make_option("--wci",        type="character", default=NA, help="WCI TSV (full network)"),
  make_option("--out_csv",    type="character", help="Output CSV")
)
args <- parse_args(OptionParser(option_list = opt))

if (is.null(args$annot_in) || is.null(args$out_csv)) {
  stop("Missing --annot_in or --out_csv", call. = FALSE)
}

# -------- helpers --------
strip_version <- function(x) sub("\\..*$", "", x %||% "")
`%||%` <- function(a, b) if (is.null(a)) b else a

# Force ID column = "gene"; fall back only if truly absent
read_metric <- function(fname, value_col, label) {
  if (is.na(fname) || !nzchar(fname) || !file.exists(fname) || file.size(fname) == 0) {
    message(sprintf("[%s] Missing or empty file -> will add column with NAs.", label))
    return(NULL)
  }
  dt <- data.table::fread(fname)
  id_col <- if ("gene" %in% names(dt)) "gene" else {
    # very defensive fallback
    for (cand in c("Gene","GENE","ensembl","ENSEMBL","id","ID")) {
      if (cand %in% names(dt)) {cand; break}
    }
  }
  if (is.null(id_col) || is.na(id_col)) {
    stop(sprintf("[%s] No gene ID column found in %s (columns=%s)",
                 label, fname, paste(names(dt), collapse=",")))
  }
  if (!(value_col %in% names(dt))) {
    stop(sprintf("[%s] Expected value column '%s' not found in %s (columns=%s)",
                 label, value_col, fname, paste(names(dt), collapse=",")))
  }
  out <- dt[, .(gene = strip_version(get(id_col)), value = get(value_col))]
  # collapse duplicates by median, just in case
  out <- out[, .(value = suppressWarnings(median(as.numeric(value), na.rm = TRUE))), by = gene]
  setnames(out, "value", label)
  out[]
}

# Ensembl → SYMBOL mapping with EnsDb v103 first, org.Hs.eg.db fallback
map_ensg_to_symbol <- function(ens_ids) {
  ens_ids <- unique(strip_version(ens_ids))
  if (length(ens_ids) == 0) return(character())

  # Try EnsDb.Hsapiens.v103 via AnnotationHub (using ANNOTATIONHUB_CACHE if set)
  symbol <- NULL
  try({
    suppressPackageStartupMessages({
      library(AnnotationHub)
      library(ensembldb)
    })
    cache <- Sys.getenv("ANNOTATIONHUB_CACHE", unset = NA)
    ah <- if (!is.na(cache) && nzchar(cache)) AnnotationHub(cache = cache) else AnnotationHub()
    q <- query(ah, "EnsDb.Hsapiens.v103")
    if (length(q) > 0) {
      edb <- q[[1]]
      symbol <- AnnotationDbi::mapIds(edb,
                                      keys   = ens_ids,
                                      keytype= "GENEID",
                                      column = "SYMBOL",
                                      multiVals = "first")
    }
  }, silent = TRUE)

  # Fallback to org.Hs.eg.db if EnsDb path failed
  if (is.null(symbol)) {
    message("[map] Falling back to org.Hs.eg.db")
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    symbol <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                    keys   = ens_ids,
                                    keytype= "ENSEMBL",
                                    column = "SYMBOL",
                                    multiVals = "first")
  }

  # Ensure a named character vector over all input keys
  symbol <- symbol %||% setNames(rep(NA_character_, length(ens_ids)), ens_ids)
  sym_vec <- unname(symbol[strip_version(ens_ids)])
  names(sym_vec) <- ens_ids
  sym_vec
}

# -------- read inputs --------
annot <- data.table::fread(args$annot_in)
if (!("Gene.Name" %in% names(annot))) {
  stop(sprintf("Input annot table lacks 'Gene.Name': %s", args$annot_in))
}

# read each metric (NULL if missing)
m_bet <- read_metric(args$betweenness, "betweenness",               "Net_Betweenness")
m_deg <- read_metric(args$degree,      "degree",                    "Net_Degree")
m_imp <- read_metric(args$impact,      "largest_component_impact",  "Net_Impact")
m_str <- read_metric(args$strength,    "strength",                  "Net_Strength")
m_wci <- read_metric(args$wci,         "wci",                       "Net_WCI")

# combine by gene
metric_list <- Filter(Negate(is.null), list(m_bet, m_deg, m_imp, m_str, m_wci))
if (length(metric_list) == 0) {
  message("[net] No metric files available; writing input with NA columns")
  out <- copy(annot)
  for (nm in c("Net_Betweenness","Net_Degree","Net_Impact","Net_Strength","Net_WCI")) {
    if (!(nm %in% names(out))) out[, (nm) := NA_real_]
  }
  out <- as.data.frame(out)
  data.table::fwrite(out, args$out_csv)
  message(sprintf("Wrote: %s", args$out_csv))
  quit(save = "no")
}

metrics <- Reduce(function(x,y) merge(x, y, by = "gene", all = TRUE), metric_list)

# Map Ensembl → SYMBOL
metrics[, gene_novers := strip_version(gene)]
sym_map <- map_ensg_to_symbol(metrics$gene_novers)
metrics[, SYMBOL := unname(sym_map[gene_novers])]

# Keep one row per SYMBOL (median across duplicates)
metrics_sym <- metrics[!is.na(SYMBOL) & nzchar(SYMBOL),
  .(
    Net_Betweenness = suppressWarnings(median(as.numeric(Net_Betweenness), na.rm=TRUE)),
    Net_Degree      = suppressWarnings(median(as.numeric(Net_Degree), na.rm=TRUE)),
    Net_Impact      = suppressWarnings(median(as.numeric(Net_Impact), na.rm=TRUE)),
    Net_Strength    = suppressWarnings(median(as.numeric(Net_Strength), na.rm=TRUE)),
    Net_WCI         = suppressWarnings(median(as.numeric(Net_WCI), na.rm=TRUE))
  ),
  by = SYMBOL]

# Merge into annotated by Gene.Name (SYMBOL)
setDT(annot)
out <- merge(annot, metrics_sym, by.x = "Gene.Name", by.y = "SYMBOL", all.x = TRUE, sort = FALSE)

# Try to place the new columns after Intogen_* if present
col_order <- names(out)
anchor <- which(col_order %in% c("Intogen_CancerType","Intogen_Cohort","Intogen_Driver_role"))
if (length(anchor)) {
  anchor <- max(anchor)
  keep_front <- col_order[seq_len(anchor)]
  net_cols   <- c("Net_Betweenness","Net_Degree","Net_Impact","Net_Strength","Net_WCI")
  rest       <- setdiff(col_order, c(keep_front, net_cols))
  out <- out[, c(keep_front, intersect(net_cols, names(out)), rest), with = FALSE]
}

# Write
out <- as.data.frame(out, stringsAsFactors = FALSE)
data.table::fwrite(out, args$out_csv)
message(sprintf("Wrote: %s", args$out_csv))
