#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
  library(msigdbr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- CLI ----
opt <- OptionParser()
opt <- add_option(opt, "--input_csv",  type="character", help="Input CSV (has Gene.Name, Borda_* columns)")
opt <- add_option(opt, "--out_csv",    type="character", help="Output CSV with GO_* and Cancer_Hallmark added")
opt <- add_option(opt, "--species",    type="character", default="Homo sapiens", help="Species for msigdbr (default: Homo sapiens)")
opt <- add_option(opt, "--include_iea",action="store_true", default=TRUE, help="Include IEA evidence for GO (default: TRUE)")
opt <- add_option(opt, "--hallmark",   type="character", default="msigdbr",
                  help="Hallmark source: 'msigdbr' (default) or path to a .gmt file")

args <- parse_args(opt)
stopifnot(!is.null(args$input_csv), !is.null(args$out_csv))

# ---- Load input ----
dt <- fread(args$input_csv)
if (!"Gene.Name" %in% names(dt)) {
  stop("Input file must contain a 'Gene.Name' column.")
}
symbols <- unique(na.omit(as.character(dt$Gene.Name)))

# ---- GO annotations (MF/BP/CC) from org.Hs.eg.db + GO.db ----
message("[GO] Retrieving SYMBOL -> GO mappings...")
go_cols <- c("GO","ONTOLOGY","EVIDENCE")
go_map <- AnnotationDbi::select(org.Hs.eg.db, keys = symbols, keytype = "SYMBOL", columns = go_cols)
setDT(go_map)

# Filter evidence if requested
if (!isTRUE(args$include_iea) && "EVIDENCE" %in% names(go_map)) {
  go_map <- go_map[EVIDENCE != "IEA" | is.na(EVIDENCE)]
}

# Map GO IDs -> term names
go2term <- AnnotationDbi::select(GO.db,
                                 keys    = unique(na.omit(go_map$GO)),
                                 keytype = "GOID",
                                 columns = "TERM")
setDT(go2term)
setnames(go2term, c("GOID","TERM"), c("GO","TERM"))
go_map <- merge(go_map, go2term, by = "GO", all.x = TRUE)

# Collapse per ontology
collapse_terms <- function(x) {
  x <- sort(unique(na.omit(x)))
  if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}

go_mf <- go_map[ONTOLOGY == "MF", .(GO_MF = collapse_terms(TERM)), by = SYMBOL]
go_bp <- go_map[ONTOLOGY == "BP", .(GO_BP = collapse_terms(TERM)), by = SYMBOL]
go_cc <- go_map[ONTOLOGY == "CC", .(GO_CC = collapse_terms(TERM)), by = SYMBOL]

# Merge GO columns onto dt by gene symbol
setDT(dt)
if ("SYMBOL" %in% names(dt)) setnames(dt, "SYMBOL", "Gene.Name") # just in case
dt <- merge(dt, go_mf, by.x = "Gene.Name", by.y = "SYMBOL", all.x = TRUE, sort = FALSE)
dt <- merge(dt, go_bp, by.x = "Gene.Name", by.y = "SYMBOL", all.x = TRUE, sort = FALSE)
dt <- merge(dt, go_cc, by.x = "Gene.Name", by.y = "SYMBOL", all.x = TRUE, sort = FALSE)

# ---- Cancer Hallmark annotations (MSigDB Hallmark H) ----
to_upper <- function(x) {
  y <- as.character(x); y[!is.na(y)] <- toupper(y[!is.na(y)]); y
}

hk_tbl <- NULL
if (identical(tolower(args$hallmark), "msigdbr")) {
  message("[Hallmark] Using msigdbr Hallmark sets (category H)...")
  hm <- msigdbr(species = args$species, category = "H")
  if (nrow(hm)) {
    hk_tbl <- as.data.table(hm[, c("gs_name","gene_symbol")])
    setnames(hk_tbl, c("gs_name","gene_symbol"), c("gs_name","symbol"))
  }
} else {
  # Simple GMT reader (name \t desc \t gene1 \t gene2 ...)
  gmt_path <- args$hallmark
  message(sprintf("[Hallmark] Reading GMT from: %s", gmt_path))
  if (!file.exists(gmt_path)) stop("Hallmark GMT file not found: ", gmt_path)
  lines <- readLines(gmt_path, warn = FALSE)
  parsed <- lapply(strsplit(lines, "\t", fixed = TRUE), function(parts) {
    if (length(parts) < 3) return(NULL)
    list(gs_name = parts[[1]], genes = parts[-c(1,2)])
  })
  parsed <- Filter(Negate(is.null), parsed)
  if (length(parsed)) {
    hk_tbl <- rbindlist(lapply(parsed, function(rec) {
      data.table(gs_name = rec$gs_name, symbol = rec$genes)
    }))
  }
}

if (!is.null(hk_tbl) && nrow(hk_tbl)) {
  # Match case-insensitively on symbols
  hk_tbl[, symbol_upper := to_upper(symbol)]
  sym_map <- data.table(Gene.Name = symbols, symbol_upper = to_upper(symbols))
  hk_tbl <- merge(hk_tbl, sym_map, by = "symbol_upper", all.x = TRUE)
  hk_tbl <- hk_tbl[!is.na(Gene.Name)]
  hallmark_per_gene <- hk_tbl[, .(Cancer_Hallmark = paste(sort(unique(gs_name)), collapse = ";")),
                              by = Gene.Name]
  dt <- merge(dt, hallmark_per_gene, by = "Gene.Name", all.x = TRUE, sort = FALSE)
} else {
  dt[, Cancer_Hallmark := NA_character_]
}

# ---- Sort by Borda_Rank (then Borda_Score) if present ----
if ("Borda_Rank" %in% names(dt)) {
  setorderv(dt, c("Borda_Rank","Borda_Score"), order = c(1, -1), na.last = TRUE)
} else if ("Borda_Score" %in% names(dt)) {
  setorder(dt, -Borda_Score, na.last = TRUE)
}

# ---- Write ----
fwrite(dt, args$out_csv)
message(sprintf("Wrote: %s", args$out_csv))
