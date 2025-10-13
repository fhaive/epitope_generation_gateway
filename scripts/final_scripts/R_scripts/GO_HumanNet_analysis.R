library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(biomaRt)
library(AnnotationDbi)
library(GO.db)
library(org.Hs.eg.db)
library(httr)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
humannet_file_rel <- args[1]
combined_output_rel <- args[2]
not_in_humannet_output_rel <- args[3]
humannet_output_rel <- args[4]
membrane_ensembl_ids_output_rel <- args[5]
# Derive project root
full_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", full_args[grep("^--file=", full_args)])
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "../../.."))

# Debug: Print paths
cat("Project root:", project_root, "\n")
cat("Relative HumanNet file:", humannet_file_rel, "\n")
cat("Relative combined output file:", combined_output_rel, "\n")
cat("Relative not_in_humannet output file:", not_in_humannet_output_rel, "\n")

# Resolve absolute paths
humannet_file <- normalizePath(file.path(project_root, humannet_file_rel))
combined_output <- normalizePath(file.path(project_root, combined_output_rel), mustWork = FALSE)
not_in_humannet_output <- normalizePath(file.path(project_root, not_in_humannet_output_rel), mustWork = FALSE)
humannet_output <- normalizePath(file.path(project_root, humannet_output_rel), mustWork = FALSE)
membrane_ensembl_ids_output <- normalizePath(file.path(project_root, membrane_ensembl_ids_output_rel), mustWork = FALSE)
# Debug: Print absolute paths
cat("Absolute HumanNet file:", humannet_file, "\n")
cat("Absolute combined output file:", combined_output, "\n")
cat("Absolute not_in_humannet output file:", not_in_humannet_output, "\n")

# Check if HumanNet file exists
if (!file.exists(humannet_file)) {
  stop("HumanNet file does not exist: ", humannet_file)
}

# Ensure output directory exists
output_dir <- dirname(combined_output)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Get all descendant GO terms of 'membrane'
roots <- c("GO:0005886", "GO:0005887", "GO:0046658", "GO:0009897", "GO:0009986")
offs <- unique(unlist(AnnotationDbi::mget(roots, GOCCOFFSPRING, ifnotfound = NA)))
go_membrane <- unique(c(roots, offs))

# Download and read GAF file
url <- "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz"
gaf <- tempfile(fileext = ".gaf.gz")
download.file(url, gaf, quiet = TRUE)

# Read GAF file (17 columns, skip comments)
cols <- paste0("V", 1:17)
ann <- read_tsv(gaf, comment = "!", col_names = cols, show_col_types = FALSE)

# Filter for membrane GO IDs and trusted evidence codes
keep_ec <- c("EXP", "IDA", "IMP", "IPI", "IGI", "IEP", "TAS", "IC")
hits <- ann[ann$V5 %in% go_membrane & ann$V7 %in% keep_ec & ann$V4 != "NOT", ]

# Collapse to gene list
genes <- unique(hits$V3)

# Cache file for biomaRt results
cache_file <- file.path(output_dir, "biomart_cache.rds")

# Map gene symbols to Ensembl IDs with error handling, caching, and mirror retries
if (file.exists(cache_file)) {
  cat("Loading cached biomaRt results from", cache_file, "\n")
  mapping <- readRDS(cache_file)
} else {
  # List of hosts to try
  hosts <- c("https://dec2024.archive.ensembl.org", "https://uswest.ensembl.org", 
             "https://useast.ensembl.org", "https://asia.ensembl.org")
  mapping <- NULL
  for (host in hosts) {
    cat("Attempting connection to", host, "\n")
    tryCatch({
      # Set SSL verification off
      set_config(config(ssl_verifypeer = FALSE))
      # Try to list marts
      marts <- listMarts(host = host)
      mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = host)
      mapping <- getBM(
        attributes = c("hgnc_symbol", "ensembl_gene_id"),
        filters = "hgnc_symbol",
        values = genes,
        mart = mart
      )
      # Save to cache
      saveRDS(mapping, cache_file)
      cat("Saved biomaRt results to", cache_file, "\n")
      break  # Exit loop on success
    }, error = function(e) {
      cat("Failed to list marts at", host, ": ", e$message, "\n")
      # Fallback: Direct dataset access
      tryCatch({
        mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = host)
        mapping <<- getBM(
          attributes = c("hgnc_symbol", "ensembl_gene_id"),
          filters = "hgnc_symbol",
          values = genes,
          mart = mart
        )
        saveRDS(mapping, cache_file)
        cat("Saved biomaRt results to", cache_file, "\n")
        return(NULL)  # Exit loop
      }, error = function(e2) {
        cat("Direct access failed at", host, ": ", e2$message, "\n")
        return(NULL)  # Continue to next host
      })
    })
  }
  if (is.null(mapping)) {
    stop("Failed to connect to any Ensembl host or retrieve data with biomaRt")
  }
}
membrane_ensembl_ids <- unique(mapping$ensembl_gene_id)

# Load HumanNet-XN
humannet_raw <- read_tsv(
  humannet_file,
  col_names = c("EntrezID1", "EntrezID2", "LLS"),
  comment = "#",
  col_types = "ccd"
)

# Map Entrez to Ensembl
entrez_vec <- unique(c(humannet_raw$EntrezID1, humannet_raw$EntrezID2))
map_df <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = entrez_vec,
  keytype = "ENTREZID",
  columns = c("ENSEMBL")
)

# Debug: Print column names and first few rows of map_df
cat("Columns in map_df:", colnames(map_df), "\n")
cat("First few rows of map_df:\n")
print(head(map_df))

# Check for ENSEMBL column and rename accordingly
if (!"ENSEMBL" %in% colnames(map_df)) {
  stop("ENSEMBL column not found in map_df. Available columns: ", paste(colnames(map_df), collapse = ", "))
}
map_df <- map_df %>%
  filter(!is.na(ENSEMBL)) %>%
  distinct(ENTREZID, ENSEMBL)

# Attach Ensembl IDs to HumanNet edges
humannet <- humannet_raw %>%
  left_join(map_df, by = c("EntrezID1" = "ENTREZID"), relationship = "many-to-many") %>%
  dplyr::rename(ensembl1 = ENSEMBL) %>%
  left_join(map_df, by = c("EntrezID2" = "ENTREZID"), relationship = "many-to-many") %>%
  dplyr::rename(ensembl2 = ENSEMBL) %>%
  filter(!is.na(ensembl1) & !is.na(ensembl2)) %>%
  dplyr::select(ensembl1, ensembl2)

# Debug: Print first few rows of humannet after joins
cat("First few rows of humannet after joins:\n")
print(head(humannet))

# Get HumanNet Ensembl IDs
humannet_ensembl_ids <- unique(c(humannet$ensembl1, humannet$ensembl2))

# Combine and compare gene lists
combined_ensembl_ids <- unique(c(humannet_ensembl_ids, membrane_ensembl_ids))
not_in_humannet <- setdiff(membrane_ensembl_ids, humannet_ensembl_ids)
cat("First few rows of membrane ensembl ids:")
print(head(membrane_ensembl_ids))

# Save outputs
write_tsv(data.frame(ensembl_id = combined_ensembl_ids), combined_output, col_names = FALSE)
write_tsv(data.frame(ensembl_id = not_in_humannet), not_in_humannet_output, col_names = FALSE)
write_rds(humannet,humannet_output,compress="none")
write_rds(membrane_ensembl_ids,membrane_ensembl_ids_output,compress="none")

