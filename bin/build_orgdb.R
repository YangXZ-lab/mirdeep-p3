#!/usr/bin/env Rscript

# =============================================================================
# Build OrgDb package from eggNOG-mapper annotations and KEGG JSON
# =============================================================================
# Usage:
#   Rscript build_orgdb.R -i emapper.annotations --kojson ko00001.json -o /output/dir
# =============================================================================

# --- 0. Load required packages -----------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(stringr)
  library(jsonlite)
  library(AnnotationForge)
})

options(stringsAsFactors = FALSE)

# --- 1. Parse command-line arguments -----------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to eggNOG-mapper annotations file (tab-separated)",
              metavar = "FILE"),
  make_option(c("--kojson"), type = "character", default = NULL,
              help = "Path to KEGG ko00001.json file",
              metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = ".",
              help = "Output directory for OrgDb package and kegg_info.RData [default: %default]",
              metavar = "DIR")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Build an OrgDb package from eggNOG-mapper and KEGG data.")
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input) || is.null(opt$kojson)) {
  print_help(opt_parser)
  stop("Both --input and --kojson must be specified.", call. = FALSE)
}

input_file  <- opt$input
ko_json     <- opt$kojson
output_dir  <- opt$output

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

message("Input annotations: ", input_file)
message("KEGG JSON: ", ko_json)
message("Output directory: ", output_dir)

# --- 2. Read and prepare emapper annotations ---------------------------------
emapper <- read.table(input_file, header = TRUE, sep = "\t", quote = "",
                      stringsAsFactors = FALSE, check.names = FALSE)
emapper[emapper == ""] <- NA

# Gene information
gene_info <- emapper %>%
  dplyr::select(GID = query, GENENAME = Preferred_name) %>%
  na.omit() %>%
  dplyr::distinct()

# GO terms
gos <- emapper %>%
  dplyr::select(query, GOs) %>%
  na.omit()

# Helper to split comma-separated values into long format
split_to_long <- function(x, col_name) {
  ids <- str_split(x[2], ",", simplify = FALSE)[[1]]
  df <- data.frame(GID = rep(x[1], length(ids)),
                   term = ids,
                   stringsAsFactors = FALSE)
  colnames(df)[2] <- col_name
  return(df)
}

# Create gene2go table
gos_matrix <- as.matrix(gos)
gene2go_list <- apply(gos_matrix, 1, function(r) split_to_long(r, "GO"))
gene2go <- do.call(rbind, gene2go_list)
gene2go <- gene2go %>%
  filter(GO != "-") %>%
  mutate(EVIDENCE = "IEA") %>%
  dplyr::select(GID, GO, EVIDENCE) %>%
  dplyr::distinct()

# Create gene2ko table
gene2ko <- emapper %>%
  dplyr::select(GID = query, Ko = KEGG_ko) %>%
  na.omit() %>%
  filter(Ko != "-")

ko_matrix <- as.matrix(gene2ko)
gene2ko_list <- apply(ko_matrix, 1, function(r) split_to_long(r, "Ko"))
gene2ko <- do.call(rbind, gene2ko_list) %>%
  dplyr::select(GID, Ko) %>%
  mutate(Ko = gsub("ko:", "", Ko)) %>%
  dplyr::distinct()

# --- 3. KEGG pathway parsing from JSON ---------------------------------------
kegg_rdata <- file.path(output_dir, "kegg_info.RData")

# If the RData file does not exist, generate it from the JSON.
if (!file.exists(kegg_rdata)) {
  message("Parsing KEGG JSON and saving to ", kegg_rdata)
  
  pathway2name <- tibble(Pathway = character(), Name = character())
  ko2pathway   <- tibble(Ko = character(), Pathway = character())
  
  kegg <- fromJSON(ko_json)
  
  for (a in seq_along(kegg[["children"]][["children"]])) {
    A <- kegg[["children"]][["name"]][[a]]
    for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
      B <- kegg[["children"]][["children"]][[a]][["name"]][[b]]
      for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
        pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
        pathway_id   <- str_match(pathway_info, "ko[0-9]{5}")[1]
        pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>%
          str_replace("[0-9]{5} ", "")
        pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
        
        kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
        kos      <- str_match(kos_info, "K[0-9]*")[,1]
        ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
      }
    }
  }
  
  save(pathway2name, ko2pathway, file = kegg_rdata)
} else {
  message("Loading existing kegg_info.RData from ", kegg_rdata)
}

load(kegg_rdata)

# Map KO to pathway
gene2pathway <- gene2ko %>%
  left_join(ko2pathway, by = "Ko") %>%
  dplyr::select(GID, Pathway) %>%
  na.omit() %>%
  dplyr::distinct()

# --- 4. Build OrgDb package --------------------------------------------------
makeOrgPackage(
  gene_info   = gene_info,
  go          = gene2go,
  ko          = gene2ko,
  pathway     = gene2pathway,
  version     = "1.46.0",
  maintainer  = "mirdeep <mirdeep@ibcas.ac.cn>",
  author      = "mirdeep <mirdeep@ibcas.ac.cn>",
  outputDir   = output_dir,
  tax_id      = 999999,
  genus       = "My",
  species     = "org",
  goTable     = "go"
)

message("OrgDb package successfully created in ", output_dir)

# --- 5. Output pathway2name and pathway2gene tables --------------------------
# Clean pathway names (remove BR:ko... suffix) and drop any rows with NA
pathway2name$Name <- gsub(" \\[BR:ko[0-9]{5}\\]", "", pathway2name$Name)
pathway2name <- na.omit(pathway2name)

# Extract pathway-gene mapping
pathway2gene <- gene2pathway[, c("Pathway", "GID")]

# Write to output directory
write.table(pathway2name,
            file = file.path(output_dir, "pathway2name"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pathway2gene,
            file = file.path(output_dir, "pathway2gene"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("pathway2name and pathway2gene written to ", output_dir)
