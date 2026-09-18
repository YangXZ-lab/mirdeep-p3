#!/usr/bin/env Rscript

# =========================================================================
# KEGG & GO Enrichment Analysis (using clusterProfiler)
# =========================================================================
# Description:
#   1. Install a local OrgDb package from source.
#   2. Read pre‑prepared pathway annotation files (pathway2gene, pathway2name).
#   3. Read a gene list.
#   4. Run KEGG enrichment (enricher) and GO enrichment (enrichGO).
#   5. Output:
#       - KEGG dotplot (kegg.pdf/png/svg)
#       - GO dotplot (go.pdf/png/svg)
#       - GO barplot by category (GO_bar.pdf/png/svg)
#       - GO DAG plot (GO_dag.pdf/png/svg)
#       - GO results Excel file (ego_list.xlsx)
#
# Exit codes:
#   0 = enrichment results were produced
#   3 = ran successfully but no significant KEGG/GO term was found (empty result)
#   1 = error
#
# Usage:
#   Rscript enrich_analysis.R -i <orgdb_dir> -f <annotation_dir> -o <output_dir> -g <gene_file>
#                             [--kegg_pvalue <p>] [--kegg_qvalue <q>]
#                             [--go_pvalue <p>]   [--go_qvalue <q>]
#                             [--goterm <n>]
#
# Arguments:
#   -i, --input     Path to the org.Morg.eg.db source directory
#   -f, --file      Directory containing "pathway2gene" and "pathway2name"
#   -o, --output    Output directory (created if absent)
#   -g, --gene      Gene list file (one column, no header)
#   --kegg_pvalue   p‑value cutoff for KEGG (default: 1 = keep all)
#   --kegg_qvalue   q‑value cutoff for KEGG (default: 1 = keep all)
#   --go_pvalue     p‑value cutoff for GO (default: 0.05)
#   --go_qvalue     q‑value cutoff for GO (default: 1 = keep all)
#   --goterm        Number of top terms per GO category in barplot (default: all)
# =========================================================================

# ---- Load necessary packages (with installation checks) ----
required_pkgs <- c("optparse", "clusterProfiler", "ggplot2", "openxlsx", "enrichplot")
installed <- required_pkgs %in% installed.packages()[,"Package"]
if (any(!installed)) {
  stop("The following packages are missing: ",
       paste(required_pkgs[!installed], collapse = ", "),
       "\nPlease install them before running this script.")
}

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(ggplot2)
  library(openxlsx)
  library(enrichplot)
  library(GOSemSim)
})

# ---- Define command‑line options ----
option_list <- list(
  make_option(c("-i", "--input"),  type = "character", default = NULL,
              help = "Path to the org.Morg.eg.db source directory [required]",
              metavar = "character"),
  make_option(c("-f", "--file"),   type = "character", default = NULL,
              help = "Directory containing 'pathway2gene' and 'pathway2name' [required]",
              metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output directory [required]",
              metavar = "character"),
  make_option(c("-g", "--gene"),   type = "character", default = NULL,
              help = "Gene list file (one column, no header) [required]",
              metavar = "character"),
  make_option("--kegg_pvalue", type = "double", default = 1,
              help = "p‑value cutoff for KEGG (default: 1 = keep all)",
              metavar = "number"),
  make_option("--kegg_qvalue", type = "double", default = 1,
              help = "q‑value cutoff for KEGG (default: 1 = keep all)",
              metavar = "number"),
  make_option("--go_pvalue", type = "double", default = 0.05,
              help = "p‑value cutoff for GO (default: 0.05)",
              metavar = "number"),
  make_option("--go_qvalue", type = "double", default = 1,
              help = "q‑value cutoff for GO (default: 1 = keep all)",
              metavar = "number"),
  make_option("--goterm", type = "integer", default = NULL,
              help = "Number of top terms per GO category in barplot (default: all)",
              metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "KEGG and GO enrichment analysis with custom OrgDb")
opt <- parse_args(opt_parser)

# ---- Validate required arguments ----
if (is.null(opt$input))  stop("Missing required argument: -i/--input")
if (is.null(opt$file))   stop("Missing required argument: -f/--file")
if (is.null(opt$output)) stop("Missing required argument: -o/--output")
if (is.null(opt$gene))   stop("Missing required argument: -g/--gene")

orgdb_dir      <- opt$input
annotation_dir <- opt$file
output_dir     <- opt$output
gene_file      <- opt$gene
kegg_pval      <- opt$kegg_pvalue
kegg_qval      <- opt$kegg_qvalue
go_pval        <- opt$go_pvalue
go_qval        <- opt$go_qvalue
goterm         <- ifelse(is.null(opt$goterm), Inf, opt$goterm)   # Inf means all

# ---- Create output directory ----
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Install and load the local OrgDb (job-local library) --------------------
# The OrgDb is rebuilt on every run, so it must NOT go into the shared conda R
# library: concurrent jobs would deadlock on 00LOCK-* and the base environment
# would be polluted.  MIRDEEP_R_LIB is exported by the Python driver; when this
# script is run standalone we fall back to <output_dir>/Rlib.
orgdb_pkg <- "org.Morg.eg.db"
job_lib   <- Sys.getenv("MIRDEEP_R_LIB", unset = "")
if (!nzchar(job_lib)) job_lib <- file.path(output_dir, "Rlib")
dir.create(job_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(job_lib, .libPaths()))
message("Job-local R library: ", job_lib)

stale_lock <- file.path(job_lib, paste0("00LOCK-", orgdb_pkg))
if (dir.exists(stale_lock)) {
  message("Removing stale lock: ", stale_lock)
  unlink(stale_lock, recursive = TRUE, force = TRUE)
}

# NOTE: check the JOB library only.  requireNamespace() would also see a stale
# org.Morg.eg.db left in the shared conda library, silently skip the install and
# then analyse with out-of-date annotations.
job_lib_pkgs <- rownames(installed.packages(lib.loc = job_lib))
if (!(orgdb_pkg %in% job_lib_pkgs)) {
  message("Installing OrgDb from: ", orgdb_dir)
  install.packages(orgdb_dir, repos = NULL, type = "source", lib = job_lib)
  job_lib_pkgs <- rownames(installed.packages(lib.loc = job_lib))
  if (!(orgdb_pkg %in% job_lib_pkgs)) {
    stop("Failed to install ", orgdb_pkg, " into ", job_lib,
         " - check the install output above")
  }
} else {
  message("OrgDb already present in ", job_lib, " - skipping installation")
}
suppressPackageStartupMessages(library(org.Morg.eg.db))


# ---- 2. Read pathway annotation files (for KEGG) ----
pathway2gene_file <- file.path(annotation_dir, "pathway2gene")
pathway2name_file <- file.path(annotation_dir, "pathway2name")

if (!file.exists(pathway2gene_file)) stop("File not found: ", pathway2gene_file)
if (!file.exists(pathway2name_file)) stop("File not found: ", pathway2name_file)

pathway2gene <- read.table(pathway2gene_file, header = TRUE, sep = "\t")
pathway2name <- read.table(pathway2name_file, header = TRUE, sep = "\t")

# ---- 3. Read gene list (shared between KEGG and GO) ----
gene <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)
gene_list <- gene$V1

if (length(gene_list) == 0) {
  warning("The gene list is empty. Please check your input file.")
}

# ---- 4. KEGG enrichment analysis ----
message("Running KEGG enrichment...")
ekp <- enricher(gene_list,
                TERM2GENE    = pathway2gene,
                TERM2NAME    = pathway2name,
                pvalueCutoff  = kegg_pval,
                qvalueCutoff  = kegg_qval,
                pAdjustMethod = "BH",
                minGSSize     = 1,
                maxGSSize     = 5000)

# ---- 5. KEGG dotplot ----
if (!is.null(ekp) && nrow(ekp) > 0) {
  kegg_plot <- dotplot(ekp) +
    theme_bw(base_size = 12) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(title = "KEGG Pathway Enrichment") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, "kegg.pdf"), plot = kegg_plot, width = 12, height = 16)
  ggsave(file.path(output_dir, "kegg.png"), plot = kegg_plot, width = 12, height = 16, dpi = 300)
  ggsave(file.path(output_dir, "kegg.svg"), plot = kegg_plot, width = 12, height = 16)
} else {
  message("No significant KEGG terms found, skipping plot.")
}

# ---- 6. GO enrichment analysis ----
message("Running GO enrichment...")
gene_list_GO <- list(test = gene_list)

ego_list <- lapply(gene_list_GO, function(x) {
  enrichGO(gene          = x,
           OrgDb         = org.Morg.eg.db,
           pAdjustMethod = "none",
           keyType       = "GID",
           ont           = "ALL",
           pvalueCutoff  = go_pval,
           qvalueCutoff  = go_qval,
           minGSSize     = 10,
           maxGSSize     = 2500)
})

# enrichGO() returns NULL when none of the genes could be mapped to a GO term.
# That is an empty result, not an error: guard it instead of dereferencing NULL
# with `x@result`, which aborted the whole script with exit code 1.
ego_valid <- ego_list[!vapply(ego_list, is.null, logical(1))]

# ---- 7. Export GO results to Excel ----
excel_path <- file.path(output_dir, "ego_list.xlsx")
if (length(ego_valid) > 0) {
  ego_result_list <- lapply(ego_valid, function(x) x@result)
  write.xlsx(ego_result_list, file = excel_path, rowNames = FALSE)
  message("GO results saved to ", excel_path)
} else {
  message("No GO enrichment result to export (no gene mapped to a GO term) - ",
          "skipping ", excel_path)
}

# ---- 8. GO dotplot (for the first element, if significant) ----
ego_first <- ego_list[[1]]
if (!is.null(ego_first) && nrow(ego_first) > 0) {
  go_dot <- dotplot(ego_first) +
    theme_bw(base_size = 12) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(title = "GO Enrichment") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, "go.pdf"), plot = go_dot, width = 12, height = 16)
  ggsave(file.path(output_dir, "go.png"), plot = go_dot, width = 12, height = 16, dpi = 300)
  ggsave(file.path(output_dir, "go.svg"), plot = go_dot, width = 12, height = 16)
} else {
  message("No significant GO terms found, skipping dotplot.")
}

# ---- 9. GO Barplot (Biological Process, Cellular Component, Molecular Function) ----
if (!is.null(ego_first) && nrow(ego_first) > 0) {
  go_res <- ego_first@result
  
  if (!"ONTOLOGY" %in% colnames(go_res)) {
    stop("ONTOLOGY column not found in GO results. Cannot create barplot.")
  }
  
  # Keep only BP, CC, MF
  go_res <- go_res[go_res$ONTOLOGY %in% c("BP", "CC", "MF"), ]
  if (nrow(go_res) == 0) {
    message("No BP/CC/MF terms found, skipping GO barplot.")
  } else {
    # Limit to top N per category if goterm is specified
    if (!is.infinite(goterm)) {
      go_res <- do.call(rbind, by(go_res, go_res$ONTOLOGY, function(df) {
        df <- df[order(df$Count, decreasing = TRUE), ]
        head(df, goterm)
      }, simplify = FALSE))
    }
    
    # Order Description by Count descending for plotting
    go_res$Description <- factor(go_res$Description,
                                 levels = go_res$Description[order(go_res$Count, decreasing = TRUE)])
    
    category_colors <- c("BP" = "#73A3CF", "CC" = "#AAC88A", "MF" = "#F3B17A")
    
    go_bar <- ggplot(go_res, aes(x = Description, y = Count, fill = ONTOLOGY)) +
      geom_col(width = 0.7) +
      facet_wrap(~ ONTOLOGY, scales = "free_x", nrow = 1) +
      scale_fill_manual(values = category_colors) +
      labs(title = "GO Enrichment BarPlot",
           x = "GO Term",
           y = "Gene Number") +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            strip.background = element_rect(fill = "grey90"),
            legend.position = "none")
    
    ggsave(file.path(output_dir, "GO_bar.pdf"), plot = go_bar, width = 24, height = 16)
    ggsave(file.path(output_dir, "GO_bar.png"), plot = go_bar, width = 24, height = 16, dpi = 300)
    ggsave(file.path(output_dir, "GO_bar.svg"), plot = go_bar, width = 14, height = 16)
  }
} else {
  message("GO enrichment result is empty, skipping barplot.")
}

# ---- 10. GO DAG plot (using enrichplot::goplot) ----
# if (!is.null(ego_first) && nrow(ego_first) > 0) {
#   message("Generating GO DAG plot...")
#   # Limit to at most 20 nodes for clarity; can be adjusted
#   n_terms <- min(20, nrow(ego_first))
#   dag_plot <- goplot(ego_first, showCategory = n_terms) +
#     labs(title = "GO DAG (top significant terms)") +
#     theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# 
#   ggsave(file.path(output_dir, "GO_dag.pdf"), plot = dag_plot, width = 14, height = 14)
#   ggsave(file.path(output_dir, "GO_dag.png"), plot = dag_plot, width = 14, height = 14, dpi = 300)
#   ggsave(file.path(output_dir, "GO_dag.svg"), plot = dag_plot, width = 14, height = 14)
# } else {
#   message("GO enrichment results empty, skipping DAG plot.")
# }

# ---- 11. Report the outcome through the exit code ----------------------------
# Exit codes: 0 = produced results, 3 = ran fine but nothing significant,
#             1 = error (stop() above).
n_kegg_terms <- if (is.null(ekp)) 0L else nrow(ekp)
n_go_terms   <- sum(vapply(ego_list, function(x) if (is.null(x)) 0L else nrow(x),
                           integer(1)))

message("Enrichment finished: ", n_kegg_terms, " KEGG term(s), ",
        n_go_terms, " GO term(s).")

if (n_kegg_terms == 0L && n_go_terms == 0L) {
  message("No significant KEGG or GO terms found for this gene list - ",
          "this is an empty result, not an error (exit code 3).")
  quit(save = "no", status = 3)
}

message("All enrichment analyses completed. Results are in: ", output_dir)
