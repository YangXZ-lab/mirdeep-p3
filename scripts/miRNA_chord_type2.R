#!/usr/bin/env Rscript

# =========================================================================
# miRNA ↔ Pathway Chord Diagram (Global GO + KEGG enrichment)
# =========================================================================
# Description:
#   1. Read miRNA–target gene pairs.
#   2. Run GO and KEGG enrichment on the union of all target genes.
#   3. For each miRNA and each significant term, count how many target genes
#      of that miRNA participate in the term.
#   4. Draw chord diagrams where left semicircle = miRNAs, right = pathways,
#      and chord thickness = number of overlapping genes.
#
# Usage:
#   Rscript miRNA_chord_global.R -i <input> --orgdb <orgdb_dir> -f <annot_dir> -o <out_dir>
#                                [--go_pvalue 0.05] [--go_qvalue 1]
#                                [--kegg_pvalue 0.05] [--kegg_qvalue 1]
# =========================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(circlize)
  library(openxlsx)
  library(AnnotationDbi)
})

# ---- Command‑line options (updated) ----------------------------------------
option_list <- list(
  make_option(c("-i", "--input"),  type = "character", default = NULL,
              help = "Input file (miRNA <tab> gene, no header) [required]"),
  make_option("--orgdb", type = "character", default = NULL,
              help = "Path to org.Morg.eg.db source directory [required]"),
  make_option(c("-f", "--file"),   type = "character", default = NULL,
              help = "Directory containing 'pathway2gene' and 'pathway2name' [required]"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output directory [required]"),
  make_option("--go_pvalue",   type = "double", default = 0.05,
              help = "p‑value cutoff for GO [default: 0.05]"),
  make_option("--go_qvalue",   type = "double", default = 1,
              help = "q‑value cutoff for GO [default: 1 (off)]"),
  make_option("--kegg_pvalue", type = "double", default = 0.05,
              help = "p‑value cutoff for KEGG [default: 0.05]"),
  make_option("--kegg_qvalue", type = "double", default = 1,
              help = "q‑value cutoff for KEGG [default: 1 (off)]"),
  make_option("--go_topn",   type = "integer", default = 10,
              help = "Top N GO terms for chord plot [default: 10]"),
  make_option("--kegg_topn", type = "integer", default = 10,
              help = "Top N KEGG terms for chord plot (trumps p/q when fewer are significant) [default: 10]")
)

opt <- parse_args(OptionParser(option_list = option_list,
                               description = "miRNA–pathway chord diagram (global enrichment)"))

if (is.null(opt$input))  stop("Missing -i/--input")
if (is.null(opt$orgdb))  stop("Missing --orgdb")
if (is.null(opt$file))   stop("Missing -f/--file")
if (is.null(opt$output)) stop("Missing -o/--output")

input_file     <- opt$input
orgdb_dir      <- opt$orgdb
annotation_dir <- opt$file
output_dir     <- opt$output
go_pval        <- opt$go_pvalue
go_qval        <- opt$go_qvalue
kegg_pval      <- opt$kegg_pvalue
kegg_qval      <- opt$kegg_qvalue
go_topn        <- opt$go_topn
kegg_topn      <- opt$kegg_topn

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Read miRNA ↔ gene pairs --------------------------------------------
message("Reading input: ", input_file)
mir_gene_raw <- read.table(input_file, header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE, comment.char = "")
if (ncol(mir_gene_raw) < 2) stop("Input must contain at least two columns (miRNA, gene).")
colnames(mir_gene_raw)[1:2] <- c("miRNA", "gene")
mir_gene <- unique(mir_gene_raw[, c("miRNA", "gene")])

mirna_list  <- sort(unique(mir_gene$miRNA))
mir2genes   <- split(mir_gene$gene, mir_gene$miRNA)
all_genes   <- unique(mir_gene$gene)
message(sprintf("miRNAs: %d  |  unique target genes: %d", length(mirna_list), length(all_genes)))

# ---- 2. Install and load OrgDb ---------------------------------------------
message("Installing OrgDb: ", orgdb_dir)
install.packages(orgdb_dir, repos = NULL, type = "source")
library(org.Morg.eg.db)

# ---- 3. KEGG annotation files ----------------------------------------------
pathway2gene_file <- file.path(annotation_dir, "pathway2gene")
pathway2name_file <- file.path(annotation_dir, "pathway2name")
if (!file.exists(pathway2gene_file)) stop("Missing: ", pathway2gene_file)
if (!file.exists(pathway2name_file)) stop("Missing: ", pathway2name_file)

pathway2gene <- read.table(pathway2gene_file, header = TRUE, sep = "\t")
pathway2name <- read.table(pathway2name_file, header = TRUE, sep = "\t")
colnames(pathway2gene)[1:2] <- c("pathway", "gene")

# ---- 4. Global GO enrichment -----------------------------------------------
message("Running GO enrichment on all target genes ...")
ego <- enrichGO(gene          = all_genes,
                OrgDb         = org.Morg.eg.db,
                keyType       = "GID",
                ont           = "ALL",
                pvalueCutoff  = go_pval,
                qvalueCutoff  = go_qval,
                pAdjustMethod = "BH",
                minGSSize     = 1,
                maxGSSize     = 5000)

go_result  <- ego@result
go_result  <- go_result[order(go_result$pvalue), ]   # keep all, sort by p
# Extract gene lists for all terms
go_genes_all <- lapply(go_result$geneID, function(x) unlist(strsplit(x, "/")))
names(go_genes_all) <- go_result$Description

# ---- 5. Global KEGG enrichment ---------------------------------------------
message("Running KEGG enrichment on all target genes ...")
ekp <- enricher(gene          = all_genes,
                TERM2GENE     = pathway2gene,
                TERM2NAME     = pathway2name,
                pvalueCutoff  = kegg_pval,
                qvalueCutoff  = kegg_qval,
                pAdjustMethod = "BH",
                minGSSize     = 1)

kegg_result <- ekp@result
kegg_result <- kegg_result[order(kegg_result$pvalue), ]
kegg_genes_all <- lapply(kegg_result$geneID, function(x) unlist(strsplit(x, "/")))
names(kegg_genes_all) <- kegg_result$Description

# ---- 6. Build miRNA × pathway matrix (overlap counts) ----------------------
build_mirna_term_matrix <- function(mir2genes, term2genes, term_names) {
  m <- matrix(0, nrow = length(mir2genes), ncol = length(term_names),
              dimnames = list(names(mir2genes), term_names))
  for (i in seq_along(mir2genes)) {
    g <- mir2genes[[i]]
    for (j in seq_along(term_names)) {
      tj <- term2genes[[ term_names[j] ]]
      if (!is.null(tj)) m[i, j] <- length(intersect(g, tj))
    }
  }
  m
}

# ---- 7. Select top terms for chord plots -----------------------------------
# GO: top go_topn by number of connected miRNAs (or all if fewer)
go_sig_desc <- go_result$Description
go_term_genes <- go_genes_all[go_sig_desc]
go_full_mat <- build_mirna_term_matrix(mir2genes, go_term_genes, go_sig_desc)
# rank by colSums(mat > 0)
go_scores <- colSums(go_full_mat > 0)
go_selected <- head(names(sort(go_scores, decreasing = TRUE)), go_topn)
go_mat <- go_full_mat[, go_selected, drop = FALSE]

# KEGG: top kegg_topn terms by p‑value (even if not passing strict thresholds)
kegg_desc <- kegg_result$Description
kegg_term_genes <- kegg_genes_all[kegg_desc]
kegg_full_mat <- build_mirna_term_matrix(mir2genes, kegg_term_genes, kegg_desc)
kegg_selected <- head(kegg_desc, kegg_topn)
kegg_mat <- kegg_full_mat[, kegg_selected, drop = FALSE]

# ---- 8. Chord diagram drawing (legend at bottom) --------------------------
draw_chord <- function(mat, type, outdir) {
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    message("No links for ", type, " — skipping.")
    return(invisible(NULL))
  }
  if (sum(mat) == 0) {
    message("All zero connections for ", type, " — skipping.")
    return(invisible(NULL))
  }
  
  # Filter rows with zero connections
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    message("No connections after row filter — skipping.")
    return(invisible(NULL))
  }
  
  mirnas <- rownames(mat)
  terms  <- colnames(mat)
  
  # Colours
  mir_cols   <- setNames(rainbow(length(mirnas), s = 0.8, v = 0.8), mirnas)
  term_cols  <- setNames(rand_color(length(terms), luminosity = "bright"), terms)
  grid_col   <- c(mir_cols, term_cols)
  
  # Gap between left (miRNA) and right (pathway) semicircles
  gap <- rep(2, length(mirnas) + length(terms))
  gap[length(mirnas)] <- 15
  
  # Chord colour matrix
  link_colors <- matrix(NA, nrow(mat), ncol(mat), dimnames = dimnames(mat))
  for (i in seq_along(mirnas)) {
    for (j in seq_along(terms)) {
      if (mat[i, j] > 0) link_colors[i, j] <- mir_cols[mirnas[i]]
    }
  }
  
  plot_chord <- function(ext) {
    if (ext == "pdf") {
      pdf(file.path(outdir, paste0(type, "_chord.pdf")), width = 20, height = 16)
    } else {
      png(file.path(outdir, paste0(type, "_chord.png")), width = 20, height = 16,
          units = "in", res = 300)
    }
    
    # Top: chord diagram (3/4 height), Bottom: legend (1/4 height)
    layout(matrix(c(1, 2), nrow = 2), heights = c(3, 1))
    
    # ---- Chord diagram ----
    circos.clear()
    circos.par(gap.after = gap, start.degree = 90)
    suppressWarnings(
      chordDiagram(mat,
                   order                = c(mirnas, terms),
                   grid.col             = grid_col,
                   col                  = link_colors,
                   directional          = 1,
                   direction.type       = c("diffHeight", "arrows"),
                   link.arr.type        = "big.arrow",
                   diffHeight           = mm_h(2),
                   annotationTrack      = c("name", "grid"),
                   annotationTrackHeight = c(0.03, 0.05),
                   scale                = TRUE,
                   reduce               = 0)
    )
    title(main = paste(type, "miRNA-Pathway Regulatory Network"),
          cex.main = 1.2, line = -2)
    
    # ---- Legend for pathways (bottom) ----
    par(mar = c(1, 1, 1, 1))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    legend("center",
           legend = terms,
           fill   = term_cols[terms],
           border = NA,
           bty    = "n",
           cex    = 0.7,
           ncol   = 4,
           title  = paste(type, "Pathways"),
           title.adj = 0,
           x.intersp = 0.5,
           text.width = 0.22)
    dev.off()
  }
  
  plot_chord("pdf")
  plot_chord("png")
  
  message(type, " chord diagram saved to ", outdir)
}

# ---- 9. Draw both diagrams -------------------------------------------------
draw_chord(go_mat,   "GO",   output_dir)
draw_chord(kegg_mat, "KEGG", output_dir)

# ---- 10. Export matrices and enrichment tables -----------------------------
if (!is.null(go_mat)) {
  write.table(as.data.frame(go_full_mat), file.path(output_dir, "GO_count_matrix.tsv"),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  write.xlsx(list(GO_Enrichment = go_result),
             file.path(output_dir, "GO_enrichment.xlsx"), rowNames = FALSE)
}
if (!is.null(kegg_mat)) {
  write.table(as.data.frame(kegg_full_mat), file.path(output_dir, "KEGG_count_matrix.tsv"),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  write.xlsx(list(KEGG_Enrichment = kegg_result),
             file.path(output_dir, "KEGG_enrichment.xlsx"), rowNames = FALSE)
}

message("All done. Output directory: ", output_dir)
