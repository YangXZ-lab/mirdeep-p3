#!/usr/bin/env Rscript

# =========================================================================
# Sample-sample expression consistency heatmap (miRNA correlation matrix)
# Output: PDF, PNG, SVG
# =========================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ComplexHeatmap)
  library(circlize)      # for colorRamp2
  library(grid)          # for grid.text
})

# -------------------------------------------------------------------------
# Parse command-line arguments
# -------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input expression file (tab-separated, header, row.names=1)", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = "./",
              help = "Output directory for figures [default: %default]", metavar = "DIR")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Sample correlation heatmap using ComplexHeatmap")
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("--input must be specified.", call. = FALSE)
}

if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# -------------------------------------------------------------------------
# Load and preprocess expression data
# -------------------------------------------------------------------------
expr_raw <- read.table(opt$input, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
expr_mat <- as.matrix(expr_raw)

# log2 transformation (add 1 to avoid log(0))
expr_log <- log2(expr_mat + 1)

# Transpose: rows = samples, columns = miRNAs
expr_t <- t(expr_log)

# -------------------------------------------------------------------------
# Calculate sample correlation matrix
# -------------------------------------------------------------------------
cor_mat <- cor(t(expr_t), method = "pearson")

message("Correlation matrix dimensions: ", nrow(cor_mat), " x ", ncol(cor_mat))
message("Samples: ", paste(rownames(cor_mat), collapse = ", "))

# -------------------------------------------------------------------------
# Color mapping: blue (low cor) -> white -> red (high cor)
# Range: min(cor) to 1, typical for biological replicates (often 0.8–1.0)
# You can adjust the lower bound if needed.
# -------------------------------------------------------------------------
cor_min <- min(cor_mat)
col_fun <- colorRamp2(c(cor_min, 1), c("#2166AC", "#B2182B"))
# For a white midpoint, use three colors:
# col_fun <- colorRamp2(c(cor_min, 0.9, 1), c("blue", "white", "red"))

# -------------------------------------------------------------------------
# Define the heatmap
# -------------------------------------------------------------------------
ht <- Heatmap(cor_mat,
              name = "Correlation",
              col = col_fun,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", cor_mat[i, j]), x, y,
                          gp = gpar(fontsize = 10, col = "black"))
              },
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
              clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
              clustering_method_rows = "complete",
              clustering_method_columns = "complete",
              row_names_gp = gpar(fontsize = 12),
              column_names_gp = gpar(fontsize = 12),
              heatmap_legend_param = list(title = "Pearson\nCorrelation",
                                          legend_height = unit(4, "cm")),
              column_title = "Sample-to-sample expression correlation\n(miRNA profiles)",
              column_title_gp = gpar(fontsize = 14, fontface = "bold"),
              border_gp = gpar(col = "grey80", lwd = 0.5))

# -------------------------------------------------------------------------
# Helper function to draw heatmap (to be called inside each graphics device)
# -------------------------------------------------------------------------
draw_heatmap <- function() {
  draw(ht, heatmap_legend_side = "right", padding = unit(c(2, 2, 2, 2), "mm"))
}

# -------------------------------------------------------------------------
# Save in multiple formats
# -------------------------------------------------------------------------
# PDF
pdf(file.path(opt$output, "sample_correlation_heatmap.pdf"), width = 8, height = 7)
draw_heatmap()
dev.off()
message("PDF saved to ", file.path(opt$output, "sample_correlation_heatmap.pdf"))

# PNG
png(file.path(opt$output, "sample_correlation_heatmap.png"), width = 8, height = 7, units = "in", res = 300)
draw_heatmap()
dev.off()
message("PNG saved to ", file.path(opt$output, "sample_correlation_heatmap.png"))

# SVG
svg(file.path(opt$output, "sample_correlation_heatmap.svg"), width = 8, height = 7)
draw_heatmap()
dev.off()
message("SVG saved to ", file.path(opt$output, "sample_correlation_heatmap.svg"))
