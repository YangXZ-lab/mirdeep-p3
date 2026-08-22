#!/usr/bin/env Rscript

# =========================================================================
# Heatmap of differentially expressed miRNAs using ComplexHeatmap
# --input : expression matrix (tab-separated, header, row.names=1)
# --deg   : DEG result file (optional) to filter DE miRNAs
# --output: output directory (PDF, PNG, SVG)
# Rows: DE miRNAs (clustered), Columns: samples (original order)
# =========================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ComplexHeatmap)
  library(circlize)    # for colorRamp2
})

# -------------------------------------------------------------------------
# Parse command-line arguments
# -------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Expression matrix (tab-separated, header, row.names=1)", metavar = "FILE"),
  make_option(c("--deg"), type = "character", default = NULL,
              help = "DEG result table (optional). Filter miRNAs where 'change' != 'NOT'",
              metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = "./",
              help = "Output directory [default: %default]", metavar = "DIR")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Heatmap of miRNA expression (DE miRNAs if --deg provided)")
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("--input is required.", call. = FALSE)
}

if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# -------------------------------------------------------------------------
# Load expression matrix
# -------------------------------------------------------------------------
expr_raw <- read.table(opt$input, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
expr_mat <- as.matrix(expr_raw)
message("Loaded expression matrix: ", nrow(expr_mat), " miRNAs × ", ncol(expr_mat), " samples")

# -------------------------------------------------------------------------
# Optional DEG filtering
# -------------------------------------------------------------------------
if (!is.null(opt$deg)) {
  deg <- read.table(opt$deg, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  if (!"change" %in% colnames(deg)) {
    stop("DEG file must contain a 'change' column.", call. = FALSE)
  }
  de_mirnas <- deg$miRNA[deg$change != "NOT"]
  if (length(de_mirnas) == 0) {
    stop("No DE miRNAs found (change != 'NOT') in the DEG file.", call. = FALSE)
  }
  # Keep only miRNAs present in the expression matrix
  common <- intersect(de_mirnas, rownames(expr_mat))
  if (length(common) == 0) {
    stop("None of the DE miRNAs were found in the expression matrix row names.", call. = FALSE)
  }
  expr_mat <- expr_mat[common, , drop = FALSE]
  message("Filtered to ", length(common), " DE miRNAs.")
} else {
  message("No DEG file provided, using all miRNAs.")
}

# -------------------------------------------------------------------------
# Row-wise Z-score scaling (for visual comparison)
# -------------------------------------------------------------------------
expr_scaled <- t(scale(t(expr_mat)))
# Remove rows that become all NA (e.g., constant expression)
na_rows <- apply(expr_scaled, 1, function(x) all(is.na(x)))
if (any(na_rows)) {
  message("Removing ", sum(na_rows), " miRNA(s) with no variation.")
  expr_scaled <- expr_scaled[!na_rows, , drop = FALSE]
}

# -------------------------------------------------------------------------
# Color mapping: soft blue-white-orange/red
# -------------------------------------------------------------------------
max_abs <- max(abs(expr_scaled), na.rm = TRUE)
# Soft, less saturated gradient
col_fun <- colorRamp2(
  breaks = seq(-max_abs, max_abs, length = 5),
  colors = c("#4393c3", "#d1e5f0", "white", "#fddbc7", "#d6604d")
)

# -------------------------------------------------------------------------
# Row clustering, columns keep original order
# -------------------------------------------------------------------------
if (nrow(expr_scaled) > 1) {
  row_dist <- dist(expr_scaled)
  row_hclust <- hclust(row_dist, method = "complete")
  cluster_rows <- row_hclust
  show_row_dend <- TRUE
} else {
  cluster_rows <- FALSE
  show_row_dend <- FALSE
}

col_order <- colnames(expr_scaled)

# -------------------------------------------------------------------------
# Heatmap definition
# -------------------------------------------------------------------------
ht <- Heatmap(
  expr_scaled,
  name = "Z-score",
  col = col_fun,
  cluster_rows = cluster_rows,
  row_dend_side = "left",
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = max(5, 8 - 0.2 * nrow(expr_scaled))),
  row_title = "miRNAs",
  cluster_columns = FALSE,
  column_order = col_order,
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 9),
  column_title = "Samples",
  border = TRUE,
  use_raster = TRUE,
  raster_quality = 2,
  heatmap_legend_param = list(
    title = "Row Z-score",
    legend_direction = "horizontal",
    title_position = "topcenter"
  )
)

# -------------------------------------------------------------------------
# Export to PDF, PNG, SVG
# -------------------------------------------------------------------------
save_heatmap <- function(filename, fmt) {
  width_in  <- max(8, ncol(expr_scaled) * 0.5)
  height_in <- max(6, nrow(expr_scaled) * 0.25)
  if (fmt == "pdf") {
    pdf(filename, width = width_in, height = height_in)
  } else if (fmt == "png") {
    png(filename, width = width_in * 100, height = height_in * 100, res = 150)
  } else if (fmt == "svg") {
    svg(filename, width = width_in, height = height_in)
  }
  draw(ht, heatmap_legend_side = "bottom")
  dev.off()
}

formats <- c("pdf", "png", "svg")
for (fmt in formats) {
  outfile <- file.path(opt$output, paste0("heatmap.", fmt))
  save_heatmap(outfile, fmt)
  message("Heatmap saved: ", outfile)
}
