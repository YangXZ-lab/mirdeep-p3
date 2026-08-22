#!/usr/bin/env Rscript

# =========================================================================
# PCA analysis of miRNA expression matrix
# Output: scatter plot (PDF, PNG, SVG) and scree plot (PDF, PNG, SVG)
# =========================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(ggrepel)
})

# -------------------------------------------------------------------------
# Parse command-line arguments
# -------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input expression file (tab-separated, header, row.names=1)", metavar = "FILE"),
  make_option(c("-r", "--replicate"), type = "integer", default = NULL,
              help = "Number of replicates per case (2 cases total)", metavar = "INT"),
  make_option(c("-o", "--output"), type = "character", default = "./",
              help = "Output directory for figures [default: %default]", metavar = "DIR")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "PCA analysis for miRNA expression data")
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$replicate)) {
  print_help(opt_parser)
  stop("Both --input and --replicate must be specified.", call. = FALSE)
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

# Transpose: rows = samples, columns = genes
expr_t <- t(expr_log)

# -------------------------------------------------------------------------
# Filter out genes (columns) with zero variance
# -------------------------------------------------------------------------
gene_sd <- apply(expr_t, 2, sd)
zero_var <- gene_sd == 0
if (any(zero_var)) {
  message(sprintf("Removing %d genes with zero variance (constant expression across samples).",
                  sum(zero_var)))
  expr_t <- expr_t[, !zero_var, drop = FALSE]
}

# -------------------------------------------------------------------------
# Check sample count
# -------------------------------------------------------------------------
expected_samples <- 2 * opt$replicate
if (nrow(expr_t) != expected_samples) {
  stop(sprintf(
    "Expected %d samples (2 cases × %d replicates), but found %d samples: %s.",
    expected_samples, opt$replicate,
    nrow(expr_t), paste(rownames(expr_t), collapse = ", ")
  ), call. = FALSE)
}

# -------------------------------------------------------------------------
# Build sample metadata
# -------------------------------------------------------------------------
n_rep <- opt$replicate
sample_info <- data.frame(
  sample = rownames(expr_t),
  case   = rep(c("case1", "case2"), each = n_rep),
  replicate = rep(seq_len(n_rep), 2),
  stringsAsFactors = FALSE
)

# -------------------------------------------------------------------------
# Perform PCA
# -------------------------------------------------------------------------
pca_res <- prcomp(expr_t, center = TRUE, scale. = TRUE)
var_exp <- summary(pca_res)$importance[2, ] * 100

pca_scores <- as.data.frame(pca_res$x)
pca_scores$sample <- rownames(pca_scores)

plot_df <- merge(sample_info, pca_scores, by = "sample")

# Determine if confidence ellipses can be drawn (need >= 4 points per group)
can_draw_ellipse <- all(table(plot_df$case) >= 4)
if (!can_draw_ellipse) {
  message("Note: Each case has fewer than 4 points, confidence ellipses will be omitted.")
}

# -------------------------------------------------------------------------
# PCA scatter plot (PC1 vs PC2)
# -------------------------------------------------------------------------
p <- ggplot(plot_df, aes(x = PC1, y = PC2, color = case, shape = factor(replicate))) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(aes(label = sample), size = 3.2, show.legend = FALSE) +
  labs(
    title = "PCA of miRNA Expression Profiles",
    x = paste0("PC1 (", round(var_exp[1], 1), "% variance)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "% variance)"),
    color = "Case",
    shape = "Replicate"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Add ellipses only if group size >= 4
if (can_draw_ellipse) {
  p <- p + stat_ellipse(aes(group = case), level = 0.95,
                        linetype = "dashed", linewidth = 0.8)
}

# -------------------------------------------------------------------------
# Save PCA plot in multiple formats
# -------------------------------------------------------------------------
formats <- c("pdf", "png", "svg")
for (fmt in formats) {
  outfile <- file.path(opt$output, paste0("PCA_scatter.", fmt))
  ggsave(outfile, plot = p, width = 8, height = 7, device = fmt)
  message("PCA scatter plot saved to ", outfile)
}

# -------------------------------------------------------------------------
# Scree plot
# -------------------------------------------------------------------------
scree_df <- data.frame(
  PC = factor(names(var_exp), levels = names(var_exp)),
  Variance = var_exp
)

p_scree <- ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
  geom_text(aes(label = paste0(round(Variance, 1), "%")), vjust = -0.5, size = 3) +
  labs(title = "Scree Plot", x = "Principal Component", y = "Variance Explained (%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Save scree plot in multiple formats
for (fmt in formats) {
  outfile <- file.path(opt$output, paste0("scree_plot.", fmt))
  ggsave(outfile, plot = p_scree, width = 7, height = 5, device = fmt)
  message("Scree plot saved to ", outfile)
}
