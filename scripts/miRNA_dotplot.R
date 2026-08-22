#!/usr/bin/env Rscript

# =========================================================================
# Individual miRNA expression dot plot (one figure per miRNA)
# Supports filtering by --miRNA or -f
# =========================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(reshape2)
})

# -------------------------------------------------------------------------
# Parse command-line arguments
# -------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input expression file (tab-separated, header, row.names=1)", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = "./",
              help = "Output directory for figures [default: %default]", metavar = "DIR"),
  make_option(c("--miRNA"), type = "character", default = NULL,
              help = "Comma-separated list of miRNA names to plot (e.g., 'miR-21-5p,miR-155-5p')", metavar = "STR"),
  make_option(c("-f", "--file"), type = "character", default = NULL,
              help = "File containing miRNA names, one per line", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Individual miRNA expression dot plot (separate figures)")
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("--input must be specified.", call. = FALSE)
}

if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# -------------------------------------------------------------------------
# Resolve target miRNA list
# -------------------------------------------------------------------------
target_mirnas <- NULL

# From --miRNA (comma-separated string)
if (!is.null(opt$miRNA)) {
  mir_list <- trimws(unlist(strsplit(opt$miRNA, ",")))
  mir_list <- mir_list[mir_list != ""]   # remove empty strings
  target_mirnas <- unique(mir_list)
  message("Using ", length(target_mirnas), " miRNA(s) from --miRNA: ",
          paste(target_mirnas, collapse = ", "))
}

# From -f file (one per line)
if (!is.null(opt$file)) {
  if (!file.exists(opt$file)) {
    stop("Specified miRNA file does not exist: ", opt$file, call. = FALSE)
  }
  file_mirnas <- readLines(opt$file, warn = FALSE)
  file_mirnas <- trimws(file_mirnas)
  file_mirnas <- file_mirnas[file_mirnas != ""]
  if (is.null(target_mirnas)) {
    target_mirnas <- unique(file_mirnas)
  } else {
    # Merge and deduplicate if both given (union)
    target_mirnas <- unique(c(target_mirnas, file_mirnas))
  }
  message("Total target miRNAs after including file: ", length(target_mirnas))
}

# -------------------------------------------------------------------------
# Load and preprocess expression data
# -------------------------------------------------------------------------
expr_raw <- read.table(opt$input, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
expr_mat <- as.matrix(expr_raw)

# Log2 transformation
expr_log <- log2(expr_mat + 1)

# Check column count (must be even)
nc <- ncol(expr_log)
if (nc %% 2 != 0) {
  stop("Number of samples must be even (two cases with equal replicates).", call. = FALSE)
}
n_rep <- nc / 2

sample_names <- colnames(expr_log)
case1_samples <- sample_names[1:n_rep]
case2_samples <- sample_names[(n_rep + 1):nc]

message("Case1 samples: ", paste(case1_samples, collapse = ", "))
message("Case2 samples: ", paste(case2_samples, collapse = ", "))

# -------------------------------------------------------------------------
# Filter miRNA if specified, else all
# -------------------------------------------------------------------------
if (!is.null(target_mirnas)) {
  # Keep only miRNAs that exist in the data
  valid_mirnas <- intersect(target_mirnas, rownames(expr_log))
  if (length(valid_mirnas) == 0) {
    stop("None of the specified miRNAs were found in the expression matrix.", call. = FALSE)
  }
  if (length(valid_mirnas) < length(target_mirnas)) {
    missing <- setdiff(target_mirnas, rownames(expr_log))
    warning("The following miRNAs were not found and will be ignored: ",
            paste(missing, collapse = ", "), call. = FALSE)
  }
  expr_log <- expr_log[valid_mirnas, , drop = FALSE]
  message("Plotting ", nrow(expr_log), " miRNA(s).")
} else {
  message("Plotting all ", nrow(expr_log), " miRNAs.")
}

# -------------------------------------------------------------------------
# Define plotting function for a single miRNA
# -------------------------------------------------------------------------
plot_mirna <- function(mirna_name, expr_data, case1_s, case2_s) {
  # Extract expression vector for this miRNA
  expr <- expr_data[mirna_name, ]
  
  # Build a data frame
  df <- data.frame(
    sample = names(expr),
    expression = as.numeric(expr),
    stringsAsFactors = FALSE
  )
  df$case <- ifelse(df$sample %in% case1_s, "case1", "case2")
  df$case <- factor(df$case, levels = c("case1", "case2"))
  
  # Create plot
  p <- ggplot(df, aes(x = case, y = expression)) +
    geom_point(position = position_jitter(width = 0.1), size = 3, alpha = 0.8, color = "steelblue") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "red", linewidth = 0.5) +
    labs(title = mirna_name,
         x = "Case",
         y = expression(log[2](expression + 1))) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold.italic", hjust = 0.5),
      axis.text.x = element_text(face = "bold")
    )
  return(p)
}

# -------------------------------------------------------------------------
# Iterate over miRNAs and save individual figures
# -------------------------------------------------------------------------
for (mir in rownames(expr_log)) {
  p <- plot_mirna(mir, expr_log, case1_samples, case2_samples)
  
  # Save as PDF
  pdf_file <- file.path(opt$output, paste0("miRNA_", mir, "_expression.pdf"))
  ggsave(pdf_file, plot = p, width = 5, height = 5, device = "pdf")
  
  # Save as PNG
  png_file <- file.path(opt$output, paste0("miRNA_", mir, "_expression.png"))
  ggsave(png_file, plot = p, width = 5, height = 5, device = "png", dpi = 300)
  
  # Save as SVG
  svg_file <- file.path(opt$output, paste0("miRNA_", mir, "_expression.svg"))
  ggsave(svg_file, plot = p, width = 5, height = 5, device = "svg")
  
  message("Saved plots for: ", mir)
}

message("All done. Output directory: ", opt$output)
