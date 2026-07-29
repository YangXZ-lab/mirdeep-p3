#!/usr/bin/env Rscript

# =========================================================================
# miRNA alignment score classification and visualization
# 
# Input:  -i <score file>   (tab-separated, header, columns: query, reference,
#                              mismatch_penalty, position_penalty, length_penalty,
#                              total_score)
# Output: -o <output dir>   (tables in TSV, combined boxplots in PDF)
# =========================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(gridExtra)
})

# ---------- command line arguments ----------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input score file (TSV)", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = ".",
              help = "Output directory [default: %default]", metavar = "DIR")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Classify miRNA alignments and generate summary plots")
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("--input is required.", call. = FALSE)
}

if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# ---------- helper functions ----------
extract_family <- function(name) {
  # Family is the part before the first hyphen (e.g., "MIR10190" from "MIR10190-isoform1-27")
  strsplit(name, "-")[[1]][1]
}

# Safe distance calculations
calc_diff <- function(high, low) { high - low }

# Create a named list for empty stats
empty_stats <- function() {
  list(max = NA_real_, min = NA_real_)
}

# ---------- read data ----------
cat("Reading input file:", opt$input, "\n")
df <- read.table(opt$input, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                 check.names = FALSE, comment.char = "")
colnames(df) <- c("query", "reference", "mismatch_penalty", "position_penalty",
                  "length_penalty", "total_score")

# Add family columns
df$query_family  <- sapply(df$query, extract_family)
df$ref_family    <- sapply(df$reference, extract_family)
# Flag self-alignment
df$is_self <- df$query == df$reference

# ---------- classify each query ----------
queries <- unique(df$query)
cat("Total unique queries:", length(queries), "\n")

# Initialize storage for categories
cat1 <- list()   # only self-alignments
cat21 <- list()  # only same-family (after removing self)
cat22 <- list()  # only different-family (after removing self)
cat23 <- list()  # both same and different family (after removing self)

for (q in queries) {
  q_data <- df[df$query == q, ]
  # Separate self and non-self
  non_self <- q_data[!q_data$is_self, ]
  self_count <- nrow(q_data[q_data$is_self, ])
  
  # If no non-self alignments, it belongs to category 1
  if (nrow(non_self) == 0) {
    cat1[[q]] <- list(name = q, family = q_data$query_family[1])
    next
  }
  
  # Families present in non-self
  fams <- unique(non_self$ref_family)
  query_fam <- q_data$query_family[1]
  same_fam <- query_fam %in% fams
  other_fam <- any(fams != query_fam)
  
  if (same_fam && !other_fam) {
    # Only same family
    same_scores <- non_self[non_self$ref_family == query_fam, "total_score"]
    max_same <- max(same_scores)
    min_same <- if (length(same_scores) == 1) 0 else min(same_scores)
    cat21[[q]] <- list(name = q, family = query_fam,
                        max_same = max_same, min_same = min_same,
                        n_same = length(same_scores))
  } else if (!same_fam && other_fam) {
    # Only different families
    diff_scores <- non_self$total_score  # all are different family
    max_diff <- max(diff_scores)
    min_diff <- if (length(diff_scores) == 1) 0 else min(diff_scores)
    diff_fams <- sort(unique(non_self$ref_family))
    diff_fams_str <- paste(diff_fams, collapse = ",")
    cat22[[q]] <- list(name = q, family = query_fam,
                        diff_families = diff_fams_str,
                        max_diff = max_diff, min_diff = min_diff,
                        n_diff = length(diff_scores))
  } else {
    # Both same and different
    same_scores <- non_self[non_self$ref_family == query_fam, "total_score"]
    diff_non_self <- non_self[non_self$ref_family != query_fam, ]
    diff_scores <- diff_non_self$total_score
    diff_fams <- sort(unique(diff_non_self$ref_family))
    diff_fams_str <- paste(diff_fams, collapse = ",")
    
    max_same <- max(same_scores)
    min_same <- if (length(same_scores) == 1) max_same else min(same_scores)
    max_diff <- max(diff_scores)
    min_diff <- if (length(diff_scores) == 1) max_diff else min(diff_scores)
    
    cat23[[q]] <- list(name = q, family = query_fam,
                        max_same = max_same, min_same = min_same,
                        diff_families = diff_fams_str,
                        max_diff = max_diff, min_diff = min_diff,
                        n_same = length(same_scores), n_diff = length(diff_scores))
  }
}

# ---------- Write tables ----------
# Category 1
if (length(cat1) > 0) {
  cat1_df <- do.call(rbind, lapply(cat1, function(x) data.frame(name = x$name, family = x$family, stringsAsFactors = FALSE)))
  write.table(cat1_df, file = file.path(opt$output, "category1_only_self.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Category 2.1
if (length(cat21) > 0) {
  cat21_df <- do.call(rbind, lapply(cat21, function(x) data.frame(name = x$name, family = x$family,
                                                                   max_same = x$max_same, min_same = x$min_same,
                                                                   stringsAsFactors = FALSE)))
  write.table(cat21_df, file = file.path(opt$output, "category21_same_family.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Category 2.2
if (length(cat22) > 0) {
  cat22_df <- do.call(rbind, lapply(cat22, function(x) data.frame(name = x$name, family = x$family,
                                                                   diff_families = x$diff_families,
                                                                   max_diff = x$max_diff, min_diff = x$min_diff,
                                                                   stringsAsFactors = FALSE)))
  write.table(cat22_df, file = file.path(opt$output, "category22_diff_family.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Category 2.3
if (length(cat23) > 0) {
  cat23_df <- do.call(rbind, lapply(cat23, function(x) data.frame(name = x$name, family = x$family,
                                                                   max_same = x$max_same, min_same = x$min_same,
                                                                   diff_families = x$diff_families,
                                                                   max_diff = x$max_diff, min_diff = x$min_diff,
                                                                   stringsAsFactors = FALSE)))
  write.table(cat23_df, file = file.path(opt$output, "category23_mixed.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# ---------- Boxplot helper ----------
# Create a combined boxplot (2x2) for a category
make_combined_boxplot <- function(data_list, title_main, output_file) {
  # data_list is a list of named vectors, each a separate boxplot
  # We'll create a ggplot for each and combine with grid.arrange
  plots <- list()
  for (i in seq_along(data_list)) {
    var_name <- names(data_list)[i]
    vals <- data_list[[i]]
    if (length(vals) == 0) next
    df_temp <- data.frame(value = vals, group = var_name)
    p <- ggplot(df_temp, aes(x = group, y = value)) +
      geom_boxplot(fill = "#69b3a2", width = 0.3) +
      labs(title = var_name, y = "Score", x = "") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
    plots[[length(plots) + 1]] <- p
  }
  if (length(plots) == 0) return()
  # Arrange in a grid, maximum 4 plots
  n <- length(plots)
  ncol <- min(2, n)
  nrow <- ceiling(n / ncol)
  combined <- do.call(grid.arrange, c(plots, list(ncol = ncol, nrow = nrow, top = title_main)))
  ggsave(filename = output_file, plot = combined, device = "pdf", width = 12, height = 10)
}

# ---------- Generate plots ----------

# --- Category 2.1 plots ---
if (length(cat21) > 0) {
  # Split into two groups: single same-family alignment vs multiple
  single_align <- sapply(cat21, function(x) x$n_same == 1)
  multi_align <- !single_align
  
  # Scores for singles
  single_scores <- sapply(cat21[single_align], function(x) x$max_same)
  # For multiples: max, min, difference
  multi_max <- sapply(cat21[multi_align], function(x) x$max_same)
  multi_min <- sapply(cat21[multi_align], function(x) x$min_same)
  multi_diff <- multi_max - multi_min
  
  plot_data <- list()
  if (length(single_scores) > 0) plot_data[["Single alignment score"]] <- single_scores
  if (length(multi_max) > 0) {
    plot_data[["Multiple: max score"]] <- multi_max
    plot_data[["Multiple: min score"]] <- multi_min
    plot_data[["Multiple: max-min difference"]] <- multi_diff
  }
  if (length(plot_data) > 0) {
    make_combined_boxplot(plot_data, "Category 2.1: Same-family only",
                          file.path(opt$output, "category21_boxplot.pdf"))
  }
}

# --- Category 2.2 plots ---
if (length(cat22) > 0) {
  single_align <- sapply(cat22, function(x) x$n_diff == 1)
  multi_align <- !single_align
  
  single_scores <- sapply(cat22[single_align], function(x) x$max_diff)
  multi_max <- sapply(cat22[multi_align], function(x) x$max_diff)
  multi_min <- sapply(cat22[multi_align], function(x) x$min_diff)
  multi_diff <- multi_max - multi_min
  
  plot_data <- list()
  if (length(single_scores) > 0) plot_data[["Single alignment score"]] <- single_scores
  if (length(multi_max) > 0) {
    plot_data[["Multiple: max score"]] <- multi_max
    plot_data[["Multiple: min score"]] <- multi_min
    plot_data[["Multiple: max-min difference"]] <- multi_diff
  }
  if (length(plot_data) > 0) {
    make_combined_boxplot(plot_data, "Category 2.2: Different-family only",
                          file.path(opt$output, "category22_boxplot.pdf"))
  }
}

# --- Category 2.3 plots ---
if (length(cat23) > 0) {
  max_same <- sapply(cat23, function(x) x$max_same)
  min_same <- sapply(cat23, function(x) x$min_same)
  max_diff <- sapply(cat23, function(x) x$max_diff)
  min_diff <- sapply(cat23, function(x) x$min_diff)
  
  diff_same_max_vs_diff_max <- max_same - max_diff
  diff_same_min_vs_diff_max <- min_same - max_diff   # as described
  
  plot_data <- list()
  if (length(max_same) > 0) plot_data[["Same-family max score"]] <- max_same
  if (length(min_same) > 0) plot_data[["Same-family min score"]] <- min_same
  if (length(diff_same_max_vs_diff_max) > 0) plot_data[["(Same max) - (Diff max)"]] <- diff_same_max_vs_diff_max
  if (length(diff_same_min_vs_diff_max) > 0) plot_data[["(Same min) - (Diff max)"]] <- diff_same_min_vs_diff_max
  
  if (length(plot_data) > 0) {
    make_combined_boxplot(plot_data, "Category 2.3: Mixed families",
                          file.path(opt$output, "category23_boxplot.pdf"))
  }
}

cat("All tables and plots have been saved to", opt$output, "\n")
