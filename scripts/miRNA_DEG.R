#!/usr/bin/env Rscript

# =========================================================================
# DESeq2 differential expression analysis for miRNA count data
# Produces: results table (with miRNA column name), up/down gene lists,
#           summary, volcano plot
# =========================================================================

args <- commandArgs(trailingOnly = TRUE)
tmp_idx <- which(args == "--tmp")
if (length(tmp_idx) == 1 && tmp_idx < length(args)) {
  custom_tmp <- args[tmp_idx + 1]
  if (!dir.exists(custom_tmp)) {
    dir.create(custom_tmp, recursive = TRUE)
  }
  Sys.setenv(TMPDIR = custom_tmp)
  Sys.setenv(TMP = custom_tmp)
  Sys.setenv(TEMP = custom_tmp)
}

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
})

# -------------------------------------------------------------------------
# Command-line arguments
# -------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input count matrix (tab-separated, header, row.names=1)", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = "./",
              help = "Output directory for results [default: %default]", metavar = "DIR"),
  make_option(c("-r", "--replicate"), type = "integer", default = NULL,
              help = "Number of replicates per condition", metavar = "INT"),
  make_option(c("-p", "--pvalue"), type = "double", default = 0.05,
              help = "Adjusted p-value cutoff [default: %default]"),
  make_option(c("-l", "--logfc"), type = "double", default = 2,
              help = "Fold change cutoff (non-log) [default: %default, i.e. |log2FC| > 1]"),
  make_option("--case1", type = "character", default = "case1",
              help = "Name of baseline condition [default: %default]"),
  make_option("--case2", type = "character", default = "case2",
              help = "Name of treatment condition [default: %default]"),
  make_option("--tmp", type = "character", default = NULL,
              help = "Custom temporary directory (useful when /tmp is not writable)", metavar = "DIR")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "DESeq2-based differential expression analysis for miRNA counts")
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$replicate)) {
  print_help(opt_parser)
  stop("Both --input and --replicate must be specified.", call. = FALSE)
}

if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# -------------------------------------------------------------------------
# Override TMPDIR again (in case option is given but not caught earlier)
# -------------------------------------------------------------------------
if (!is.null(opt$tmp)) {
  if (!dir.exists(opt$tmp)) {
    dir.create(opt$tmp, recursive = TRUE)
  }
  Sys.setenv(TMPDIR = opt$tmp)
  Sys.setenv(TMP = opt$tmp)
  Sys.setenv(TEMP = opt$tmp)
}

# -------------------------------------------------------------------------
# Load and validate count matrix
# -------------------------------------------------------------------------
countData <- as.matrix(read.table(opt$input, header = TRUE, row.names = 1,
                                  sep = "\t", check.names = FALSE))
# Ensure integer counts
countData <- round(countData)
mode(countData) <- "integer"

expected_cols <- 2 * opt$replicate
if (ncol(countData) != expected_cols) {
  stop(sprintf("Expected %d columns (%d replicates × 2 conditions), but found %d.",
               expected_cols, opt$replicate, ncol(countData)), call. = FALSE)
}

# -------------------------------------------------------------------------
# Build sample metadata
# -------------------------------------------------------------------------
sample_names <- paste0(rep(c(opt$case1, opt$case2), each = opt$replicate),
                       "_", rep(1:opt$replicate, 2))
colnames(countData) <- sample_names
group <- gsub("_[0-9]+$", "", sample_names)
colData <- data.frame(row.names = sample_names,
                      condition = factor(group, levels = c(opt$case1, opt$case2)))

message("Conditions: ", opt$case1, " (baseline) vs. ", opt$case2, " (treatment)")
message("Replicates per group: ", opt$replicate)

# -------------------------------------------------------------------------
# DESeq2 workflow
# -------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)
# Use poscounts to handle prevalent zero counts in miRNA-seq data
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds)

# ---------- lfcShrink with apeglm (fallback to normal) ----------
res_names <- resultsNames(dds)
expected_coef <- paste0("condition_", opt$case2, "_vs_", opt$case1)

if (expected_coef %in% res_names) {
  target_coef <- expected_coef
  message("Using coefficient: ", target_coef)
} else {
  target_coef <- res_names[2]  # fallback to the second coefficient
  warning(sprintf("Coefficient '%s' not found. Using '%s' instead.",
                  expected_coef, target_coef))
}

if (requireNamespace("apeglm", quietly = TRUE)) {
  res <- lfcShrink(dds, coef = target_coef, type = "apeglm")
  message("lfcShrink with apeglm applied.")
} else {
  warning("Package 'apeglm' not found. Using type='normal' for lfcShrink.")
  res <- lfcShrink(dds, coef = target_coef, type = "normal")
}
# ----------------------------------------------------------------

res <- res[order(res$padj), ]   # sort by adjusted p-value
DEG <- as.data.frame(res)
DEG <- na.omit(DEG)             # remove genes with NA padj

# -------------------------------------------------------------------------
# Define significance thresholds
# -------------------------------------------------------------------------
logFC_cutoff <- log2(opt$logfc)
padj_cutoff  <- opt$pvalue

k_up   <- DEG$padj < padj_cutoff & DEG$log2FoldChange >  logFC_cutoff
k_down <- DEG$padj < padj_cutoff & DEG$log2FoldChange < -logFC_cutoff
DEG$change <- ifelse(k_up, "UP", ifelse(k_down, "DOWN", "NOT"))

n_up   <- sum(k_up)
n_down <- sum(k_down)
n_ns   <- nrow(DEG) - n_up - n_down

# -------------------------------------------------------------------------
# Write result files
# -------------------------------------------------------------------------
comp_label <- paste0(opt$case1, "_", opt$case2)

# Convert rownames to a column named "miRNA" for clean output
DEG_out <- data.frame(miRNA = rownames(DEG), DEG, check.names = FALSE, stringsAsFactors = FALSE)

# Full results table
write.table(DEG_out,
            file = file.path(opt$output, paste0(comp_label, "_res.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Up‑regulated miRNA list
write.table(data.frame(miRNA = rownames(DEG)[k_up]),
            file = file.path(opt$output, paste0(comp_label, "_up.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Down‑regulated miRNA list
write.table(data.frame(miRNA = rownames(DEG)[k_down]),
            file = file.path(opt$output, paste0(comp_label, "_down.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Summary statistics
summary_text <- c(
  paste0("Comparison: ", opt$case2, " vs ", opt$case1),
  paste0("Total miRNAs after filtering: ", nrow(DEG)),
  paste0("Up-regulated (padj < ", padj_cutoff, ", log2FC > ", round(logFC_cutoff, 2), "): ", n_up),
  paste0("Down-regulated (padj < ", padj_cutoff, ", log2FC < ", -round(logFC_cutoff, 2), "): ", n_down),
  paste0("Not significant: ", n_ns)
)
writeLines(summary_text, file.path(opt$output, paste0(comp_label, "_summary.txt")))
message(paste(summary_text, collapse = "\n"))

# -------------------------------------------------------------------------
# Volcano plot
# -------------------------------------------------------------------------
DEG$neg_log10_padj <- -log10(DEG$padj)

# Label top 10 miRNAs with smallest padj
top_genes <- rownames(DEG)[1:min(10, nrow(DEG))]
# Use NA instead of empty strings to avoid repelling blank labels
DEG$label <- ifelse(rownames(DEG) %in% top_genes, rownames(DEG), NA_character_)

# Define color mapping
change_levels <- c("DOWN", "NOT", "UP")
DEG$change <- factor(DEG$change, levels = change_levels)
color_values <- c("DOWN" = "blue", "NOT" = "grey60", "UP" = "red")

p <- ggplot(DEG, aes(x = log2FoldChange, y = neg_log10_padj, color = change)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(values = color_values, drop = FALSE) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff),
             linetype = "dashed", color = "darkgrey", linewidth = 0.5) +
  geom_hline(yintercept = -log10(padj_cutoff),
             linetype = "dashed", color = "darkgrey", linewidth = 0.5) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20,
                  show.legend = FALSE, box.padding = 0.3, na.rm = TRUE) +
  labs(title = paste0("Volcano plot: ", opt$case2, " vs. ", opt$case1),
       x = expression(log[2]~Fold~Change),
       y = expression(-log[10]~adjusted~italic(p))) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Save in multiple formats
for (fmt in c("pdf", "png", "svg")) {
  ggsave(file.path(opt$output, paste0("volcano_", comp_label, ".", fmt)),
         plot = p, width = 8, height = 7, device = fmt)
  message("Volcano plot saved: ", comp_label, ".", fmt)
}

message("All results written to: ", opt$output)
