#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Script: basic_stat.R
# Description:
#   Generate four visualization plots from a miRNA basic-info file:
#     1. Sequence length distribution (pie chart)
#     2. First nucleotide composition by sequence length (stacked bar)
#     3. Relative nucleotide composition by position (stacked bar)
#     4. miRNA family abundance (bar plot)
#   Each plot is saved as PDF, SVG, and PNG in the specified output directory.
#
# Usage:
#   Rscript basic_stat.R -i <input_file> -o <output_dir>
# -----------------------------------------------------------------------------

# Load required packages
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(stringr)
})

# ---- Command-line arguments -------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to the input basic-info file (tab-separated) [required]",
              metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = "./output",
              help = "Output directory for saving plots [default: %default]",
              metavar = "DIR")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Visualize miRNA annotation from basic-info table.")
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file must be provided via -i/--input.")
}

input_file <- opt$input
output_dir <- opt$output

# Create output directory if needed
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---- Data loading and unified preprocessing --------------------------------
# Define expected column names
col_names <- c("MIRNA", "MIRNA_ID", "MIRNA_Family", "Species",
               "Chromosome", "Strand", "Precursor_20nt_start", "Precursor_20nt_end",
               "Precursor_20nt_seq", "Precursor_20nt_struc", "Precursor_start",
               "Precursor_end", "Precursor_seq", "Precursor_struc",
               "Mature", "Mature_start", "Mature_end", "Mature_seq",
               "Star", "Star_start", "Star_end", "Star_seq")

raw_data <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(raw_data) <- col_names

# Preprocess: trim, compute length/first base, keep only canonical DNA sequences
processed_data <- raw_data %>%
  mutate(
    Mature_seq = str_trim(Mature_seq),
    seq_length = nchar(Mature_seq),
    first_base = substr(Mature_seq, 1, 1)
  ) %>%
  filter(
    seq_length >= 1,
    !str_detect(Mature_seq, "[^ATCG]")  # exclude sequences with non-canonical bases
  )

total_sequences <- nrow(processed_data)

# ---- Color palettes ---------------------------------------------------------
pie_colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
                "#8491B4", "#91D1C2", "#DC0000")
base_colors <- c("A" = "#74CE98", "T" = "#E28B7C", "C" = "#6A9FB5", "G" = "#F2BE3C")
family_bar_color <- "#74CE98"

# ---- Helper to save a plot in three formats --------------------------------
save_plot <- function(plot_obj, prefix, w = 8, h = 6) {
  ggsave(paste0(prefix, ".pdf"), plot_obj, width = w, height = h, device = "pdf")
  ggsave(paste0(prefix, ".svg"), plot_obj, width = w, height = h, device = "svg")
  ggsave(paste0(prefix, ".png"), plot_obj, width = w, height = h, device = "png")
}

# =============================================================================
# Plot 1: Sequence length distribution (pie chart)
# =============================================================================
length_dist <- processed_data %>%
  count(seq_length, name = "count") %>%
  arrange(desc(count)) %>%
  mutate(
    seq_length = factor(seq_length),
    percent = count / sum(count),
    label = percent(percent, accuracy = 0.1)
  )

p1 <- ggplot(length_dist, aes(x = "", y = count, fill = seq_length)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 4, color = "white", fontface = "bold") +
  scale_fill_manual(values = pie_colors) +
  labs(
    title = "Statistics of sequence length distribution",
    fill = "Length (bp)",
    caption = paste0("Total sequences：", format(total_sequences, big.mark = ","))
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold",
                              margin = margin(b = 5)),
    plot.caption = element_text(hjust = 0.5, size = 11, color = "black",
                                margin = margin(t = 5)),
    legend.title = element_text(size = 12, face = "bold", margin = margin(b = 2)),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.8, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    legend.margin = margin(l = -10),
    legend.position = "right",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
    plot.background = element_rect(fill = "white", color = NA)
  )

save_plot(p1, file.path(output_dir, "length_dist"))

# =============================================================================
# Plot 2: First nucleotide composition by sequence length (stacked bar)
# =============================================================================
base_by_length <- processed_data %>%
  count(seq_length, first_base, name = "count") %>%
  group_by(seq_length) %>%
  mutate(
    total = sum(count),
    percent = count / total
  ) %>%
  ungroup() %>%
  mutate(seq_length = factor(seq_length))

p2 <- ggplot(base_by_length, aes(x = seq_length, y = percent, fill = first_base)) +
  geom_col(position = "stack", color = "white") +
  geom_text(aes(label = percent(percent, accuracy = 0.1)),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = base_colors, name = "First base") +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "First nucleotide composition by sequence length",
    x = "Sequence length (bp)",
    y = "Percentage of sequences",
    caption = paste0("Total sequences：", format(total_sequences, big.mark = ","))
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold",
                              margin = margin(b = 10)),
    plot.caption = element_text(hjust = 0.5, size = 11, color = "black",
                                margin = margin(t = 10)),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

save_plot(p2, file.path(output_dir, "first_base_dist"))

# =============================================================================
# Plot 3: Relative nucleotide composition by position
# =============================================================================
max_len <- max(processed_data$seq_length)

# Efficiently tabulate base counts at each position
position_list <- lapply(1:max_len, function(pos) {
  idx <- processed_data$seq_length >= pos
  if (sum(idx) == 0) return(NULL)
  bases <- substr(processed_data$Mature_seq[idx], pos, pos)
  counts <- table(bases)
  data.frame(
    position = pos,
    base = names(counts),
    count = as.integer(counts),
    total = sum(idx),
    stringsAsFactors = FALSE
  )
})
position_data <- do.call(rbind, position_list[lengths(position_list) > 0])

position_data <- position_data %>%
  mutate(
    percent = count / total,
    base = factor(base, levels = c("A", "T", "C", "G"))
  )

p3 <- ggplot(position_data, aes(x = position, y = percent, fill = base)) +
  geom_col(position = "stack", color = "white") +
  geom_text(
    aes(label = ifelse(percent > 0.1, percent(percent, accuracy = 1), "")),
    position = position_stack(vjust = 0.5),
    size = 3, color = "white", fontface = "bold"
  ) +
  scale_fill_manual(values = base_colors, name = "Base") +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_x_continuous(breaks = 1:max_len) +
  labs(
    title = "Relative nucleotide composition by position",
    x = "Position in sequence",
    y = "Percentage of reads",
    caption = paste0("Total sequences：", format(total_sequences, big.mark = ","))
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold",
                              margin = margin(b = 10)),
    plot.caption = element_text(hjust = 0.5, size = 11, color = "black",
                                margin = margin(t = 10)),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

save_plot(p3, file.path(output_dir, "base_dist"), w = 10)

# =============================================================================
# Plot 4: miRNA family abundance
# =============================================================================
family_counts <- processed_data %>%
  mutate(family = str_trim(MIRNA_Family)) %>%
  filter(!is.na(family), family != "") %>%
  count(family, name = "count") %>%
  arrange(desc(count))

p4 <- ggplot(family_counts, aes(x = reorder(family, -count), y = count)) +
  geom_col(fill = family_bar_color, color = "white") +
  geom_text(
    aes(label = count),
    vjust = -0.3,
    size = 3.5,
    fontface = "bold"
  ) +
  labs(
    title = "miRNA family abundance distribution",
    x = "miRNA family",
    y = "Number of sequences",
    caption = paste0("Total sequences：",
                     format(sum(family_counts$count), big.mark = ","))
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold",
                              margin = margin(b = 10)),
    plot.caption = element_text(hjust = 0.5, size = 11, color = "black",
                                margin = margin(t = 10)),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    plot.margin = margin(t = 25, r = 10, b = 10, l = 10)
  ) +
  ylim(0, max(family_counts$count) * 1.1)

save_plot(p4, file.path(output_dir, "miRNA_family"), w = 12)

# Final message
message("All plots successfully saved in ", normalizePath(output_dir))
