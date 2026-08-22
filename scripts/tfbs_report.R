#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(optparse)
  library(svglite)
  library(igraph)
  library(ggraph)
  library(ggrepel)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "input file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = ".",
              help = "output dir [default: %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file path must be provided (-i/--input)", call. = FALSE)
}

if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
  message("Create output folder: ", opt$output)
}

message("Read input file: ", opt$input)
df <- read.table(opt$input, header = FALSE, stringsAsFactors = FALSE)

df <- df %>%
  mutate(
    start = as.integer(sub("-.*", "", V4)),
    end = as.integer(sub(".*-", "", V4))
  )

# picture1: TFBS location distribution map
message("Draw TFBS location distribution map...")
p1 <- ggplot(df, aes(y = V3)) +
  geom_rect(aes(xmin = start, xmax = end, 
                ymin = as.numeric(factor(V3)) - 0.4, 
                ymax = as.numeric(factor(V3)) + 0.4, 
                fill = V2),
            alpha = 0.7) +
  geom_text(aes(x = (start + end)/2, label = ifelse(V5 == "+", "→", "←")), 
            vjust = 0.5, size = 4) +
  scale_y_discrete(limits = unique(df$V3)) +
  labs(x = "Position (bp)", y = "Target", 
       title = "TFBS Distribution on miRNA Targets", 
       fill = "TF Family") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p1_pdf <- file.path(opt$output, "TFBS_distribution.pdf")
p1_png <- file.path(opt$output, "TFBS_distribution.png")
p1_svg <- file.path(opt$output, "TFBS_distribution.svg")

ggsave(p1_pdf, p1, width = 26, height = 14)
ggsave(p1_png, p1, width = 26, height = 14, dpi = 300)  # 300dpi适合学术发表
ggsave(p1_svg, p1, width = 26, height = 14)

message("Save TFBS distribution map:")
message("  PDF: ", p1_pdf)
message("  PNG: ", p1_png)
message("  SVG: ", p1_svg)

# picture2: TFBS quantity statistics
message("\nDraw a statistical chart of TFBS quantity...")
tfbs_count <- df %>% count(V3, name = "TFBS_Number")

p2 <- ggplot(tfbs_count, aes(x = V3, y = TFBS_Number, fill = V3)) +
  geom_col(alpha = 0.7) +
  geom_text(aes(label = TFBS_Number), vjust = -0.5, size = 4) +
  labs(x = "", y = "Number of TFBS", title = "TFBS Count per miRNA Target") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)
  )

p2_pdf <- file.path(opt$output, "TFBS_count.pdf")
p2_png <- file.path(opt$output, "TFBS_count.png")
p2_svg <- file.path(opt$output, "TFBS_count.svg")

ggsave(p2_pdf, p2, width = 14, height = 26)
ggsave(p2_png, p2, width = 14, height = 26, dpi = 300)
ggsave(p2_svg, p2, width = 14, height = 26)

message("Save TFBS count chart:")
message("  PDF: ", p2_pdf)
message("  PNG: ", p2_png)
message("  SVG: ", p2_svg)

# picture3: TF family statistics
tf_family_count <- df %>% count(V2, name = "TF_Family_Number")

p3 <- ggplot(tf_family_count, aes(x = V2, y = TF_Family_Number, fill = V2)) +
  geom_col(alpha = 0.7) +
  geom_text(aes(label = TF_Family_Number), vjust = -0.5, size = 4) +
  labs(x = "", y = "Number of TFBS", title = "TFBS Count per TF Family") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)
  )

p3_pdf <- file.path(opt$output, "TF_family_count.pdf")
p3_png <- file.path(opt$output, "TF_family_count.png")
p3_svg <- file.path(opt$output, "TF_family_count.svg")

ggsave(p3_pdf, p3, width = 26, height = 14)
ggsave(p3_png, p3, width = 26, height = 14, dpi = 300)
ggsave(p3_svg, p3, width = 26, height = 14)

message("Save TF family count chart:")
message("  PDF: ", p3_pdf)
message("  PNG: ", p3_png)
message("  SVG: ", p3_svg)


# picture4: TF-miRNA interaction network
edges <- df %>%
  select(V2, V3) %>%
  distinct() %>%
  rename(TF = V2, miRNA = V3)

message("Number of unique TF-miRNA pairs: ", nrow(edges))

g <- graph_from_data_frame(edges, directed = FALSE)
V(g)$degree <- degree(g)
V(g)$type <- ifelse(grepl("MIR", V(g)$name, ignore.case = TRUE), "miRNA", "TF")

message("Plotting network...")
set.seed(42)
p4 <- ggraph(g, layout = "fr") +
  geom_edge_link(color = "grey60", alpha = 0.5) +
  geom_node_point(aes(size = degree, color = type), alpha = 0.85) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.5, max.overlaps = 20) +
  scale_color_manual(values = c("TF" = "steelblue", "miRNA" = "tomato"),
                     name = "Node Type") +
  scale_size_continuous(name = "Degree", range = c(2, 12)) +
  labs(title = "TF-miRNA Regulatory Network") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "right"
  )

p4_pdf <- file.path(opt$output, "TF_miRNA_network.pdf")
p4_png <- file.path(opt$output, "TF_miRNA_network.png")
p4_svg <- file.path(opt$output, "TF_miRNA_network.svg")

ggsave(p4_pdf, p4, width = 16, height = 16)
ggsave(p4_png, p4, width = 16, height = 16, dpi = 300)
ggsave(p4_svg, p4, width = 16, height = 16)

message("Save TF-miRNA network chart:")
message("  PDF: ", p4_pdf)
message("  PNG: ", p4_png)
message("  SVG: ", p4_svg)

message("\nAll charts generated successfully!")
