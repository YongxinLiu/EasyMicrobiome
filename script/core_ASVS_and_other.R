#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(tibble)
  library(ggplot2)
  library(dplyr)
})

# Error handling function
options(error = function() {
  cat("An error occurred:\n")
  traceback()
  quit(save = "no", status = 1)
})

# Define command-line options
option_list <- list(
  make_option(c("-i", "--otu_table"), type = "character", default = "result/otutab.txt",
              help = "Path to OTU table file [default: %default]"),
  make_option(c("-d", "--metadata"), type = "character", default = "result/metadata.txt",
              help = "Path to metadata file [default: %default]"),
  make_option(c("-t", "--taxonomy"), type = "character", default = "result/taxonomy.txt",
              help = "Path to taxonomy file [default: %default]"),
  make_option(c("-o", "--output_dir"), type = "character", default = "result/coreASVs2",
              help = "Output directory for plots and results [default: %default]"),
  make_option(c("-g", "--group"), type = "character", default = "Group",
              help = "Metadata column name for grouping [default: %default]"),
  make_option(c("-n", "--ntaxa"), type = "integer", default = 8,
              help = "Number of top taxa to display [default: %default]")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# Create output directory if it doesn't exist
if (!dir.exists(opts$output_dir)) {
  dir.create(opts$output_dir, recursive = TRUE)
}

# Read input files
otu_table <- read.table(opts$otu_table, header = TRUE, row.names = 1, sep = "\t",
                        check.names = FALSE, comment.char = "")
metadata <- read.table(opts$metadata, header = TRUE, row.names = 1, sep = "\t",
                       check.names = FALSE, comment.char = "")
taxonomy <- read.table(opts$taxonomy, header = FALSE, row.names = 1, sep = "\t",
                       check.names = FALSE, comment.char = "")

# Remove taxonomy header row if present as a row name
if ("OTUID" %in% rownames(taxonomy)) {
  taxonomy <- taxonomy[rownames(taxonomy) != "OTUID", ]
}

colnames(taxonomy) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Check matching samples between OTU table and metadata
common_samples <- intersect(colnames(otu_table), rownames(metadata))
if (length(common_samples) == 0) {
  stop("No matching sample IDs between OTU table and metadata!")
}

# Subset OTU table and metadata to common samples
otu_table <- otu_table[, common_samples]
metadata  <- metadata[common_samples, ]

# Identify core ASVs: present in at least 80% of samples
presence_threshold <- 0.8
sample_count <- ncol(otu_table)
core_asvs <- rowSums(otu_table > 0) >= (presence_threshold * sample_count)

core_otu_table <- otu_table[core_asvs, ]
other_otu_table <- otu_table[!core_asvs, ]

# Calculate relative abundance (%)
rel_abundance <- sweep(otu_table, 2, colSums(otu_table), FUN = "/") * 100

core_rel_abundance <- colSums(rel_abundance[core_asvs, ])
other_rel_abundance <- colSums(rel_abundance[!core_asvs, ])

# Summarize relative abundance across samples (mean)
relative_abundance_df <- tibble(
  group = c("Core ASVs", "Other ASVs"),
  abundance = c(mean(core_rel_abundance), mean(other_rel_abundance))
)

# Calculate absolute abundance (sum of counts) for core and other ASVs
core_abs_abundance <- sum(core_otu_table)
other_abs_abundance <- sum(other_otu_table)

absolute_abundance_df <- tibble(
  group = c("Core ASVs", "Other ASVs"),
  abundance = c(core_abs_abundance, other_abs_abundance)
)

# Plotting function for pie charts
plot_piechart <- function(df, title, filename, fill_colors = c("#D7263D", "#3F88C5")) {
  p <- ggplot(df, aes(x = "", y = abundance, fill = group)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    geom_text(aes(label = paste0(round(abundance, 1), ifelse(max(df$abundance) <= 1, "", "%"))),
              position = position_stack(vjust = 0.5), size = 5) +
    scale_fill_manual(values = fill_colors) +
    theme_void() +
    labs(title = title)
  
  ggsave(filename = filename, plot = p, width = 6, height = 6)
  cat("Plot saved to:", filename, "\n")
}

# Plot and save relative abundance pie chart
relative_plot_file <- file.path(opts$output_dir, "core_vs_other_abundance_piechart.png")
plot_piechart(relative_abundance_df, "Relative Abundance of Core vs Other ASVs", relative_plot_file)

# Plot and save absolute abundance pie chart
absolute_plot_file <- file.path(opts$output_dir, "absolute_abundance_piechart.png")
plot_piechart(absolute_abundance_df, "Absolute Abundance of Core vs Other ASVs", absolute_plot_file)

# Extract taxonomy for core ASVs
core_taxonomy <- taxonomy[rownames(core_otu_table), ]

# Save core OTU table and taxonomy
write.table(core_otu_table, file = file.path(opts$output_dir, "core_otu_table.txt"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(core_taxonomy, file = file.path(opts$output_dir, "core_taxonomy.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

cat("Core OTU table and taxonomy saved in:", opts$output_dir, "\n")

# Save other ASVs OTU table
write.table(other_otu_table, file = file.path(opts$output_dir, "other_otu_table.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

cat("Other OTU table saved in:", opts$output_dir, "\n")

cat("All processing completed successfully.\n")
