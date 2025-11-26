#!/usr/bin/env Rscript

# --- Load packages ---
if (!require("optparse")) install.packages("optparse", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")
if (!require("tidyverse")) install.packages("tidyverse", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")
if (!require("ggpmisc")) install.packages("ggpmisc", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")
if (!require("ggpubr")) install.packages("ggpubr", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")

library(optparse)
library(tidyverse)
library(ggpmisc)
library(ggpubr)

# --- Define command-line options ---
option_list <- list(
  make_option(c("-i", "--otu_table"), type = "character", default = "result/otutab.txt",
              help = "Input OTU table [default: %default]"),
  make_option(c("-a", "--metadata"), type = "character", default = "result/metadata.txt",
              help = "Metadata file [default: %default]"),
  make_option(c("-o", "--output"), type = "character", default = "result2/tax/",
              help = "Output PDF file [default: %default]"),
  make_option(c("--output_png"), type = "character", default = NULL,
              help = "Optional: Output PNG file"),
  make_option(c("--filtered_otu"), type = "character", default = NULL,
              help = "Optional: Save filtered OTU table as TSV"),
  make_option(c("--filtered_metadata"), type = "character", default = NULL,
              help = "Optional: Save filtered metadata as TSV"),
  make_option(c("--core_asvs"), type = "character", default = NULL,
              help = "Optional: Save core ASVs as TSV"),
  make_option(c("--plot_data"), type = "character", default = NULL,
              help = "Optional: Save merged plot data as TSV")
)
opts <- parse_args(OptionParser(option_list = option_list))

# --- Read input files ---
otu <- read.table(opts$otu_table, header=TRUE, row.names=1, sep="\t", check.names=FALSE, comment.char = "")
metadata <- read.table(opts$metadata, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Transpose OTU table to long format
otu_long <- otu %>%
  rownames_to_column("ASV") %>%
  pivot_longer(-ASV, names_to = "SampleID", values_to = "Abundance")

# Merge with metadata
otu_long <- left_join(otu_long, metadata, by = "SampleID")

# Check how many samples per group
cat("Number of samples per group:\n")
print(table(otu_long$Group))

# Calculate core ASVs (present in ≥80% of samples per group)
core_asvs <- otu_long %>%
  filter(Abundance > 0) %>%
  group_by(Group, ASV) %>%
  summarise(present_in = n_distinct(SampleID), .groups = "drop") %>%
  left_join(otu_long %>% select(Group, SampleID) %>% distinct() %>% count(Group, name = "total"), by = "Group") %>%
  mutate(frequency = present_in / total) %>%
  filter(frequency >= 0.8) %>%
  select(Group, ASV)

# Save core ASVs if path provided
if (!is.null(opts$core_asvs)) {
  write.table(core_asvs, opts$core_asvs, sep = "\t", row.names = FALSE, quote = FALSE)
}

# Count how many core ASVs per group
cat("Number of core ASVs per group:\n")
print(core_asvs %>% group_by(Group) %>% summarise(n_core_ASVs = n()))

# Filter core ASV abundance
core_abundance <- semi_join(otu_long, core_asvs, by = c("Group", "ASV")) %>%
  group_by(SampleID, Group) %>%
  summarise(core_abundance = sum(Abundance), .groups = "drop")

# Total abundance of all ASVs
total_abundance <- otu_long %>%
  group_by(SampleID, Group) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop")

# Combine for plotting
plot_df <- left_join(core_abundance, total_abundance, by = c("SampleID", "Group"))

# Save plot data if path provided
if (!is.null(opts$plot_data)) {
  write.table(plot_df, opts$plot_data, sep = "\t", row.names = FALSE, quote = FALSE)
}

# Save filtered OTU and metadata if requested
if (!is.null(opts$filtered_otu)) {
  write.table(otu, opts$filtered_otu, sep = "\t", quote = FALSE)
}
if (!is.null(opts$filtered_metadata)) {
  write.table(metadata, opts$filtered_metadata, sep = "\t", row.names = FALSE, quote = FALSE)
}

# Check how many points are being plotted
cat("Number of points to plot:", nrow(plot_df), "\n")

# Generate and save plots
p <- ggplot(plot_df, aes(x = log10(total_abundance + 1), y = log10(core_abundance + 1), color = Group)) +
  geom_jitter(size = 1.2, alpha = 0.8, width = 0.02, height = 0.02) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +
  facet_wrap(~Group, scales = "free") +
  stat_poly_eq(
    aes(label = paste(..rr.label..)),
    formula = y ~ x,
    parse = TRUE,
    size = 3,
    label.x.npc = "left",
    label.y.npc = 0.95,
    color = "black"
  ) +
  labs(
    title = "Core vs Total ASV Abundance (log-transformed)",
    x = expression(paste("Log"[10], " of Total ASV Abundance + 1")),
    y = expression(paste("Log"[10], " of Core ASV Abundance + 1"))
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "none",
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )

# Save plots with dimensions and units matching previous usage
ggsave(paste0(opts$output, "core_vs_total_abundance_plot.pdf"), p, width = 100, height = 100, units = "mm")
if (!is.null(opts$output_png)) {
  ggsave(paste0(opts$output_png, "core_vs_total_abundance_plot.png"), p, width = 100, height = 100, units = "mm")
}
