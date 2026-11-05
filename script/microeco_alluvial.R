#!/usr/bin/env Rscript

# ---------------------------
# Script function: Generate Alluvial Plot (Genus level) using microeco
# Input: OTU table, taxonomy table, metadata file
# Output: Alluvial plot in PDF format
# ---------------------------

options(warn = -1)  # Suppress warnings

# ---------------------------
# 1. Load required libraries
# ---------------------------
site = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"

pkg_list <- c("optparse", "microeco", "ggplot2", "RColorBrewer", "ggalluvial")
for (p in pkg_list) {
  if (!suppressWarnings(require(p, character.only = TRUE))) {
    install.packages(p, repos = site)
    library(p, character.only = TRUE)
  }
}

# ---------------------------
# 2. Command line arguments
# ---------------------------
option_list = list(
  make_option(c("-o", "--otu"), type = "character", default = "result/otutab.txt",
              help = "OTU table [default %default]"),
  make_option(c("-t", "--taxonomy"), type = "character", default = "result/taxonomy.txt",
              help = "Taxonomy table [default %default]"),
  make_option(c("-m", "--metadata"), type = "character", default = "result/metadata.txt",
              help = "Metadata file [default %default]"),
  make_option(c("-g", "--group"), type = "character", default = "Group",
              help = "Group column name in metadata [default %default]"),
  make_option(c("-n", "--ntaxa"), type = "integer", default = 8,
              help = "Number of top taxa to show [default %default]"),
  make_option(c("-r", "--taxrank"), type = "character", default = "Genus",
              help = "Taxonomic rank to analyze [default %default]"),
  make_option(c("-w", "--width"), type = "numeric", default = 100,
              help = "Figure width in mm [default %default]"),
  make_option(c("-e", "--height"), type = "numeric", default = 80,
              help = "Figure height in mm [default %default]"),
  make_option(c("-x", "--xangle"), type = "numeric", default = 30,
              help = "Angle of x-axis labels [default %default]"),
  make_option(c("-s", "--xsize"), type = "numeric", default = 3,
              help = "Font size of x-axis labels [default %default]"),
  make_option(c("-f", "--output"), type = "character", default = "result/alluvial_microeco_plot.pdf",
              help = "Output PDF file path [default %default]")
)
opts = parse_args(OptionParser(option_list = option_list))

# ---------------------------
# 3. Read input files
# ---------------------------
otu_tab <- read.table(opts$otu, header = TRUE, row.names = 1, sep = "\t", comment.char = "")
tax_tab <- read.table(opts$taxonomy, header = TRUE, row.names = 1, sep = "\t", comment.char = "")
meta_tab <- read.table(opts$metadata, header = TRUE, row.names = 1, sep = "\t", comment.char = "", stringsAsFactors = FALSE)

# ---------------------------
# 4. Create dataset object
# ---------------------------
dataset <- microtable$new(otu_table = otu_tab, tax_table = tax_tab, sample_table = meta_tab)
dataset$tidy_dataset()

# ---------------------------
# 5. Generate alluvial plot
# ---------------------------
t1 <- trans_abund$new(dataset = dataset, taxrank = opts$taxrank, ntaxa = opts$ntaxa)

p <- t1$plot_bar(
  bar_full = FALSE,
  use_alluvium = TRUE,
  clustering = TRUE,
  xtext_angle = opts$xangle,
  xtext_size = opts$xsize,
  color_values = RColorBrewer::brewer.pal(opts$ntaxa, "Set2")
)

# ---------------------------
# 6. Save plot
# ---------------------------
ggsave(paste0(opts$output, "alluvial_microeco_plot.pdf"), plot = p, width = opts$width, height = opts$height, units = "mm")

message("Alluvial plot saved to: ", opts$output)
