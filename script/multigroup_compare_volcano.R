#!/usr/bin/env Rscript

# Description:
# Perform DESeq2 pairwise comparisons for all group pairs in metadata,
# then generate a combined faceted volcano plot colored by significance status,
# with enhanced colors and gradient background.
#
# Usage example:
# Rscript CompareAllPairs_SignificanceColored.R --input result/otutab.txt --metadata result/metadata.txt --group_column Group --output result/compare_all_pairs/

options(warn = -1)  # Suppress warnings globally

# -------- 1. Load Required Packages and Install if Missing --------
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

install_if_missing("DESeq2")
install_if_missing("optparse")
install_if_missing("ggplot2")
install_if_missing("dplyr")
install_if_missing("tidyr")
install_if_missing("scales")

library(DESeq2)
library(optparse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# -------- 2. Parse Command-line Arguments --------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "result2/otutab.txt",
              help = "Path to OTU count table (tab-delimited, samples as columns) [default %default]"),
  make_option(c("-d", "--metadata"), type = "character", default = "result2/metadata.txt",
              help = "Path to metadata file (tab-delimited, samples as rows) [default %default]"),
  make_option(c("-o", "--output"), type = "character", default = "result_illumina/tax/",
              help = "Output directory for results [default %default]"),
  make_option(c("-g", "--group_column"), type = "character", default = "Group",
              help = "Column name in metadata to use as grouping factor [default %default]"),
  make_option(c("-t", "--threshold"), type = "numeric", default = 0.1,
              help = "Log2 fold change threshold for significance [default %default]"),
  make_option(c("-p", "--pvalue"), type = "numeric", default = 0.05,
              help = "Adjusted p-value cutoff for significance [default %default]")
)

opts <- parse_args(OptionParser(option_list = option_list))

# -------- 3. Validate Input Files and Output Directory --------
# if (!file.exists(opts$input)) {
#   stop("Input OTU table file not found: ", opts$input)
# }
# if (!file.exists(opts$metadata)) {
#   stop("Metadata file not found: ", opts$metadata)
# }

# if (!dir.exists(opts$output)) {
#   dir.create(opts$output, recursive = TRUE)
#   message("Created output directory: ", opts$output)
# }

# -------- 4. Load Data --------
otu <- read.table(opts$input, header = TRUE, row.names = 1, sep = "\t",
                  check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
otu[] <- lapply(otu, function(x) as.numeric(as.character(x)))
if (any(is.na(otu))) {
  stop("Non-numeric values detected in OTU table after coercion. Please check your input file.")
}

metadata <- read.table(opts$metadata, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)

if (!(opts$group_column %in% colnames(metadata))) {
  stop("Group column '", opts$group_column, "' not found in metadata.")
}

# -------- 5. Match Samples --------
common_samples <- intersect(colnames(otu), rownames(metadata))
if (length(common_samples) == 0) {
  stop("No common samples found between OTU table and metadata. Please check sample names.")
}

otu <- otu[, common_samples]
metadata <- metadata[common_samples, , drop = FALSE]

# -------- 6. DESeq2 Analysis --------
dds <- DESeqDataSetFromMatrix(countData = round(otu),
                              colData = metadata,
                              design = as.formula(paste0("~", opts$group_column)))

dds <- estimateSizeFactors(dds, type = "poscounts")

dds <- tryCatch({
  DESeq(dds, fitType = "local")
}, error = function(e) {
  stop("Error running DESeq2: ", e$message)
})

# -------- 7. Get All Pairwise Group Comparisons --------
groups <- unique(metadata[[opts$group_column]])
pairwise_comparisons <- combn(groups, 2, simplify = FALSE)

all_results <- list()

for (cmp in pairwise_comparisons) {
  message("Processing comparison: ", cmp[1], " vs ", cmp[2])
  contrast <- c(opts$group_column, cmp[1], cmp[2])
  
  res <- tryCatch({
    results(dds, contrast = contrast)
  }, error = function(e) {
    warning("Skipping comparison ", cmp[1], " vs ", cmp[2], ": ", e$message)
    return(NULL)
  })
  
  if (is.null(res)) next
  
  res_df <- as.data.frame(res)
  res_df$OTU <- rownames(res_df)
  res_df$comparison <- paste(cmp[1], "vs", cmp[2])
  
  # Filter and assign status
  res_df <- res_df %>%
    filter(!is.na(padj)) %>%
    mutate(status = case_when(
      padj < opts$pvalue & log2FoldChange > opts$threshold ~ "Sigup",
      padj < opts$pvalue & log2FoldChange < -opts$threshold ~ "Sigdown",
      TRUE ~ "NotSig"
    ))
  
  all_results[[paste(cmp, collapse = "_vs_")]] <- res_df
}

if (length(all_results) == 0) {
  stop("No valid pairwise comparisons could be performed.")
}

combined_res <- bind_rows(all_results)

# -------- 8. Save Combined Results --------
combined_results_file <- file.path(opts$output, "all_pairwise_results_combined.txt")
write.table(combined_res, file = combined_results_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Combined pairwise results saved to: ", combined_results_file)

# -------- 9. Enhanced Volcano Plot with Gradient Background and Significance Colors --------
# Calculate plot ranges for gradient background
xrange <- range(combined_res$log2FoldChange, na.rm = TRUE)
yrange <- range(-log10(combined_res$padj), na.rm = TRUE)

# Create data frame for vertical gradient rectangles
bg_df <- data.frame(
  xmin = seq(xrange[1], xrange[2], length.out = 100)[-100],
  xmax = seq(xrange[1], xrange[2], length.out = 100)[-1],
  ymin = yrange[1],
  ymax = yrange[2]
)
# Assign gradient fill colors from light blue to white left to right
bg_df$fill <- colorRampPalette(c("#cceeff", "white"))(nrow(bg_df))

# Define colors for significance status
status_colors <- c(Sigup = "red", Sigdown = "blue", NotSig = "grey70")

p <- ggplot(combined_res, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_rect(data = bg_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            inherit.aes = FALSE, alpha = 0.3) +
  scale_fill_identity() +  # use fill colors as is
  geom_point(alpha = 0.7, size = 1.8) +
  facet_wrap(~ comparison, scales = "free") +
  scale_color_manual(values = status_colors, name = "Significance") +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "lightblue", color = NA),
    strip.text = element_text(color = "black", face = "bold"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  labs(x = "Log2 Fold Change", y = "-log10 Adjusted P-value") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

volcano_plot_file <- file.path(opts$output, "all_pairwise_volcano_plots_significance_colored.pdf")
ggsave(volcano_plot_file, p, width = 14, height = 7)
message("Enhanced volcano plot with significance colors saved to: ", volcano_plot_file)

#ggsave(paste0(opts$output, "all_pairwise_volcano_plots_significance_colored.pdf"), plot = p, width = 14, height = 7, units = "mm")

#message("all_pairwise_volcano_plots saved to: ", opts$output)

#cat("✅ DESeq2 all pairwise comparisons completed successfully.\nResults and plots saved to:", opts$output, "\n")
