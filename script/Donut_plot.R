#!/usr/bin/env Rscript

options(warn = -1)

site <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN"

ensure_package <- function(package_name){
  if(!requireNamespace(package_name, quietly = TRUE)){
    install.packages(package_name, repos = site)
  }
  library(package_name, character.only = TRUE)
}

ensure_package("optparse")
ensure_package("microeco")
ensure_package("dplyr")
ensure_package("ggplot2")
ensure_package("tibble")
ensure_package("ggrepel")
ensure_package("scales")
ensure_package("RColorBrewer")

# ---- Command Line Options ----
option_list <- list(
  make_option(c("-i", "--otu_table"), type = "character", default = "result2/tax/otutab2.txt"),
  make_option(c("-d", "--metadata"), type = "character", default = "result2/tax/metadata.txt"),
  make_option(c("-t", "--taxonomy"), type = "character", default = "result2/tax/taxonomy.txt"),
  make_option(c("-o", "--output_dir"), type = "character", default = "result2/tax"),
  make_option(c("-g", "--group"), type = "character", default = "Group"),
  make_option(c("-n", "--ntaxa"), type = "integer", default = 8),
  make_option(c("-l", "--label_threshold"), type = "double", default = 1)
)
opts <- parse_args(OptionParser(option_list = option_list))

# ---- Read Input Data ----
otu_table <- read.table(opts$otu_table, header = TRUE, row.names = 1, sep = "\t",
                        check.names = FALSE, comment.char = "")
metadata <- read.table(opts$metadata, header = TRUE, row.names = 1, sep = "\t",
                       check.names = FALSE, comment.char = "")
taxonomy <- read.table(opts$taxonomy, header = FALSE, row.names = 1, sep = "\t",
                       check.names = FALSE, comment.char = "")

if ("OTUID" %in% rownames(taxonomy)) {
  taxonomy <- taxonomy[rownames(taxonomy) != "OTUID", ]
}
colnames(taxonomy) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# ---- Create microeco Object ----
mt <- microtable$new(
  otu = otu_table,
  tax = taxonomy,
  sample = metadata
)

# ---- Initialize trans_abund Object ----
t1 <- trans_abund$new(
  dataset = mt,
  taxrank = "Phylum",
  ntaxa = opts$ntaxa,
  groupmean = opts$group
)

if(is.null(t1$data_abund) || nrow(t1$data_abund) == 0){
  message("Abundance data not available. Verify your input and parameters.")
  quit(status = 1)
}

# ---- Use group info already in t1$data_abund ----
data_abund_with_group <- t1$data_abund %>%
  rename(Taxon = Taxonomy) %>%
  rename(Group = Sample)  # Sample is actually group here

group_names <- unique(data_abund_with_group$Group)

cat("Detected group names:\n")
print(group_names)

# ---- Unique color palette ----
all_taxa <- unique(t1$data_abund$Taxonomy)
set.seed(42)
custom_colors <- grDevices::rainbow(length(all_taxa), s = 0.7, v = 0.9)
names(custom_colors) <- all_taxa

# ---- Donut plot function for a single group ----
create_single_group_donut <- function(data_abund_with_group, group_name, color_palette, ntaxa = 8, label_threshold = 1) {
  data_group <- data_abund_with_group %>%
    filter(Group == group_name) %>%
    group_by(Taxon) %>%
    summarize(MeanAbundance = mean(Abundance), .groups = "drop") %>%
    arrange(desc(MeanAbundance))
  
  if (nrow(data_group) == 0) {
    warning(paste("No data for group", group_name))
    return(NULL)
  }
  
  # Keep top ntaxa taxa, lump others as "Others"
  if (nrow(data_group) > ntaxa) {
    top_taxa <- data_group$Taxon[1:ntaxa]
    data_group$Taxon <- ifelse(data_group$Taxon %in% top_taxa, data_group$Taxon, "Others")
    data_group <- data_group %>%
      group_by(Taxon) %>%
      summarize(MeanAbundance = sum(MeanAbundance), .groups = "drop") %>%
      arrange(desc(MeanAbundance))
  }
  
  data_group <- data_group %>%
    mutate(Percentage = MeanAbundance / sum(MeanAbundance) * 100) %>%
    mutate(cumulative = cumsum(Percentage) - Percentage / 2)
  
  data_group$label <- ifelse(data_group$Percentage >= label_threshold, paste0(round(data_group$Percentage, 1), "%"), "")
  
  # Assign colors for taxa present in this group
  color_use <- color_palette[names(color_palette) %in% data_group$Taxon]
  if ("Others" %in% data_group$Taxon) {
    color_use["Others"] <- "grey80"
  }
  
  ggplot(data_group, aes(x = 2, y = MeanAbundance, fill = Taxon)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    theme_minimal(base_size = 16) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    scale_fill_manual(values = color_use, drop = FALSE) +
    ggrepel::geom_label_repel(
      aes(x = 2.7, y = cumulative, label = label),
      box.padding = 0.3,
      point.padding = 0.1,
      segment.color = "grey50",
      show.legend = FALSE,
      size = 5,
      na.rm = TRUE
    ) +
    geom_text(aes(x = 0, y = sum(MeanAbundance)/2, label = group_name),
              size = 7, fontface = "bold", inherit.aes = FALSE) +
    xlim(c(0, 3)) +
    labs(title = NULL, fill = "Phylum")
}

# ---- Create output directory if not exists ----
if (!dir.exists(opts$output_dir)) dir.create(opts$output_dir, recursive = TRUE)

# ---- Save one donut plot per group, each as separate PDF file ----
for (g in group_names) {
  plot_g <- create_single_group_donut(
    data_abund_with_group,
    g,
    custom_colors,
    ntaxa = opts$ntaxa,
    label_threshold = opts$label_threshold
  )
  
  if (!is.null(plot_g)) {
    filename <- file.path(opts$output_dir, paste0("donut_plot_", g, ".pdf"))
    ggsave(filename = filename, plot = plot_g, width = 8, height = 8)
    message("Saved donut plot for group ", g, " to ", filename)
  } else {
    warning("No plot generated for group ", g)
  }
}

# ---- Combined donut plot (all samples) ----
create_combined_donut <- function(t1_object, color_palette, label_threshold = 1) {
  data_filtered <- t1_object$data_abund %>%
    rename(Taxon = Taxonomy) %>%
    group_by(Taxon) %>%
    summarize(MeanAbundance = mean(Abundance), .groups = "drop") %>%
    mutate(Percentage = MeanAbundance / sum(MeanAbundance) * 100) %>%
    arrange(desc(MeanAbundance)) %>%
    mutate(cumulative = cumsum(Percentage) - Percentage / 2)
  
  data_filtered$label <- ifelse(data_filtered$Percentage >= label_threshold, paste0(round(data_filtered$Percentage, 1), "%"), "")
  
  color_use <- color_palette[names(color_palette) %in% data_filtered$Taxon]
  
  ggplot(data_filtered, aes(x = 2, y = MeanAbundance, fill = Taxon)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    theme_minimal(base_size = 16) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    scale_fill_manual(values = color_use, drop = FALSE) +
    ggrepel::geom_label_repel(
      aes(x = 2.7, y = cumulative, label = label),
      box.padding = 0.3,
      point.padding = 0.1,
      segment.color = "grey50",
      show.legend = FALSE,
      size = 5,
      na.rm = TRUE
    ) +
    xlim(c(0.5, 3)) +
    labs(title = "Donut Plot of All Samples", fill = "Phylum")
}

combined_donut <- create_combined_donut(t1, custom_colors, label_threshold = opts$label_threshold)
ggsave(filename = file.path(opts$output_dir, "donut_plot_all_samples.pdf"),
       plot = combined_donut, width = 8, height = 8)

# ---- Combined donut plot by group in one multi-page PDF ----
pdf(file.path(opts$output_dir, "donut_plot_by_group.pdf"), width = 8, height = 8)
for (g in group_names) {
  plot_g <- create_single_group_donut(
    data_abund_with_group,
    g,
    custom_colors,
    ntaxa = opts$ntaxa,
    label_threshold = opts$label_threshold
  )
  if (!is.null(plot_g)) {
    print(plot_g)  # Important: print ggplot objects inside pdf device
  }
}
dev.off()

# ---- Radar plot ----
pdf(file.path(opts$output_dir, "radar_plot_micro.pdf"), width = 10, height = 8)
t1$plot_radar(
  color_values = custom_colors,
  values.radar = c("0%", "25%", "50%"),
  grid.min = 0,
  grid.mid = 0.25,
  grid.max = 0.5
)
dev.off()

cat("All plots saved to:", opts$output_dir, "\n")
