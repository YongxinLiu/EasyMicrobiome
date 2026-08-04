#!/usr/bin/env Rscript

# Copyright 2024 De-feng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, et al. 2021. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12: 315-330. https://doi.org/10.1007/s13238-020-00724-8

# 手动运行脚本，使用 Ctrl+Shift+H 或 Session 需要设置工作目录
# Set Work Directory - Choose Directory / To Source File Location

# 1.1 简介 Introduction #----
 
# 程序功能：处理得到的物种组成数据作为SparCC计算的输入
# Functions: Processing Metaphlan4 species abundance table
# as SparCC input files for each group

# Clean environment
rm(list = ls())

options(warn = -1)

# =========================================================
# Parameters
# =========================================================

site = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"

# Install optparse if needed
if (!suppressWarnings(
  suppressMessages(
    require("optparse",
            character.only = TRUE,
            quietly = TRUE,
            warn.conflicts = FALSE)
  )
)) {
  
  install.packages("optparse", repos = site)
  
  require(optparse)
}

# Parse arguments
option_list = list(
  
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "result12/metaphlan4/Species.txt",
    help = "Metaphlan4 species table"
  ),
  
  make_option(
    c("-g", "--group"),
    type = "character",
    default = "result12/metadata.txt",
    help = "Metadata file"
  ),
  
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "result12/metaphlan4/",
    help = "Output directory"
  ),
  
  make_option(
    c("--groups"),
    type = "character",
    default = NULL,
    help = "Groups to extract, separated by comma. Example: Cancer,Normal"
  )
)

opts = parse_args(OptionParser(option_list = option_list))

cat("=====================================\n")
cat("Parameters:\n")
print(opts)
cat("=====================================\n")

# =========================================================
# Load packages
# =========================================================

packages <- c(
  "reshape2",
  "ggplot2",
  "ggprism",
  "dplyr",
  "plyr",
  "igraph"
)

for (p in packages) {
  
  suppressWarnings(
    suppressMessages(
      library(p,
              character.only = TRUE)
    )
  )
}

# =========================================================
# Create output directory
# =========================================================

if (!dir.exists(opts$output)) {
  
  dir.create(opts$output,
             recursive = TRUE)
}

# =========================================================
# Read metadata
# =========================================================

design <- read.table(
  opts$group,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# Check Group column
if (!"Group" %in% colnames(design)) {
  
  stop("ERROR: metadata file must contain a column named 'Group'")
}

# =========================================================
# Read species abundance table
# =========================================================

df3 <- read.table(
  opts$input,
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

df_species <- df3

# =========================================================
# Aggregate taxonomy
# =========================================================

# 获取第一列列名
feature_col <- colnames(df_species)[1]

cat("Feature column:", feature_col, "\n")

data <- aggregate(
  df_species[, -1],
  by = list(df_species[[feature_col]]),
  FUN = sum
)

colnames(data)[1] <- feature_col

rownames(data) <- data[[feature_col]]

data <- data[, -1]

# =========================================================
# Convert abundance
# =========================================================

data3 <- apply(data,
               2,
               function(x) x / 100)

cat("Data dimension:\n")
print(dim(data3))

# SparCC recommended scaling
# data3 <- data3 * 100000

# =========================================================
# Format OTU table
# =========================================================

# keep numeric matrix
data3 <- as.data.frame(data3)

# species as rownames
rownames(data3) <- rownames(data)

# samples as columns
colnames(data3) <- colnames(data)

# transpose:
# rows = species
# cols = samples
otutab <- as.data.frame(data3)

# ensure numeric
otutab[] <- lapply(otutab, as.numeric)

# check
cat("OTU table dimension:\n")
print(dim(otutab))

cat("OTU table preview:\n")
print(head(otutab[,1:min(5,ncol(otutab))]))


# =========================================================
# Determine groups
# =========================================================

if (is.null(opts$groups)) {
  
  group_list <- unique(design$Group)
  
} else {
  
  group_list <- unlist(
    strsplit(opts$groups, ",")
  )
}

cat("Groups to process:\n")
print(group_list)

# =========================================================
# Function for extracting group
# =========================================================

extract_group <- function(group_name,
                          design,
                          otutab,
                          output_dir) {
  
  cat("\nProcessing group:", group_name, "\n")
  
  # subset metadata
  sub_design <- subset(
    design,
    Group %in% c(group_name)
  )
  
  sub_design$Group <- factor(
    sub_design$Group,
    levels = c(group_name)
  )
  
  # sample intersection
  idx <- rownames(sub_design) %in% colnames(otutab)
  
  sub_design <- sub_design[idx, ]
  
  # subset otu table
  sub_otutab <- otutab[, rownames(sub_design)]
  
  sub_otutab <- as.data.frame(sub_otutab)
  
  # output file
  outfile <- paste0(
    output_dir,
    "/",
    group_name,
    "_sparcc.txt"
  )
  
  write.table(
    sub_otutab,
    file = outfile,
    row.names = TRUE,
    sep = "\t",
    quote = FALSE,
    col.names = TRUE
  )
  
  cat("Output:", outfile, "\n")
  cat("Samples:", ncol(sub_otutab), "\n")
  cat("Species:", nrow(sub_otutab), "\n")
}

# =========================================================
# Run all groups
# =========================================================

for (g in group_list) {
  
  extract_group(
    group_name = g,
    design = design,
    otutab = otutab,
    output_dir = opts$output
  )
}

cat("\nAll done!\n")


