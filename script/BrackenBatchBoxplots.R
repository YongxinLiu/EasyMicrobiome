
#!/usr/bin/env Rscript

# ============================================================
# plot_species_boxplot.R
# ============================================================

# ===========================
# 1. 定义需要的包
# ===========================
packages <- c(
  "optparse",
  "data.table",
  "dplyr",
  "ggplot2",
  "ggpubr",
  "cowplot"
)

# ===========================
# 2. 检查是否安装，没有就安装
# ===========================
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste0("⚠️ Installing package: ", pkg))
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
}

# ===========================
# 3. 加载包
# ===========================
suppressPackageStartupMessages({
  lapply(packages, library, character.only = TRUE)
})


# ======================
# 1. 命令行参数
# ======================
option_list <- list(
  make_option(c("-e", "--exp"), type = "character",
              help = "Abundance table (e.g. bracken.S.txt)", metavar = "file"),
  make_option(c("-m", "--meta"), type = "character",
              help = "Metadata table", metavar = "file"),
  make_option(c("-s", "--species"), type = "character",
              help = "Species list, comma-separated (e.g. \"A,B,C\")"),
  make_option(c("--groups"), type = "character",
              default = "VfNnVa5Rs,ViNnVa5Rs,VfNnVb5Rs,ViNnVb5Rs",
              help = "Group order, comma-separated"),
  make_option(c("-o", "--out"), type = "character",
              default = "Species_boxplot.pdf",
              help = "Output PDF [default: %default]"),
  make_option(c("--sample_col"), type = "character",
              default = "SampleID",
              help = "Sample column name in metadata [default: %default]"),
  make_option(c("--group_col"), type = "character",
              default = "GroupID",
              help = "Group column name in metadata [default: %default]"),
  make_option(c("--ncol"), type = "integer",
              default = 4,
              help = "Number of columns in plot layout [default: %default]"),
  make_option(c("--pdf_width"), type = "double", default = 14,
              help = "PDF width [default: %default], unit = cm"),
  make_option(c("--pdf_height"), type = "double", default = 10,
              help = "PDF height [default: %default], unit = cm")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$exp) || is.null(opt$meta) || is.null(opt$species)) {
  stop("❌ Required arguments: -e, -m, --species\n", call. = FALSE)
}

# ======================
# 2. 解析 species（命令行输入）
# ======================
species <- unlist(strsplit(opt$species, ","))
species <- trimws(species)
if (length(species) == 0) stop("❌ No species provided.", call. = FALSE)

# ======================
# 3. 解析 groups（命令行输入）
# ======================
group_levels <- unlist(strsplit(opt$groups, ","))
group_levels <- trimws(group_levels)

# ======================
# 4. 读入丰度表
# ======================
Exp <- fread(opt$exp, sep = "\t", fill = TRUE, data.table = FALSE)
if (!"Taxonomy" %in% colnames(Exp)) {
  stop("❌ Abundance table must contain 'Taxonomy' column.", call. = FALSE)
}
rownames(Exp) <- Exp$Taxonomy
Exp <- Exp[, -1]

# ======================
# 5. 每列相对丰度
# ======================
Exp <- Exp %>%
  mutate(across(where(is.numeric), ~ .x / sum(.x)))

Exp2 <- as.data.frame(t(Exp))
Exp2$SampleID <- rownames(Exp2)
colnames(Exp2) <- gsub("\\.", " ", colnames(Exp2))
Exp2 <- Exp2[, c(ncol(Exp2), 1:(ncol(Exp2) - 1))]

# ======================
# 6. 物种存在性检查
# ======================
missing_sp <- setdiff(species, colnames(Exp2))
if (length(missing_sp) > 0) {
  warning("⚠️ Species not found:\n", paste(missing_sp, collapse = ", "))
}
species <- intersect(species, colnames(Exp2))
if (length(species) == 0) {
  stop("❌ None of the species found in abundance table.", call. = FALSE)
}
Exp_plot <- Exp2[, c("SampleID", species)]

# ======================
# 7. 读入 metadata 并合并
# ======================
meta <- read.table(opt$meta, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (!opt$sample_col %in% colnames(meta) || !opt$group_col %in% colnames(meta)) {
  stop("❌ sample_col or group_col not found in metadata.", call. = FALSE)
}

Exp_plot <- Exp_plot[Exp_plot$SampleID %in% meta[[opt$sample_col]], ]
Exp_plot <- Exp_plot[match(meta[[opt$sample_col]], Exp_plot$SampleID), ]

Exp_plot$group <- meta[[opt$group_col]]
Exp_plot$group <- factor(Exp_plot$group, levels = group_levels)
Exp_plot <- Exp_plot[!is.na(Exp_plot$group), ]

# ======================
# 8. 作图
# ======================
colors <- c("#EB746A", "#7AA82C", "#1EB5B8", "#A07DB7")
comparisons <- combn(group_levels, 2, simplify = FALSE)

plist <- list()
for (i in seq_along(species)) {
  bar_tmp <- Exp_plot[, c(species[i], "group")]
  colnames(bar_tmp) <- c("Relative abundance", "group")
  
  p <- ggboxplot(
    bar_tmp,
    x = "group",
    y = "Relative abundance",
    color = "group",
    add = "jitter",
    palette = colors
  ) +
    ggtitle(species[i]) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    ) +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = comparisons,
      label = "p.signif"
    )
  
  plist[[i]] <- p
}

# ======================
# 9. 输出 PDF
# ======================
#pdf(opt$out, width = 14, height = 10)
pdf(opt$out, width = opt$pdf_width, height = opt$pdf_height)
plot_grid(plotlist = plist, ncol = opt$ncol)
dev.off()

cat("✅ Done! Output:", opt$out, "\n")
