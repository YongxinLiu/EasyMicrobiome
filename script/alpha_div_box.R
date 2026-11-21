# 加载必要包
library(ggplot2)
library(dplyr)
library(ggsignif)
library(scales)
library(optparse)

# 定义命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Alpha diversity table file (e.g., alpha.txt)"),
  make_option(c("-m", "--metadata"), type = "character", help = "Metadata file (e.g., metadata.txt)"),
  make_option(c("-a", "--alpha_index"), type = "character", help = "Comma-separated alpha diversity indices (e.g., shannon_e,richness)", default = "shannon_e,richness"),
  make_option(c("-g", "--group"), type = "character", help = "Column name in metadata used for grouping (e.g., Group)", default = "Group"),
  make_option(c("-o", "--out_prefix"), type = "character", help = "Output file prefix", default = "alpha")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 读取数据
diversity_data <- read.table(opt$input, header = TRUE, sep = "\t")
metadata <- read.table(opt$metadata, header = TRUE, sep = "\t")

# 合并数据
merged_data <- merge(diversity_data, metadata, by.x = "SampleID", by.y = "SampleID")

# 解析 alpha 多样性指标
alpha_indices <- unlist(strsplit(opt$alpha_index, ","))

# 分组列名
group_col <- opt$group

# 获取唯一的组名组合进行比较（最多三组默认处理）
group_levels <- unique(merged_data[[group_col]])
if (length(group_levels) < 2) {
  stop("分组列中组数不足，至少需要两个组")
}
comparisons <- combn(group_levels, 2, simplify = FALSE)

# 保存 p 值结果
p_values <- data.frame(Comparison = sapply(comparisons, paste, collapse = " vs "))

# 颜色设置（如组数超过3可扩展）
#color_palette <- c("#edae11", "#0f6657", "#c74732", "#1EB5B8", "#7AA82C", "#A07DB7")
color_palette <- c("#c74732", "#1EB5B8", "#A07DB7", "#edae11", "#0f6657",  "#7AA82C")

# 依次处理各个多样性指数
for (alpha in alpha_indices) {
  # 计算 t 检验 p 值
  #pval_vec <- sapply(comparisons, function(x) {
  #  subset_data <- subset(merged_data, merged_data[[group_col]] %in% x)
  #  t.test(as.formula(paste(alpha, "~", group_col)), data = subset_data)$p.value
  #})
  #p_values[[paste0(alpha, "_pval")]] <- pval_vec
  
  # 计算显著性标记位置
  ymax <- max(merged_data[[alpha]], na.rm = TRUE)
  offset <- (ymax - min(merged_data[[alpha]], na.rm = TRUE)) * 0.05
  y_positions <- seq(ymax + offset, by = offset * 1.5, length.out = length(comparisons))
  
  # 绘图
  p <- ggplot(merged_data, aes_string(x = group_col, y = alpha, fill = group_col)) +
    geom_boxplot(position = position_dodge(width = 0.4),
                 outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 2, aes_string(color = group_col), shape = 16, alpha = 0.7) +
    labs(y = alpha, x = "") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold")) +
    scale_fill_manual(values = color_palette[seq_along(group_levels)]) +
    scale_color_manual(values = color_palette[seq_along(group_levels)]) +
    geom_signif(comparisons = comparisons,
                map_signif_level = TRUE,
                #test = "t.test",
                test = "wilcox.test",
                tip_length = 0.03,
                y_position = y_positions) +
    scale_y_continuous(breaks = pretty_breaks(n = 5),
                       expand = expansion(mult = c(0.05, 0.2)))
  
  # 保存图像
  ggsave(paste0(opt$out_prefix, "_", alpha, ".pdf"), plot = p, width = 3, height = 4)
}

