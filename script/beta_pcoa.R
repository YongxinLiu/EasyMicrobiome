#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(vegan)
  library(ggplot2)
})

# 设置命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "OTU count 表格文件（行为OTU，列为样本）"),
  make_option(c("-m", "--metadata"), type = "character", help = "样本分组文件（包含Group列）"),
  make_option(c("-g", "--group"), type = "character", default = "Group", help = "元数据中用于分组的列名 [默认: %default]"),
  make_option(c("-o", "--output"), type = "character", default = "PCoA.pdf", help = "输出PDF文件路径")
)

# 解析命令行
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 检查必要参数
if (is.null(opt$input) || is.null(opt$metadata)) {
  print_help(opt_parser)
  stop("请提供输入的 OTU 表文件和元数据文件", call. = FALSE)
}

# 读取 OTU 表
contigs <- read.delim(opt$input, row.names = 1, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
contigs <- data.frame(t(contigs))

# 预处理（log10转化并移除稀疏物种）
contigs <- log10(contigs + 1)
contigs <- contigs[, colSums(contigs > 0) > (0.1 * nrow(contigs))]

# Bray-Curtis 距离
distance <- vegdist(contigs, method = "bray")
distance_matrix <- as.matrix(distance)
distance_matrix[is.na(distance_matrix)] <- 0

# 主坐标分析（PCoA）
pcoa <- cmdscale(distance_matrix, k = (nrow(contigs) - 1), eig = TRUE)
plot_data <- data.frame(pcoa$points)[, 1:2]
names(plot_data) <- c("PCoA1", "PCoA2")

# 计算百分比
eig <- pcoa$eig
PCOA1 <- format(100 * eig[1] / sum(eig), digits = 4)
PCOA2 <- format(100 * eig[2] / sum(eig), digits = 4)

# 读取分组文件
group <- read.delim(opt$metadata, row.names = 1, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# 合并分组信息
data <- merge(plot_data, group, by = "row.names", all = TRUE)
names(data)[1] <- "Sample"
group_col <- opt$group

# 清除空分组
data <- data[!is.na(data[[group_col]]) & data[[group_col]] != "", ]

# 设置默认颜色和形状（根据分组自动生成）
group_levels <- unique(data[[group_col]])
n_groups <- length(group_levels)
default_colors <- c("#1EB5B8", "#edae11","#A07DB7","#7AA82C", "#547aa5", "#916BBF", "#46ACC8")[1:n_groups]
default_shapes <- c(6, 21, 23, 24, 22, 25)[1:n_groups]

# 画图
p <- ggplot(data, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes_string(color = group_col, shape = group_col), size = 4) +
  scale_color_manual(values = default_colors) +
  scale_shape_manual(values = default_shapes) +
  stat_ellipse(aes_string(color = group_col), linetype = "solid", level = 0.98, size = 0.85) +
  labs(
    x = paste0("PCoA1 (", PCOA1, "%)"),
    y = paste0("PCoA2 (", PCOA2, "%)"),
    title = "PCoA"
  ) +
  #theme_bw() +
  theme_classic()+
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.background = element_rect(fill = "transparent"),
    legend.position = c(0.1, 0.1)
  ) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted")

# 保存图像
ggsave(opt$output, plot = p, device = "pdf", width = 6, height = 5.5)

# PerMANOVA 置换多元方差分析
# 距离指数 Distance matrix
dis = distance_matrix
sub_design = group
idx = rownames(sub_design) %in% rownames(dis)
sub_design = sub_design[idx,]
sub_dis = dis[rownames(sub_design),rownames(sub_design)]
dis1 <- as.dist(sub_dis)

# anonis
adonis_result <- (adonis2(dis1~Group, data = sub_design, permutations = 999))
p_val <- adonis_result$`Pr(>F)`[1]
print(paste0("PERMANOVA 分析的 p 值为 ", signif(p_val, 3),
             ifelse(p_val < 0.05, "，结果具有统计学显著性。", "，结果无统计学显著性。")))

