#!/usr/bin/env Rscript
# 关闭警告信息
options(warn = -1)

# 设置清华源
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"

# 解析命令行参数
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE)))) {
  install.packages("optparse", repos=site)
  require("optparse", character.only = TRUE)
}

# 设置命令行参数
option_list <- list(
  #make_option(c("-i", "--input"), type = "character", default = "./result2/compare/saliva-plaque2.txt", help = "组间比较结果表格"),
  make_option(c("-i", "--input"), type="character", default="result2/compare/saliva_plaque2.txt",
              help="Taxonomy composition [default %default]"),
  make_option(c("-g", "--group"), type = "character", default = "saliva-plaque", 
              help = "Column name in metadata used for grouping (e.g., Group)"),
  make_option(c("-o", "--output"), type="character", default="",
              help="Output directory [default %default]"),
  #make_option(c("-o", "--output"), type = "character", default = "./result2/compare/", help = "Valcano_plot.pdf输出PDF文件路径"),
  make_option(c("-w", "--width"), type="numeric", default=120 * 1.5,
              help="Figure width (mm) [default %default]"),
  make_option(c("-e", "--height"), type="numeric", default=80 * 1.5,
              help="Figure height (mm) [default %default]")
)

# 解析命令行
opts = parse_args(OptionParser(option_list=option_list))

if (opts$output == "") { opts$output = opts$input }

suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(ggplot2)))

group <- opts$group

# 读取数据
df <- read.table(opts$input, header = TRUE, sep = "\t")

#df <- read.table("Figure2/data/saliva-plaque2.txt", header = TRUE, sep = "\t")
#df <- read.table("Figure2/data/feces-plaque.txt", header = TRUE, sep = "\t")
#df <- read.table("Figure2/data/feces-saliva.txt", header = TRUE, sep = "\t")


# 处理P值（去掉NA值）
df$PValue <- as.numeric(df$PValue)
df$FDR <- as.numeric(df$FDR)
df$log2FC <- as.numeric(df$log2FC)
df <- df[!is.na(df$PValue), ]

# 计算 -log10(P值)
df$neg_log10_pval <- -log10(df$PValue)

# 设置显著性类别
df$Significance <- "Not Significant"
df$Significance[df$log2FC > 1 & df$FDR < 0.5] <- "Enriched"
df$Significance[df$log2FC < -1 & df$FDR < 0.5] <- "Depleted"

# 选择要标注的 ASV（仅标注 Enriched 和 Depleted）
df$Label <- ifelse(df$Significance %in% c("Enriched", "Depleted"), df$ID, NA)

# 绘制火山图
p1 <- ggplot(df, aes(x = log2FC, y = neg_log10_pval, color = Significance)) +
  geom_point(alpha = 0.8, size = 2) +  # 散点透明度和大小
  scale_color_manual(values = c("Enriched" = "#c74732", "Depleted" = "#1EB5B8", "Not Significant" = "gray")) +  # 颜色
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "black") +  # 显著性阈值线
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # 折叠变化阈值
  geom_text_repel(aes(label = Label), size = 2, box.padding = 0.3, max.overlaps = 15) +  # 防止标注重叠
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray", linetype = "dotted"),  # 横纵虚线
    panel.grid.minor = element_blank(),  # 去除次要网格线
    axis.line = element_blank(),  # 去除x和y轴线
    axis.text = element_text(color = "black"),  # x、y轴刻度文字颜色
    axis.ticks = element_blank(),  # 去除x和y轴刻度线
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # 外四周的实线
  ) +
  scale_x_continuous(limits = c(-8,8))+
  labs(title = "Volcano Plot of ASVs",
       x = "log2 Fold Change",
       y = "-log10(PValue)",
       color = "Significance")
#p1
#ggplot2::ggsave(paste("Figure2/data/Two_group_volcano_plot_amplicon3",".pdf", sep=""), 
#                p1, width=120 * 1.5, height=80 * 1.5, unit='mm')
#ggsave(paste0(opts$output, ".group.pdf"), p1, width = opts$width, height = opts$height, units = "mm")
# 保存图像
#ggsave(opts$output, plot = p1, device = "pdf", width = opts$width, height = opts$height, units = "mm")
#ggsave(opts$output, plot = p1, device = "pdf", width = opts$width, height = opts$height, units = "mm")
ggsave(paste0(opts$output, group, "_", "volcano.pdf"), plot = p1, width = opts$width, height = opts$height,units = "mm", dpi = 600)
