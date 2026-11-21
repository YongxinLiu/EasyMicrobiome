#!/usr/bin/env Rscript

# Copyright 2024 De-feng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, et al. 2021. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12: 315-330. https://doi.org/10.1007/s13238-020-00724-8

# 手动运行脚本，使用 Ctrl+Shift+H 或 Session 需要设置工作目录
# Set Work Directory - Choose Directory / To Source File Location

# Clean enviroment object
rm(list=ls()) 

# 1.1 简介 Introduction #----

# 程序功能：MAGs完整性和污染率关系图
# Functions: The relationship between MAGs completeness and contamination


options(warn = -1) # Turn off warning

# 1.2 参数 Parameters #----
# 设置清华源加速下载
# (Optional) Set up Tsinghua Mirror to speed up download
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析包是否安装，没安装则安装，然后加载
# Determine whether the command line parsing package is installed, install it if it is not installed, then load
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}

# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="result/checkm2/taxonomy_merge2.txt",
                help="MAGs taxonomy tables [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/checkm2/",
                help="Output for relationship between MAGs completeness and contamination [default %default]") 
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
print("You are using the following parameters:")
print(opts)



# Install related packages
if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggplot2","RColorBrewer","dplyr","gridExtra","cowplot")) # ,"vegan"
}
# load related packages
suppressWarnings(suppressMessages(library("ggplot2")))
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("RColorBrewer")))
# suppressWarnings(suppressMessages(library("vegan")))
suppressWarnings(suppressMessages(library("gridExtra")))
suppressWarnings(suppressMessages(library("cowplot")))



# 读取数据并选择列
#df2 <- read.table("taxonomy_merge2.txt", header = TRUE)[, c(1:4, 6)]
df2 <- read.table(opts$input, header = TRUE)[, c(1:4, 6)]

#df2$Group <- "All"
# 排序 Phylum 列
df2 <- df2 %>%
  arrange(factor(Phylum, levels = unique(Phylum)))

# Check and convert to numeric if necessary
df2$Completeness <- as.numeric(as.character(df2$Completeness))
df2$Contamination <- as.numeric(as.character(df2$Contamination))

# Filter out NAs (if any)
df2 <- df2 %>%
  filter(!is.na(Completeness) & !is.na(Contamination))

# Check unique values
print(unique(df2$Completeness))  # Check unique values
print(unique(df2$Contamination))  # Check unique values
print(nrow(df2))  # Check number of rows


# 创建散点图
p1 <- ggplot(data = df2, aes(x = Completeness, y = Contamination, color = Phylum, size = Genome_Size)) +
  geom_point(alpha = 1) +
  geom_smooth(size = 1.3, color = "#4c6389", method = "loess") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size = 1, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 0.8)
  ) +
  scale_colour_manual(values = c("#f7fcfd", "#e5f5f9",  "#99d8c9", "#66c2a4", "#41ae76", "#238b45", "#006d2c",
                                 "#00441b", "#e7e1ef", "#d4b9da", "#c994c7", "#df65b0", "#e7298a", "#ce1256",
                                 "#980043", "#67001f", "#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8",
                                 "#807dba", "#6a51a3", "#54278f", "#3f007d", "#ffffcc", "#ffeda0", "#fed976",
                                 "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026")) +
  scale_size_continuous(range = c(1, 6)) +  # Adjust size range (e.g., 1 to 6)
  labs(
    y = "Contamination (%)",
    x = "Completeness (%)"
  )# +
  #theme(legend.position = "bottom")  # 隐藏图例，便于组合时不占用空间

p1_legend <- get_legend(p1)

# 创建 Completeness 方向的峰峦图 (顶部)
density_x <- ggplot(data = df2, aes(x = Completeness, fill = Phylum)) +
  geom_density(alpha = 0.5, color = "#969696") +
  scale_fill_manual(values = c("#f7fcfd", "#e5f5f9",  "#99d8c9", "#66c2a4", "#41ae76", "#238b45", "#006d2c",
                               "#00441b", "#e7e1ef", "#d4b9da", "#c994c7", "#df65b0", "#e7298a", "#ce1256",
                               "#980043", "#67001f", "#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8",
                               "#807dba", "#6a51a3", "#54278f", "#3f007d", "#ffffcc", "#ffeda0", "#fed976",
                               "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026")) +
  #theme_void() +  # 移除所有轴和刻度
  theme_classic()+
  scale_x_continuous(expand = c(0,0))+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20))  # 隐藏图例

# 创建 Contamination 方向的峰峦图 (右侧)
density_y <- ggplot(data = df2, aes(x = Contamination, fill = Phylum)) +
  geom_density(alpha = 0.5, color = "#969696") +
  scale_fill_manual(values = c("#f7fcfd", "#e5f5f9",  "#99d8c9", "#66c2a4", "#41ae76", "#238b45", "#006d2c",
                               "#00441b", "#e7e1ef", "#d4b9da", "#c994c7", "#df65b0", "#e7298a", "#ce1256",
                               "#980043", "#67001f", "#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8",
                               "#807dba", "#6a51a3", "#54278f", "#3f007d", "#ffffcc", "#ffeda0", "#fed976",
                               "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026")) +
  #theme_void() +  # 移除所有轴和刻度
  theme_classic()+
  scale_x_continuous(expand = c(0,0))+
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = 20),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 20)) +
  coord_flip()  # 翻转坐标轴，使其在右侧显示

# 组合图形：密度图放在上面和右边，散点图在中间
# combined_plot <- plot_grid(
#   density_x, NULL, p1, density_y,
#   ncol = 2, nrow = 2,
#   rel_heights = c(0.25, 1),  # 调整密度图和散点图的相对高度
#   rel_widths = c(1, 0.25)  # 调整密度图和散点图的相对宽度
# )

# combined_plot <- plot_grid(
#   density_x, NULL, p1, density_y,
#   ncol = 2, nrow = 2,
#   rel_heights = c(0.25, 1),  # 调整密度图和散点图的相对高度
#   rel_widths = c(1, 0.25)  # 调整密度图和散点图的相对宽度
# )

combined_plot <- ggdraw(plot_grid(
  plot_grid(density_x, p1+theme(legend.position = 'none'), ncol = 1, rel_widths = c(1, 2), rel_heights = c(1, 2)),
  plot_grid(p1_legend, density_y, ncol = 1, rel_heights = c(1, 2)), rel_widths = c(2, 1)))


# 打印组合图
#print(combined_plot)

# 保存图像
#save_plot("combined_plot.pdf", combined_plot, base_width = 9, base_height = 8)
ggsave(paste(opts$output, "MAGs_completeness_contamination.pdf", sep=""), combined_plot, width = 89*3, height = 59*4, unit = "mm")

