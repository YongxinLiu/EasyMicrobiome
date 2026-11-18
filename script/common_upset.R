#!/usr/bin/env Rscript

# Copyright 2026 De-feng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, et al. 2021. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12: 315-330. https://doi.org/10.1007/s13238-020-00724-8

# 手动运行脚本，使用 Ctrl+Shift+H 或 Session 需要设置工作目录
# Set Work Directory - Choose Directory / To Source File Location

# Clean enviroment object
rm(list=ls()) 

# 1.1 简介 Introduction #----

# 程序功能：功能通路注释结果分组比较UpSet图（ComplexUpset版本）
# Functions: UpSet diagram for group comparison of functional pathway annotation results

options(warn = -1) # Turn off warning

# 1.2 参数 Parameters #----
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"

if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}

if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="result/eggnog/d3.data4venn_hope8.txt",
                help="Unfiltered OTU table [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/eggnog/",
                help="Output MAGs UpSet plot for different databases  [default %default]") 
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
print("You are using the following parameters:")
print(opts)


# 安装并加载 ComplexUpset（如未安装）
if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
  install.packages("ComplexUpset", repos = site)
}
suppressWarnings(suppressMessages(library(ComplexUpset)))
suppressWarnings(suppressMessages(library(ggplot2)))

# 读入数据
Data4Pic = read.table(opts$input, header=T, row.names=1)
Data4Pic[Data4Pic > 0] = 1

# ComplexUpset 需要数据框，且行名不再需要转换为数字
# 将二值化的矩阵转为数据框
Data4Pic_df <- as.data.frame(Data4Pic)

# 设定要用于交集分析的列（即所有列）
sets_to_intersect <- colnames(Data4Pic_df)

# 绘图并保存
pdf(
  file = paste(opts$output,
               "MAGs_proteins_UpSet_single_dataset.pdf",
               sep = ""),
  width = 9,
  height = 5,
  pointsize = 8,
  onefile = FALSE
)

p7 <- ComplexUpset::upset(
  Data4Pic_df,
  intersect = sets_to_intersect,
  # 按交集度（degree）排序，对应原 order.by = "degree"
  sort_intersections_by = "degree",
  # 显示空交集
  min_size = 0,
  keep_empty_groups = TRUE,
  # 矩阵与柱状图的高度比例，对应 mb.ratio = c(0.6, 0.4)
  height_ratio = 0.4,
  # 集合名称标签
  name = "Set size",
  # 集合大小条形图颜色（红色）
  set_sizes = (
    upset_set_size(geom = geom_bar(fill = "#E41A1C", width = 0.6)) +
      theme(axis.text.x = element_text(size = 10))
  ),
  # 交集矩阵：点大小与线宽
  matrix = intersection_matrix(
    geom = geom_point(size = 2, stroke = 0.6)
  ),
  # 交集大小柱状图的样式
  base_annotations = list(
    'Intersection size' = intersection_size(
      counts = FALSE,
      text = list(size = 3),
      fill = "#80cdc1"
    ) + theme(
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 11)
    )
  ),
  # 整体主题调整
  themes = upset_modify_themes(
    list(
      'intersections_matrix' = theme(
        axis.text.y = element_text(size = 11, colour = "black")
      ),
      'overall_sizes' = theme(
        axis.text.x = element_text(size = 10)
      )
    )
  ),
  # 调整列宽比例，给集合名称留出空间
  width_ratio = 0.25
)

print(p7)

dev.off()