#!/usr/bin/env Rscript

# Copyright 2024-2026 Defeng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Bai, et al. 2025. EasyMetagenome: A User‐Friendly and Flexible Pipeline for Shotgun Metagenomic Analysis in Microbiome Research. iMeta 4: e70001. https://doi.org/10.1002/imt2.70001

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 更新
# 2024/11/12：增加利用MaAsLin2差异分析结果绘制火山图代码
# 2025/11/27：规范代码

# 1.1 程序功能描述和主要步骤

# 程序功能：火山图展示组间差异
# Functions: Species difference illustrated using volcano plot

options(warn = -1) # Turn off warning


## 设置输入输出文件和参数

# 修改下面`default=`后面的文件和参数。
#
# 输入文件为MaAsLin2组间差异分析得到的两个结果文件(metaphlan4/MaAsLin2_overall_difference.csv)+趋势信息(metaphlan4/MaAsLin2_enriched_depleted.csv)
#
# 输入文件"-i", "--input"，metaphlan4/MaAsLin2_overall_difference.csv; MaAsLin2差异分析结果汇总信息；
#
# 实验设计"-d", "--trend"，默认`metaphlan4/MaAsLin2_enriched_depleted.csv`，差异分析结果趋势信息；
#
# 图片宽"-w", "--width"，默认89 mm，根据图像布局可适当增大或缩小
#
# 图片高"-e", "--height"，默认59 mm，根据图像布局可适当增大或缩小
#
# 分组列名"-o", "--output"，默认为输出目录；


# 1.2 依赖包安装

site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
a = rownames(installed.packages())

# install CRAN
install_CRAN <- c("ggplot2", "BiocManager", "optparse", "dplyr", "magrittr")
for (i in install_CRAN) {
  if (!i %in% a)
    install.packages(i, repos = site)
  require(i,character.only=T)
  a = rownames(installed.packages())
}

# install bioconductor
install_bioc <- c("ggplot2", "multcompView")
for (i in install_bioc) {
  if (!i %in% a)
    BiocManager::install(i, update = F)
  a = rownames(installed.packages())
}

# install github
if (!"amplicon" %in% a){
  devtools::install_github("microbiota/amplicon")
}

if (!"ggvolcano" %in% a){
  devtools::install_github("BioSenior/ggvolcano", force = TRUE)
}


# 1.3 解析命令行
# 设置清华源加速下载
# (Optional) Set up Tsinghua Mirror to speed up download
# site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析包是否安装，没安装则安装，然后加载
# Determine whether the command line parsing package is installed, install it if it is not installed, then load
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="metaphlan4/MaAsLin2_overall_difference.csv",
                help="Unfiltered OTU table [default %default]"),
    make_option(c("-d", "--trend"), type="character", default="metaphlan4/MaAsLin2_enriched_depleted.csv",
                help="metadata file or metadata [default %default]"),
    make_option(c("-o", "--output"), type="character", default="metaphlan4/",
                help="Output quantile value for filter feature table [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=89,
                help="Width of figure [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=59,
                help="Height of figure [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
print("You are using the following parameters:")
print(opts)


# # Install related packages
# if (FALSE){
#   source("https://bioconductor.org/biocLite.R")
#   biocLite(c("ggplot2","magrittr","dplyr")) 
# }
# # load related packages
# suppressWarnings(suppressMessages(library("ggplot2")))
# suppressWarnings(suppressMessages(library("dplyr")))
# suppressWarnings(suppressMessages(library("magrittr")))

# 2. 依赖关系检查、安装和加载

# 依赖包列表
package_list <- c(
  "ggplot2", "BiocManager", "optparse", "dplyr", "magrittr", "ggvolcano"
)

# 批量安装和加载
for (p in package_list) {
  # 如果未安装，则安装
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  }
  
  # 批量加载，抑制警告和消息
  suppressWarnings(
    suppressMessages(
      library(p, character.only = TRUE)
    )
  )
}


# 读取输入文件

# 读取趋势表格
data_species01 <- read.csv(opts$trend, row.names = 1)
data_species01$species <- rownames(data_species01)
data_species01 <- data_species01[, -c(1:4)]

# 读取差异结果汇总表
data_MWAS <- read.csv(opts$input, row.names = 1)
data_MWAS$species <- rownames(data_MWAS)

data_species02 <- merge(data_species01, data_MWAS, by = "species")

# devtools::install_github("BioSenior/ggvolcano", force = TRUE)
# library(ggVolcano)
data_vol <- data_species02
data_vol = as.data.frame(data_vol)

data_vol2 <- data_vol
data_vol2$padj2 <- -log10(data_vol2$FDR)

# logFC = 0.5
# P.Value = 0.05
# library(ggplot2)
p_volcano1 <- ggplot(data = data_vol2, aes(x = -Beta, y = padj2)) +
  geom_point(alpha = 0.4, size = 2.0, aes(color = NPC.association)) + 
  ylab("-log10(Pvalue)") +
  scale_color_manual(values = c("#74add1","#a60026", "grey")) + 
  geom_vline(xintercept = 0, lty = 4, col = "black", lwd = 0.4) + 
  geom_hline(yintercept = -log10(0.75), lty = 4, col = "black", lwd = 0.4) + 
  labs(x = "Coef. (by MaAsLin2)", y= bquote(atop(-Log[10]~italic(FDR))))+
  theme_bw()

# add labels
# library(dplyr)
# select top 5 enriched species
up_data1 <- filter(data_vol2, data_vol2$NPC.association == "enriched")
up_data2 <- arrange(up_data1, desc(up_data1$padj2))
up_data_5 <- up_data2[1:5, ]

# select top 25 depleted species
down_data1 <- filter(data_vol2, data_vol2$NPC.association == "depleted")
down_data2 <- arrange(down_data1, desc(down_data1$padj2))
down_data_5 <- down_data2[1:5, ]

# using geom_text_repel() to add labels
library(ggrepel)
p_volcano2 <- p_volcano1 +  
  geom_text_repel(data = up_data_5, aes(x = -Beta, 
                                        y = padj2, 
                                        label = up_data_5$Feature), size = 2, max.overlaps = 20) + 
  geom_text_repel(data = down_data_5, aes(x = -Beta,
                                          y = padj2,
                                          label = down_data_5$Feature), size = 2, max.overlaps = 20)+
  theme(legend.position = c(0.84, 0.85),panel.grid = element_blank())
#ggsave(paste0(opts$output,"Difference_valcano.pdf"), p_volcano2, width=59 * 1.5, height=50 * 1.5, unit='mm')
ggsave(paste0(opts$output,"Difference_valcano.pdf"), p_volcano2, width=opts$width, height=opts$height, unit='mm')
