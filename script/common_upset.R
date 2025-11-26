#!/usr/bin/env Rscript

# Copyright 2024 De-feng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, et al. 2021. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12: 315-330. https://doi.org/10.1007/s13238-020-00724-8

# 手动运行脚本，使用 Ctrl+Shift+H 或 Session 需要设置工作目录
# Set Work Directory - Choose Directory / To Source File Location

# Clean enviroment object
rm(list=ls()) 

# 1.1 简介 Introduction #----

# 程序功能：功能通路注释结果分组比较UpSet图
# Functions: UpSet diagram for group comparison of functional pathway annotation results


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
    make_option(c("-i", "--input"), type="character", default="result/eggnog/d3.data4venn_hope8.txt",
                help="Unfiltered OTU table [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/eggnog/",
                help="Output MAGs UpSet plot for different databases  [default %default]") 
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
print("You are using the following parameters:")
print(opts)



# Install related packages
if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggplot2","colorspace","RColorBrewer","UpSetR","VennDiagram",
             "formattable")) # ,"vegan"
}
# load related packages
suppressWarnings(suppressMessages(library("ggplot2")))
suppressWarnings(suppressMessages(library("colorspace")))
suppressWarnings(suppressMessages(library("RColorBrewer")))
# suppressWarnings(suppressMessages(library("vegan")))
suppressWarnings(suppressMessages(library("UpSetR")))
suppressWarnings(suppressMessages(library("VennDiagram")))
suppressWarnings(suppressMessages(library("formattable")))


# 读入数据
Data4Pic = read.table(opts$input, header=T, row.names=1)
Data4Pic[Data4Pic>0] = 1

# 为方便代码兼容R3.6和R4.0版本，再进行一次数据处理
# 目的是将'Data4Pic'的行名转换为'numeric'的行号
row.names(Data4Pic) <- 1:nrow(Data4Pic)

# 上方柱状图显示否一分类级别的比例
pdf(file=paste(opts$output, "MAGs_proteins_UpSet_single_dataset.pdf", sep=""), width=9, height=5, pointsize=8)
p7 <- UpSetR::upset(Data4Pic, sets = colnames(Data4Pic), mb.ratio = c(0.6, 0.4), #order.by = "freq",
                    order.by = "degree",
                    empty.intersections = TRUE,
                    nsets = 4, number.angles = 0, point.size = 2, line.size = 0.6, 
                    sets.bar.color=brewer.pal(4,"Set1"),
                    shade.color="#80cdc1",
                    sets.x.label = "Set size", text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5))
p7
dev.off()


