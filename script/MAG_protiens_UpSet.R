#!/usr/bin/env Rscript

# Copyright 2024 De-feng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, et al. 2021. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12: 315-330. https://doi.org/10.1007/s13238-020-00724-8

# 手动运行脚本，使用 Ctrl+Shift+H 或 Session 需要设置工作目录
# Set Work Directory - Choose Directory / To Source File Location

# Clean enviroment object
rm(list=ls()) 

# 1.1 简介 Introduction #----

# 程序功能：MAGs不同功能数据库注释结果比较
# Functions: Comparison of annotation results of MAGs in different functional databases


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
    make_option(c("-i", "--input"), type="character", default="result/eggnog/d3.data4venn_hope4.txt",
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



#Data4Pic = read.table("data/d3.data4venn_hope4.txt", header=T, row.names=1)
Data4Pic = read.table(opts$input, header=T, row.names=1)
Data4Pic[Data4Pic>0] = 1

# 使用Queries参数对集合图进行修饰
#pdf(file="results/p3.GenusUpsetIndiv112.pdf", width=9, height=5, pointsize=8)
pdf(file=paste(opts$output, "MAGs_proteins_UpSet.pdf", sep=""), width=9, height=5, pointsize=8)
queries<-list(list(query=intersects,
                   params=list("COGs"),
                   active=T,
                   color="#E41A1C"),
              list(query=intersects,
                   params=list("KEGG"),
                   active=T,
                   color="#377EB8"),
              list(query=intersects,
                   params=list("ECs"),
                   active=T,
                   color="#984EA3"),
              list(query=intersects,
                   params=list("GOs"),
                   active=T,
                   color="#4DAF4A")
)
p3 <- UpSetR::upset(Data4Pic, sets = c("COGs","KEGG","ECs","GOs"), mb.ratio = c(0.6, 0.4), 
                    scale.intersections ="log10",
                    order.by = "freq",
                    queries = queries,
                    nintersects=NA,
                    decreasing = T,
                    empty.intersections = TRUE,
                    nsets = 4, number.angles = 0, point.size = 4, line.size = 1,
                    main.bar.color = "gray23",
                    sets.bar.color = c("#E41A1C", "#377EB8", "#984EA3", "#4DAF4A"),
                    #sets.bar.color = "gray23",
                    #sets.bar.color = c("#006bba", "#d4550e", "#00a6dc", "#67a330"),
                    #mainbar.y.label = "Number of shared genus",
                    sets.x.label = "Set size (x105)", text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5)
)
p3
dev.off()

