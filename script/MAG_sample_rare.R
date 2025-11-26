#!/usr/bin/env Rscript

# Copyright 2024 De-feng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, et al. 2021. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12: 315-330. https://doi.org/10.1007/s13238-020-00724-8

# 手动运行脚本，使用 Ctrl+Shift+H 或 Session 需要设置工作目录
# Set Work Directory - Choose Directory / To Source File Location

# Clean enviroment object
rm(list=ls()) 

# 1.1 简介 Introduction #----

# 程序功能：MAG数量和样本数量的关系稀疏曲线
# Functions: The rare curve for the relationship between MAGs and samples


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
    make_option(c("-i", "--input"), type="character", default="binning/result/coverm/abundance.tsv",
                help="Unfiltered OTU table [default %default]"),
    make_option(c("-t", "--taxonomy"), type="character", default="binning/result/coverm/taxonomy_MAG.txt",
                help="metadata file or metadata [default %default]"),
    make_option(c("-o", "--output"), type="character", default="binning/result/coverm/",
                help="Output quantile value for filter feature table [default %default]") 
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
print("You are using the following parameters:")
print(opts)

# Version 1.0, Based on ASV table and taxonomy, output well info (purity, counts and taxonomy) and candidate wells of non-redundancy ASV
# 版本 1.0, 基于ASV表和7级物种注释文件，输出每个孔的信息，筛选每个孔的信息()，以及非冗余ASV的修行孔，纯度优先，数据量其次排序的Top 5


# Install related packages
if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggplot2","grid","scales","dplyr")) # ,"vegan"
}
# load related packages
suppressWarnings(suppressMessages(library("ggplot2")))
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("scales")))
# suppressWarnings(suppressMessages(library("vegan")))
suppressWarnings(suppressMessages(library("grid")))

# Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=7),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=7),
                   text=element_text(family="sans", size=7))

# 2.1 读文件 Load #----

# Load ASV table
# Public file 1. "result/ASV_table.txt"  raw reads count of each ASV in each sample
otu_table = read.delim(opts$input, row.names= 1,  header=T, sep="\t")

# Load ASV metadata
# Public file 2. "result/taxonomy_8.txt"  taxonomy for each ASV, tab seperated
#taxonomy = read.delim("result/taxonomy_8.txt", row.names= 1,header=T, sep="\t")
taxonomy = read.delim(opts$taxonomy, row.names= 1,header=T, sep="\t")
taxonomy$Full=paste(taxonomy$Phylum,taxonomy$Class,taxonomy$Order,taxonomy$Family,taxonomy$Genus,taxonomy$Species,sep = ";")

# Extract only those samples in common between the two tables
idx = rownames(otu_table) %in% rownames(taxonomy)
otu_table = otu_table[idx,]
taxonomy = taxonomy[rownames(otu_table),]


# 2.2 稀疏曲线 Rarefraction curve #----

sample_rare  =  function(df, count_cutoff = 3, length = 30, rep = 30){  
  result = data.frame(sample = c(0), richness = c(0))
  count=df
  # Set otu table to binary for easy calculate richness
  count[count < count_cutoff] = 0
  count[count >= count_cutoff] = 1
  # Sample number
  n = dim(count)[2]
  # 设置X轴各取样数量
  x = unique(as.integer(seq(1, n, length.out = length)))
  for (i in x){
    m = choose(n, i)
    if (m > rep){m = rep}
    # loop list and calculate richness
    for(j in 1:m) {
      idx = sample(n, i)
      temp = count[, idx, drop = F]
      # subset non-zero row
      temp1=temp[rowSums(temp)>0, , drop=F]
      # row number is accumulative OTUs
      result = rbind(result, c(i, dim(temp1)[1]))
    }
  }
  # remove start 0,0
  result = result[-1,]
  # factor draw as each box
  result$sample = factor(result$sample, levels = unique(result$sample))
  return(result)
}
rare_box = sample_rare(otu_table, count_cutoff = 1, length = 30, rep = 30)
p = ggplot(rare_box,aes(x=sample,y=richness, fill=sample))+
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.5, size = 0.1)+main_theme+ xlab("Sample numbers")+ylab("MAG numbers") + 
  theme(legend.position = "NA",axis.text.x = element_text(angle = 45,vjust=1, hjust=1))
p
ggsave(paste(opts$output, "MAG_rare_curve.pdf", sep=""), p, width = 69*2, height = 79, unit = "mm")
ggsave(paste(opts$output, "MAG_rare_curve.png", sep=""), p, width = 69*2, height = 79, unit = "mm")

print(paste("Rarefaction curve in ",opts$output, "MAG_rare_curve.pdf/png", sep=""))

