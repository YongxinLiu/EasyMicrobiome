#!/usr/bin/env Rscript

# Copyright 2024 De-feng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, et al. 2021. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12: 315-330. https://doi.org/10.1007/s13238-020-00724-8

# 手动运行脚本，使用 Ctrl+Shift+H 或 Session 需要设置工作目录
# Set Work Directory - Choose Directory / To Source File Location

# Clean enviroment object
rm(list=ls()) 

# 1.1 简介 Introduction #----

# 程序功能：处理Metaphlan4得到的物种组成数据作为SparCC计算的输入
# Functions: Processing of species composition data obtained from Metaphlan4 as input for SparCC calculations


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
    make_option(c("-i", "--input"), type="character", default="result12/metaphlan4/Species.txt",
                help="Metaphlan4 species table"),
    make_option(c("-g", "--group"), type="character", default="result12/metadata.txt",
                help="Unfiltered OTU table [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result12/metaphlan4/",
                help="Output file for SparCC analysis in different groups  [default %default]") 
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
print("You are using the following parameters:")
print(opts)


# Install related packages
if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("reshape2","ggplot2","ggprism","dplyr","plyr",
             "igraph")) # ,"vegan"
}
# load related packages
suppressWarnings(suppressMessages(library("reshape2")))
suppressWarnings(suppressMessages(library("ggplot2")))
suppressWarnings(suppressMessages(library("ggprism")))
# suppressWarnings(suppressMessages(library("vegan")))
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("plyr")))
suppressWarnings(suppressMessages(library("igraph")))


# metadata 
#design <- read.table(file = "result12/metaphlan4/metadata.txt", sep = "\t", header = T, row.names=1)
design <- read.table(opts$group, sep = "\t", header = T, row.names=1)
#df3 <- read.table(file = "result12/metaphlan4/Species.txt", sep = "\t", header = T, check.names = FALSE)
df3 <- read.table(opts$input, sep = "\t", header = T, check.names = FALSE)
df_species <- df3

# sum of Species
data<-aggregate(.~ Taxonomy,data=df_species,sum)
rownames(data) = data$Taxonomy
data = data[, -1]

#If species are too much, please select part of them according to relative abundace
#species_selected <- read.table(file = "data/Species_selected.txt", sep = "\t", header = T, check.names = FALSE)
#rownames(species_selected) <- species_selected$ID
#data <- data[rownames(data)%in%rownames(species_selected), ]

data3 = apply(data, 2, function(x)x/100)
#sum(data3$H001F)

# Check data
dim(data3)
data3 = data3 * 100000
OTU.table.filtered.colnames <- colnames(data3)
OTU.table.filtered.sparcc <- cbind(rownames(data3), data3)
colnames(OTU.table.filtered.sparcc) <- c("OTU_id", OTU.table.filtered.colnames)
OTU.table.filtered.sparcc2 <- t(OTU.table.filtered.sparcc)
OTU.table.filtered.sparcc2 <- OTU.table.filtered.sparcc2[-1,]
OTU.table.filtered.sparcc2 <- as.data.frame(OTU.table.filtered.sparcc2)
#OTU.table.filtered.sparcc2$group <- rownames(OTU.table.filtered.sparcc2)
#OTU.table.filtered.sparcc2$group = as.character(OTU.table.filtered.sparcc2$group)
#OTU.table.filtered.sparcc2$group = gsub("[0-9]","", OTU.table.filtered.sparcc2$group)
otutab <- as.data.frame(t(OTU.table.filtered.sparcc2))

# Select by manual set group
# NPC group
if (TRUE){
  sub_design = subset(design, Group %in% c("Cancer")) 
  sub_design$Group  = factor(sub_design$Group, levels=c("Cancer"))
}
idx = rownames(sub_design) %in% colnames(otutab)
sub_design_Cancer = sub_design[idx,]
sub_otutab_Cancer = otutab[,rownames(sub_design_Cancer)]
sub_otutab_Cancer <- as.data.frame(sub_otutab_Cancer)
#write.table(sub_otutab_npc, file = "results/Species_sparcc_npc_p01_11R_adjusted_1031.txt", row.names = T, sep = "\t", quote = T, col.names = T)
write.table(sub_otutab_Cancer, file = paste(opts$output, "Cancer_sparcc.txt", sep=""), row.names = T, sep = "\t", quote = T, col.names = T)

# Healthy group
if (TRUE){
  sub_design = subset(design, Group %in% c("Normal")) 
  sub_design$Group  = factor(sub_design$Group, levels=c("Normal"))
}
idx = rownames(sub_design) %in% colnames(otutab)
sub_design_Normal = sub_design[idx,]
sub_otutab_Normal = otutab[,rownames(sub_design_Normal)]
sub_otutab_Normal = as.data.frame(sub_otutab_Normal)
sub_otutab_Normal = as.data.frame(sub_otutab_Normal)
#write.table(sub_otutab_healthy, file = "results/Species_sparcc_npc_h01_11R_adjusted_1031.txt", row.names = T, sep = "\t", quote = T, col.names = T)
write.table(sub_otutab_Normal, file = paste(opts$output, "Normal_sparcc.txt", sep=""), row.names = T, sep = "\t", quote = T, col.names = T)
