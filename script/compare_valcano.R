#!/usr/bin/env Rscript

# Copyright 2024 De-feng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, et al. 2021. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12: 315-330. https://doi.org/10.1007/s13238-020-00724-8

# 手动运行脚本，使用 Ctrl+Shift+H 或 Session 需要设置工作目录
# Set Work Directory - Choose Directory / To Source File Location

# Clean enviroment object
rm(list=ls()) 

# 1.1 简介 Introduction #----

# 程序功能：运用MaAsLin2方法计算组间差异
# Functions: Difference analysis using MaAsLin2


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
    make_option(c("-i", "--input"), type="character", default="result12/metaphlan4/MaAsLin2_overall_difference.csv",
                help="Unfiltered OTU table [default %default]"),
    make_option(c("-d", "--trend"), type="character", default="result12/metaphlan4/MaAsLin2_enriched_depleted.csv",
                help="metadata file or metadata [default %default]"),
    make_option(c("-o", "--output"), type="character", default="metaphlan4/",
                help="Output quantile value for filter feature table [default %default]") 
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
print("You are using the following parameters:")
print(opts)


# Install related packages
if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggplot2","magrittr","dplyr")) 
}
# load related packages
suppressWarnings(suppressMessages(library("ggplot2")))
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("magrittr")))

data_species01 <- read.csv(opts$trend, row.names = 1)
data_species01$species <- rownames(data_species01)
data_species01 <- data_species01[, -c(1:4)]

data_MWAS <- read.csv(opts$input, row.names = 1)
data_MWAS$species <- rownames(data_MWAS)

data_species02 <- merge(data_species01, data_MWAS, by = "species")

#devtools::install_github("BioSenior/ggvolcano", force = TRUE)
library(ggVolcano)
data_vol <- data_species02
data_vol = as.data.frame(data_vol)

data_vol2 <- data_vol
data_vol2$padj2 <- -log10(data_vol2$FDR)

#logFC = 0.5
#P.Value = 0.05
library(ggplot2)
p_volcano1 <- ggplot(data = data_vol2, aes(x = -Beta, y = padj2)) +
  geom_point(alpha = 0.4, size = 2.0, aes(color = NPC.association)) + 
  ylab("-log10(Pvalue)") +
  scale_color_manual(values = c("#74add1","#a60026", "grey")) + 
  geom_vline(xintercept = 0, lty = 4, col = "black", lwd = 0.4) + 
  geom_hline(yintercept = -log10(0.75), lty = 4, col = "black", lwd = 0.4) + 
  labs(x = "Coef. (by MaAsLin2)", y= bquote(atop(-Log[10]~italic(FDR))))+
  theme_bw()

# add labels
library(dplyr)
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
ggsave(paste0(opts$output,"Difference_valcano.pdf"), p_volcano2, width=59 * 1.5, height=50 * 1.5, unit='mm')
#p_volcano2

