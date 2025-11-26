#!/usr/bin/env Rscript

# Copyright 2024 De-feng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, et al. 2021. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12: 315-330. https://doi.org/10.1007/s13238-020-00724-8

# 手动运行脚本，使用 Ctrl+Shift+H 或 Session 需要设置工作目录
# Set Work Directory - Choose Directory / To Source File Location

# Clean enviroment object
rm(list=ls()) 

# 1.1 简介 Introduction #----

# 程序功能：不同层级功能通路和KO基因桑基图展示
# Functions: Sankey diagram display of functional pathways and KO genes at different levels



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
    make_option(c("-i", "--input"), type="character", default="result/humann3/data_sankey7.txt",
                help="MAGs taxonomy tables [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/humann3/",
                help="Output Sankey diagram display of functional pathways and KO genes at different levels [default %default]") 
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
print("You are using the following parameters:")
print(opts)


# Install related packages
if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggsankey","ggplot2","ggalluvial"))
}
# load related packages
suppressWarnings(suppressMessages(library("ggsankey")))
suppressWarnings(suppressMessages(library("ggplot2")))
suppressWarnings(suppressMessages(library("ggalluvial")))


#df01 <- read.table(file = "data/data_sankey7.txt", sep = "\t", header = T, check.names = FALSE)
df01 <- read.table(file = opts$input, sep = "\t", header = T, check.names = FALSE)
data <- df01
df <- to_lodes_form(data[,1:ncol(data)],
                    axes = 1:ncol(data),
                    id = "value")

# Set color
col<- rep(c('#3690c0', '#f16913', '#238b45', '#ff6f81', '#fc9272', '#ffc2c0','#8c96c6',
            '#bfd3e6', '#fae6f0', '#eb6fa6', '#ff88b5', '#00b1a5',"#ffa68f","#ffca75","#ccebc5","#7bccc4",
            "#6baed6","#2171b5","#c6dbef","#448c99","#67a9cf","#b8d8c9","#88419d","#d4b9da","#fee6ce",
            "#8f9898","#bfcfcb"), 6)


# Sankey diagram
p3 <- ggplot(df, aes(x = x, fill=stratum, label=stratum,
                     stratum = stratum, alluvium  = value), width = 0.1)+
  geom_flow(width = 0.1,
            curve_type = "sine",
            alpha = 0.6,
            color = 'white',
            size = 0.05)+
  geom_stratum(width = 0.1, color = "white")+
  geom_text(stat = 'stratum', size = 3.5, color = 'black')+
  scale_fill_manual(values = col)+
  theme_void()+
  theme(legend.position = 'none',
        text = element_text(size = 18))
#p3

pdf(file=paste(opts$output, "MAGs_KEGG_Sankey.pdf", sep=""), height = 5.2, width = 11)
p3
dev.off()

