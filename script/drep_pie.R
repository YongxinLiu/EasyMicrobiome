#!/usr/bin/env Rscript

# Copyright 2016-2022 Yong-Xin Liu <yxliu@genetics.ac.cn / metagenome@126.com>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, Yang Bai. 
# A practical guide to amplicon and metagenomic analysis of microbiome data. 
# Protein Cell 2021(12) 5:315-330 doi: 10.1007/s13238-020-00724-8

# 手动运行脚本请，需要设置工作目录
# 使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 1.1 程序功能描述和主要步骤

# 程序功能：读取文件表中的某列，绘制直方图
# Functions: Calculate one line to plot histogram
# 主要步骤Main steps: 
# - 读取输出表和元数据
# - 保存图片

options(warn = -1) # Turn off warning

# 1.2 解析命令行
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T) 
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="result/drep/0.95/stat.txt",
                help="Feature table [default %default]"),
    make_option(c("-d", "--metadata"), type="character", default="metadata.txt",
                help="Experiment metadata or sample metadata [default %default]"),
    make_option(c("-n", "--group"), type="character", default="cluster_members",
                help="Column name of group [default %default]"),
    make_option(c("-t", "--transform"), type="character", default="FALSE",
                help="Column name of group [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output directory and prefix [default %default]")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
  
  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$output==""){
    opts$output=paste(opts$input, "_", opts$group, "_pie", sep = "")}
  
  # 调试参数区，完成后请注释掉
  # setwd("/mnt/m3/yongxin/rice/binning")
  # opts$input="result/drep/ANI95_stat.txt"
  # opts$transform="log10"

    # 显示输入输出确认是否正确
  print(paste("Feature table: ", opts$input,  sep = ""))
  print(paste("Metadata: ", opts$metadata,  sep = ""))
  print(paste("Group name: ", opts$group,  sep = ""))
  print(paste("transform: ", opts$transform,  sep = ""))
  print(paste("Output filename: ", opts$output, sep = ""))
}

# 0. 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
package_list <- c("dplyr","ggplot2")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 1. 读取文件并绘图
# Load require packages
library(ggplot2)
library(dplyr)

# Prepare data
df = read.table(opts$input, row.names = 1, header = T, sep = "\t")
df = df[,c(opts$group),drop=F]
MAG1=dim(df[df==1, , drop=F])[1]
MAG2=dim(df[df==2, , drop=F])[1]
MAG2plus=dim(df[df>2, , drop=F])[1]

# Create Data
data <- data.frame(
  group=c("Single","Double","> 2 genomes"),
  value=c(MAG1,MAG2,MAG2plus)
)

# Compute the position of labels
data <- data %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
p = ggplot(data, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  geom_text(aes(y = ypos, label = group), color = "white", size=5) +
  scale_fill_brewer(palette="Set1")
ggsave(paste(opts$output,".pdf", sep=""), p, width = 89, height = 59, units = "mm")
write.table(data, file=paste(opts$output,".txt",sep=""), append = F, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)