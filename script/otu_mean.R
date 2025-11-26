#!/usr/bin/env Rscript

# Copyright 2016-2021 Yong-Xin Liu <yxliu@genetics.ac.cn / metagenome@126.com>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, Yang Bai. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein Cell 2021(12) 5:315-330 doi: 10.1007/s13238-020-00724-8

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 1.1 程序功能描述和主要步骤

# 程序功能：OTU按实验设计分组，计算均值
# Functions: Calculate mean OTU abundance for each group
# 主要步骤Main steps: 
# - 读取输出表和元数据
# - 交叉筛选，按降序排序
# - (可选)重新标准化
# - 按丰度筛选
# - 计算分组求均值
# - 保存结果


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
    make_option(c("-i", "--input"), type="character", default="otutab.txt",
                help="Feature table [default %default]"),
    make_option(c("-T", "--thre"), type="numeric", default=0,
                help="Threshold of abundance, such as 0.001 [default %default]"),
    make_option(c("-s", "--scale"), type="logical", default=FALSE,
                help="Normalize by itself [default %default]"),
    make_option(c("-z", "--zoom"), type="numeric", default=100,
                help="Normalize to 100 [default %default]"),
    make_option(c("-d", "--metadata"), type="character", default="metadata.txt",
                help="Experiment metadata or sample metadata [default %default]"),
    make_option(c("-a", "--all"), type="character", default="FALSE",
                help="All sample average or summary column [default %default]"),
    make_option(c("-n", "--group"), type="character", default="GroupID",
                help="Column name of group [default %default]"),
    make_option(c("-t", "--type"), type="character", default="mean",
                help="mean, sum [default %default]"),
    make_option(c("-o", "--output"), type="character", default="otutab_mean.txt",
                help="output directory and prefix [default %default]")
  )
  opts <- parse_args(OptionParser(option_list=option_list))

  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$output==""){
    opts$output=paste(opts$input, "_", opts$group, sep = "")}

  # 调试参数区，完成后请注释掉
  # setwd("/mnt/m1/yongxin/rice/MAG/genome/isolate")
  # opts$input="result/dbcan2/sum.Level1.raw.txt"
  # opts$metadata="result/itol/all_tax.txt"
  # opts$group="Phylum"
  # opts$type="sum"
  # opts$output="result/dbcan2/sum.Level1.sum.Phylum.txt"

  # 显示输入输出确认是否正确
  print(paste("Feature table: ", opts$input,  sep = ""))
  print(paste("Metadata: ", opts$metadata,  sep = ""))
  print(paste("Group name: ", opts$group,  sep = ""))
  print(paste("Abundance threshold: ", opts$thre,  sep = ""))
  print(paste("Calculate type: ", opts$type,  sep = ""))
  print(paste("Output filename: ", opts$output, sep = ""))
}



# 0. 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
package_list <- c("dplyr")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 1. 读取OTU表
otutab = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)
# head(otutab)

# 读取实验设计
metadata = read.table(opts$metadata, header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F) 

# 将选定的分组列统一命名为group
metadata$group=metadata[,opts$group]

# 交叉筛选
idx = rownames(metadata) %in% colnames(otutab)
# table(idx)
metadata = metadata[idx,,drop=F]
otutab = otutab[,rownames(metadata)]
# 按丰度排序
idx = order(rowMeans(otutab), decreasing = T)
otutab = otutab[idx,]

# 标准化，并筛选高丰度菌均值最小百万分之一0.0001%
# 条件判断是否标准化
if (opts$scale){
  norm = t(otutab)/colSums(otutab,na=T)*opts$zoom
}else{
  norm = t(otutab)
}
# rowSums(norm)
idx = colMeans(norm) > opts$thre
HA = norm[,idx]
# dim(HA)
# rowSums(HA)

# 按group合并
merge=cbind(HA, metadata[,c("group"),drop=F])
if (opts$type == "sum"){
  HA_group_mean = merge %>% group_by(group) %>% summarise_all(sum)
}else {
  HA_group_mean = merge %>% group_by(group) %>% summarise_all(mean)
}

if (opts$all){
  HA_t = as.data.frame(cbind(c("All", round(colMeans(HA), digits = 6)),t(HA_group_mean)), stringsAsFactors = F)
}else {
  HA_t = as.data.frame(t(HA_group_mean), stringsAsFactors = F)
}

rownames(HA_t)[1] = "OTUID"

write.table(HA_t, file=paste(opts$output, "", sep = ""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=F)

