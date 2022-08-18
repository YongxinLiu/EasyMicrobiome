#!/usr/bin/env Rscript

# Copyright 2016-2022 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).


# 1. 分析前准备：帮助、参数、依赖包和读取文件

# 命令行运行为当前目录；Rstudio手动运行脚本前需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 1.1 程序功能描述和主要步骤

# 程序功能：转换Salmon计数表(含小数)为OTU类整数表
# Script functions: Format salmon count table to OTU-like table
# Main steps:
# - Read data table result/salmon/gene.count
# - table + 0.5, then int
# - Save in result/otutab.txt

# 程序使用示例
# USAGE
# Default
# # 显示脚本帮助
# Rscript ./script/format_salmon2otutab.R -h
# # 完整参数，输出文件名为alpha指数类型
# Rscript ./script/format_salmon2otutab.R -i result/salmon/gene.count \
# -o result/otutab.txt
# options(warn = -1) # Turn off warning



# 1.2 解析命令行
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
option_list <- list(
  make_option(c("-i", "--input"), type="character", default="result/salmon/gene.count",
              help="Input reads count file; such as OTU table [default %default]"),
  make_option(c("-t", "--taxonomy"), type="character", default="result/taxonomy.txt",
              help="Input taxonomy file [default %default]"),
  make_option(c("-d", "--design"), type="character", default="result/metadata.tsv",
              help="Experiment design or sample metadata [default %default]"),
  make_option(c("-T", "--threshold"), type="numeric", default=0.1,
              help="Threshold of abundance percentage, such as 0.1% [default %default]"),
  make_option(c("-n", "--group"), type="character", default="Group",
              help="Column name of group [default %default]"),
  make_option(c("-o", "--output"), type="character", default="result/otutab.txt",
              help="Output directory and filename [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

# 显示输入输出确认是否正确
print(paste("The input feature table is ", opts$input,  sep = ""))
# print(paste("Input taxonomy file: ", opts$taxonomy,  sep = ""))
# print(paste("Experiment design or sample metadata: ", opts$design, sep = ""))
# print(paste("Threshold of abundance percentage: ", opts$threshold,  sep = ""))
# print(paste("Column name of group: ", opts$group,  sep = ""))
print(paste("Output directory and filename: ", opts$output,  sep = ""))


# 1. 读取OTU表
otutab = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)

# 2. 取整
otutab = otutab + 0.5
otutab = round(otutab)

# 3. 保存
write.table("OTUID\t", file=opts$output, append = F, sep="\t", quote=F, eol = "", row.names=F, col.names=F)
suppressWarnings(write.table(otutab, file=opts$output, append = T, sep="\t", quote=F, row.names=T, col.names=T))
