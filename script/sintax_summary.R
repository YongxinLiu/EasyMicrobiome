#!/usr/bin/env Rscript

# Copyright 2016-2020 Tong Chen <chentong_biology@163.com>

# If used this script, please cited:
# Tong Chen, Yong-Xin Liu, Luqi Huang. 2022. ImageGP: An easy-to-use data visualization web server for scientific researchers. iMeta 1: e5. https://doi.org/10.1002/imt2.5
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, Yang Bai. 2021. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12: 315-330. https://doi.org/10.1007/s13238-020-00724-8


# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

#----1. 参数 Parameters#----

#----1.1 功能描述 Function description#----

# 程序功能：物种组成弦图
# Functions: Taxonomy circlize

options(warn = -1) # Turn off warning


## 设置输入输出文件和参数

# 修改下面`default=`后面的文件和参数。
#
# 输入文件为特征表(otutab.txt)+分组信息(metadata.tsv)
#
# 输入文件"-i", "--input"，otutab.txt; 特征表
#
# 实验设计"-d", "--design"，默认`metadata.tsv`，可手动修改文件位置；
#
# 分组列名"-n", "--group"，默认将metadata.tsv中的Group列作为分组信息，可修改为任意列名；
#
# 分组列名"-c", "--compare_pair"，默认将比较metadata.tsv中的Group列的前两个值，建议手动设定；
#
# 分组列名"-t", "--threhold"，丰度筛选阈值，默认千分之1
#
# 图片宽"-w", "--width"，默认89 mm，根据图像布局可适当增大或缩小
#
# 图片高"-e", "--height"，默认59 mm，根据图像布局可适当增大或缩小


#----1.2 参数缺少值 Default values#----
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-t", "--taxonomy"), type="character", default="result/taxonomy.txt",
                help="Taxonomy table [default %default]"),
    make_option(c("-i", "--otutab"), type="character", default="result/otutab_rare.txt",
                help="Feature table [default %default]"),
    make_option(c("-r", "--rank"), type="character", default="p",
                help="Rank value like p c o f g [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/tax/sum_p.txt",
                help="Output ; [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
# 设置输出文件缺省值，如果为空，则为输入+pcoa.pdf

suppressWarnings(dir.create(dirname(opts$output), showWarnings = F))


#----1.3. 加载包 Load packages#----

suppressWarnings(suppressMessages(library(dplyr)))


#----2. 读取文件 Read files#----

taxonomy_file = opts$taxonomy
otutab_file   = opts$otutab
rank = toupper(opts$rank)


# taxonomy_file = "result/taxonomy.txt"
# otutab_file   = "result/otutab_rare.txt"

#----2.1 实验设计 Metadata#----
taxonomy = read.table(taxonomy_file, header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)

rank_column = match.arg(rank, colnames(taxonomy))

#----2.1 特征表 Feature table#----
otutab = read.table(otutab_file, header=T, row.names=1, sep="\t", comment.char="", quote="")

taxonomy_otutab = merge(taxonomy, otutab, by=0)

# head(taxonomy_otutab)

summary_abundance <- function(data, group){
  a1 = data %>% group_by(across(all_of(group))) %>% summarise_if(is.numeric, sum)
  a1 <- as.data.frame(a1)
  rownames(a1) = a1[[1]]
  a1 = a1[,-1]
  a_sum = apply(a1,2,sum)
  a1 = t(t(a1) / a_sum) * 100
}

rank_data = summary_abundance(taxonomy_otutab, rank_column)

rank_data = rank_data[order(apply(rank_data, 1, sum), decreasing = T) ,]

#----3.2 保存表格 Saving#----
filename = paste0(opts$output)
write.table("ID\t", file = filename,
            append = F, quote = F, eol = "", row.names = F, col.names = F)
suppressWarnings(write.table(rank_data, file = filename,
                             append = T, quote = F, sep = '\t', row.names = T, col.names = T))
