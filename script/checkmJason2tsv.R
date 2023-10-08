#!/usr/bin/env Rscript
# 
# Copyright 2016-2022 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo & Yang Bai. (2021). 
# A practical guide to amplicon and metagenomic analysis of microbiome data. 
# Protein & Cell 12, 315-330, doi: https://doi.org/10.1007/s13238-020-00724-8

# Modified from https://git.list.lu/eScience/ICoVeR/blob/8d530cc1db29d71ed0c433e95cbd7ca92f036848/R.postprocessing/ParseCheckMResults.R

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 1.1 程序功能描述和主要步骤

# 程序功能：格式化jason格式的checkM结果为表格
# Functions: Fomrat checkm jason to table seperated value

# 程序使用示例
# # 显示脚本帮助
# Rscript ~/db/script/checkmJason2tsv.R -h
# Rscript ~/db/script/checkmJason2tsv.R -i temp/checkM/storage/bin_stats_ext.tsv -o temp/checkM/storage/bin_stats_ext.txt

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
  make_option(c("-i", "--input"), type="character", default="temp/checkM/storage/bin_stats_ext.tsv",
              help="OTU [default %default]"),
  make_option(c("-o", "--output"), type="character", default="",
              help="output directory or prefix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))
# 自动设置默认值
if(opts$output==""){opts$output=gsub(".tsv", ".txt", opts$input)}
# 显示输入输出确认是否正确
print(paste("Input file: ", opts$input,  sep = ""))
print(paste("Output file: ", opts$output, sep = ""))


checkm.raw <- read.csv(opts$input, sep="\t", header = F)
checkm.raw[,2] <- gsub("'", "\"", checkm.raw[,2])
checkm.parsed <- as.data.frame(t(rbind(sapply(checkm.raw[,2], jsonlite::fromJSON))))
rownames(checkm.parsed) <- checkm.raw$V1

# IF checkm.file == bin_stats_ext THEN
ignore <- c("GCN0", "GCN1", "GCN2", "GCN3", "GCN4", "GCN5+")
checkm.parsed <- checkm.parsed[ , -which(names(checkm.parsed) %in% ignore)]
checkm.parsed[ , c(1:4, 6:27)] <- apply(checkm.parsed[ , c(1:4, 6:27)], 2, function(x) as.numeric(x))

# ELSE IF checkm.file == bin_stats.analyze.tsv THEN
# checkm.parsed[ , c(1:15)] <- apply(checkm.parsed[ , c(1:15)], 2, function(x) as.numeric(x))

# 筛选重要结果
colnames(checkm.parsed)
checkm.parsed = checkm.parsed[,c("Genome size", "Longest contig", "GC", "Completeness", "# contigs", "Coding density", "Mean contig length", "N50 (contigs)", "Contamination", "# predicted genes")]
idx = order(checkm.parsed$Contamination, decreasing = T)
checkm.parsed = checkm.parsed[idx,]
# Write the results as csv
# Error in write.table(checkm.parsed, file = gsub(".tsv", ".csv", checkm.file),  : 种类'list'目前没有在'EncodeElement'里实现
# write.csv(checkm.parsed, file = gsub(".tsv", ".csv", checkm.file))
write.table(paste("ID\t", sep=""), file=opts$output,append = F, quote = F, eol = "", row.names = F, col.names = F)
write.table(checkm.parsed, file=opts$output, append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)

# 统计污染>5的数量，53/586=9%污染>5%
idx = checkm.parsed$Contamination > 5
table(idx)

# Or view the table within R directly
# View(checkm.parsed)
