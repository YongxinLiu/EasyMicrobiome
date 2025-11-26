#!/usr/bin/env Rscript
# Copyright 2016-2021 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, Yang Bai. (2021). 
# A practical guide to amplicon and metagenomic analysis of microbiome data. Protein Cell 12 
# doi: https://doi.org/10.1007/s13238-020-00724-8


#---------0. Prepare #----

#---------0.1 Functions and examples #----

# Functions: Merge tables
# Main steps: 
# - Merge table by row names
# - Fill the NA as 0

# Examples
# merge2tableFill0.R -1 table1.txt -2 table2.txt -n KO -o merged_table.txt


#---------0.2 Parameters #----
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T) 
}

if (TRUE){
  option_list = list(
    make_option(c("-i", "--input1"), type="character", default="table1.txt",
                help="Input table [default %default]"),
    make_option(c("-j", "--input2"), type="character", default="table2.txt",
                help="Input table [default %default]"),
    make_option(c("-n", "--name"), type="character", default="ID",
                help="Colname of rownames [default %default]"),
    make_option(c("-o", "--output"), type="character", default="merged_table.txt",
                help="output table [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
print(paste0("The input table 1 ", opts$input1))
print(paste0("The input table 2 ", opts$input2))

#---------0.3 Packages #----
# package_list = c("dplyr","ggplot2","agricolae") 
# for(p in package_list){
#   if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
#     install.packages(p, repos=site)
#     suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
#   }
# }

#---------1 Read input #----
df1 = read.table(opts$input1, header=T, row.names=1, sep="\t", comment.char="") 
dim(df1)
df2 = read.table(opts$input2, header=T, row.names=1, sep="\t", comment.char="") 
dim(df2)

#---------2 Calculate #----
# Merge by row names, keep all ID
df = merge(df1, df2, by="row.names", all = T)
# Reset name
rownames(df) = df$Row.names
# Fill 0 for miss data
df[is.na(df)]=0
# Remove duplicate rownames
df = df[,-1]

#---------3 Write output #----
write.table(paste0(opts$name, "\t"), file = opts$output, append = F, quote = F, eol = "", row.names = F, col.names = F)
suppressWarnings(write.table(round(df,3), file=opts$output, append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))
print(paste("Merge finished! Output file ", opts$output, sep = ""))
# Summary output
dim(df)
colSums(df)[1:3]

#---------S1. Test Codes #----
# Read 2 tables as dataframe
# df1 = read.table("/mnt/m1/yongxin/medicago/lyr4_201230/result/eggnog/sum.KO.raw.txt", header=T, row.names=1, sep="\t", comment.char="") 
# df2 = read.table("/mnt/m1/yongxin/medicago/lyr4_210115b12/result/eggnog/sum.KO.raw.txt", header=T, row.names=1, sep="\t", comment.char="") 
# 
# # Merge by row names, keep all ID
# df = merge(df1, df2, by="row.names", all = T)
# # Reset name
# rownames(df) = df$Row.names
# # Fill 0 for miss data
# df[is.na(df)]=0
# # Remove duplicate rownames
# df = df[,-1]
# # Check point: summary each sample
# colSums(df)
# 
# # Save table
# write.table("ID\t", file=paste("/mnt/m1/yongxin/medicago/medAll/result/eggnog/sum.KO.txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
# suppressWarnings(write.table(round(df,3), file=paste("/mnt/m1/yongxin/medicago/medAll/result/eggnog/sum.KO.txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

