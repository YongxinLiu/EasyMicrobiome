#!/usr/bin/env Rscript

# Copyright 2016-2024 De-feng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo & Yang Bai. (2021). A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12, doi: https://doi.org/10.1007/s13238-020-00724-8


# 1. 分析前准备：帮助、参数、依赖包和读取文件

# 命令行运行为当前目录；Rstudio手动运行脚本前需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录


# 1.1 程序功能描述和主要步骤

# 程序功能：Metaphlan4相对丰度表计算Alpha多样性
# 主要步骤：
# - 读取Metaphlan4物种组成表 result12/metaphlan4/taxonomy.tsv
# - 读取样本分组信息 result12/metadata.txt

# # 程序使用示例 USAGE
# # 显示脚本帮助 help
# Rscript /db/script/metaphlan4_alpha.R -h
# # 默认读取result12/metaphlan4/taxonomy.tsv
# # 完整参数：-i输入metaphlan4文件；-g输入metadata文件
# # -t 分类级别，可选Kingdom/Phylum/Class/Order/Family/Genus/Species，界门纲目科属种，推荐种
# # -o输出图表前缀，默认根据输入文件、分类级别自动生成

# Rscript /db/script/metaphlan4_alpha.R \
#   -i result12/metaphlan4/taxonomy.tsv \
#   -g result12/metadata.txt \
#   -t 7 \# Species
#   -o alpha_diversity.txt

options(warn = -1) # Turn off warning

# 1.2 解析命令行
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T) 
}


# 解析参数-h显示帮助细腻些
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="result12/metaphlan4/taxonomy.tsv", help="Metaphlan4 relative abundance [default %default]"),
    make_option(c("-g", "--metadata"), type="character", default="result12/metadata.txt", help="Metaphlan4 [default %default]"),
    make_option(c("-t", "--taxonomy"), type="numeric", default="7", help="Taxonomy level [default %default]"),
    make_option(c("-o", "--output"), type="character", default="", help="Output Metaphlan4 alpha diversity filename [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  
  # 调置如果无调设置输出，根据其它参数设置默认输出
  prefix = gsub("taxonomy.tsv$", "", opts$input, perl = T)
  if (opts$output==""){opts$output=paste0(prefix, "alpha")}
  
  # 显示输入输出确认是否正确
  # Metaphlan2物种组成表
  print(paste("The input file: ", opts$input,  sep = ""))
  # Metadata信息
  print(paste("The metadata file: ", opts$metadata,  sep = ""))
  # 绘制的分类级别, 默认为种
  print(paste("Taxonomy level: ", opts$taxonomy, ". Default if Species", sep = ""))
  # 输出文件名，不填则为输入目录+heatmap+taxonomy
  print(paste("Output Metaphlan4 alpha diversity filename: ", opts$output, sep = ""))
}


# 1.3 安装CRAN来源常用包
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("tidyr","reshape2","vegan")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}


# 1.4 读取输入文件

# 读取metaphlan4 taxonomy.tsv文件
# 默认的quote会跳过2/3的数据，导致行减少产生NA，改默认值为空
taxonomy = read.table(opts$input, header=T, sep="\t", quote = "", row.names=NULL, comment.char="")
print(paste0("All taxonomy annotations are ", dim(taxonomy)[1], " lines!"))
# 去除NA，否则无法计算
taxonomy = na.omit(taxonomy)
# 显示样本总数据，有冗余
# colSums(taxonomy)
i = opts$taxonomy

# 读取metaphlan4 metadata.txt文件
metadata = read.table(opts$metadata, header=T, sep="\t", row.names=1, quote = "", comment.char="")

# 2. 计算过程
# Results from metaphlan4
# otutable <- read.table("data/taxonomy.tsv",header=T,sep='\t',stringsAsFactors = F)
otutable <- separate(taxonomy, ID, c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Taxonomy"),sep="\\|",extra = "drop", fill = "right")
otutable <- otutable[-which(is.na(otutable$Taxonomy)),]

otutable = data.frame(otutable,stringsAsFactors = F) 
otutable[,9:ncol(otutable)] = as.data.frame(lapply(otutable[,9:ncol(otutable)],as.numeric))
colnames(metadata)[1] <- 'Group'

level = cbind(otutable[,i],otutable[,9:ncol(otutable)])
level = melt(level,id.vars= colnames(level)[1],
             measure.vars = colnames(level[,2:ncol(level)]),
             variable.name = "sample",value.name = "relative_abundance")
level = dcast(level, otutable[, i] ~ sample, fun.aggregate = sum)
#RA transposition
level = t(level)
colnames(level) = level[1,]
level = level[-1,]
level = data.frame(level,stringsAsFactors = F) 
colnames(level) = gsub('\\.','',colnames(level))
level1 = apply(level,2,as.numeric)
rownames(level1)=rownames(level)
level=level1

#Calculate diversity
level_diversity = data.frame(Sample_ID = colnames(otutable[9:ncol(otutable)]),
                             observed_species=specnumber(level),
                             shannon=vegan::diversity(level, index="shannon"),
                             simpson=vegan::diversity(level, index="simpson"),
                             invsimpson=vegan::diversity(level, index="invsimpson"),
                             Pielou_evenness=vegan::diversity(level,
                                                              index="shannon")/log(specnumber(level)))

#Export diversity table
write.table(level_diversity, file=paste(opts$output, ".txt", sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
