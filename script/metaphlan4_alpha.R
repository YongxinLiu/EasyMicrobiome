#!/usr/bin/env Rscript

# Copyright 2024-2026 Defeng Bai <baidefeng@caas.cn> 

# If used this script, please cited:
# Bai, et al. 2025. EasyMetagenome: A User‐Friendly and Flexible Pipeline for Shotgun Metagenomic Analysis in Microbiome Research. iMeta 4: e70001. https://doi.org/10.1002/imt2.70001

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 更新
# 2024/11/12：增加利用Metaphlan4相对丰度表计算Alpha的功能
# 2025/11/25：规范代码

# 1.1 程序功能描述和主要步骤

# 程序功能：Metaphlan4相对丰度表计算Alpha多样性
# Functions: Get Alpha indices using Metaphlan4 relative abundance table

options(warn = -1) # Turn off warning


## 设置输入输出文件和参数

# 修改下面`default=`后面的文件和参数。
#
# 输入文件为原始alpha多样性(metaphlan4/taxonomy.tsv)+分组信息(metadata.txt)
#
# 输入文件"-i", "--input"，metaphlan4/taxonomy.tsv; 微生物相对丰度表格
#
# 实验设计"-g", "--metadata"，默认`metadata.txt`，可手动修改文件位置；
#
# 分组列名"-t", "--taxonomy"，默认'7'，代表物种级别；
#
# 分组列名"-o", "--output"，默认为输出目录；


# 1.2 依赖包安装

site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
a = rownames(installed.packages())

# install CRAN
install_CRAN <- c("optparse","dplyr","tidyr","reshape2","vegan")
for (i in install_CRAN) {
  if (!i %in% a)
    install.packages(i, repos = site)
  require(i,character.only=T)
  a = rownames(installed.packages())
}

# install bioconductor
install_bioc <- c("ggplot2", "multcompView")
for (i in install_bioc) {
  if (!i %in% a)
    BiocManager::install(i, update = F)
  a = rownames(installed.packages())
}

# install github
if (!"amplicon" %in% a){
  devtools::install_github("microbiota/amplicon")
}


# 1.3 解析命令行
# 解析参数-h显示帮助细腻些
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="metaphlan4/taxonomy.tsv", help="Metaphlan4 relative abundance [default %default]"),
    make_option(c("-g", "--metadata"), type="character", default="metadata.txt", help="Metaphlan4 [default %default]"),
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


# 2. 依赖关系检查、安装和加载

# 依赖包列表
package_list <- c(
  "tidyr","reshape2","vegan"
)

# 批量安装和加载
for (p in package_list) {
  # 如果未安装，则安装
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  }
  
  # 批量加载，抑制警告和消息
  suppressWarnings(
    suppressMessages(
      library(p, character.only = TRUE)
    )
  )
}


# 3. 读取输入文件

# 读取metaphlan4 taxonomy.tsv文件
# 默认的quote会跳过2/3的数据，导致行减少产生NA，改默认值为空
taxonomy = read.table(opts$input, header=T, sep="\t", quote = "", row.names=NULL, comment.char="")
print(paste0("All taxonomy annotations are ", dim(taxonomy)[1], " lines!"))

# 去除NA，否则无法计算
taxonomy = na.omit(taxonomy)

# 显示样本总数据，有冗余
i = opts$taxonomy

# 读取实验设计
metadata = read.table(opts$metadata, header=T, sep="\t", row.names=1, quote = "", comment.char="")


# 4. Alpha多样性指数计算
otutable <- separate(taxonomy, ID, c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Taxonomy"),sep="\\|",extra = "drop", fill = "right")
otutable <- otutable[-which(is.na(otutable$Taxonomy)),]

otutable = data.frame(otutable,stringsAsFactors = F) 
otutable[,9:ncol(otutable)] = as.data.frame(lapply(otutable[,9:ncol(otutable)],as.numeric))
colnames(metadata)[1] <- 'Group'

# 假设 otutable 已经处理好，i 是你想用的列索引（例如 Genus）
# 构建长格式数据
level <- otutable %>%
  select(!!sym(colnames(otutable)[i]), 9:ncol(otutable)) %>%
  pivot_longer(
    cols = 2:ncol(.),
    names_to = "sample",
    values_to = "relative_abundance"
  )

# 转换回宽格式
level <- level %>%
  group_by(!!sym(colnames(otutable)[i]), sample) %>%
  summarise(relative_abundance = sum(relative_abundance), .groups = "drop") %>%
  pivot_wider(
    names_from = sample,
    values_from = relative_abundance,
    values_fill = 0
  )

# reshape2中函数使用有问题
# level = cbind(otutable[,i],otutable[,9:ncol(otutable)])
# 
# level = melt(level,id.vars= colnames(level)[1],
#              measure.vars = colnames(level[,2:ncol(level)]),
#              variable.name = "sample",value.name = "relative_abundance")
# 
# level = dcast(level, otutable[, i] ~ sample, fun.aggregate = sum)

# RA transposition
level = t(level)
colnames(level) = level[1,]
level = level[-1,]
level = data.frame(level,stringsAsFactors = F) 
colnames(level) = gsub('\\.','',colnames(level))
level1 = apply(level,2,as.numeric)
rownames(level1)=rownames(level)
level=level1

# Calculate diversity
level_diversity = data.frame(Sample_ID = colnames(otutable[9:ncol(otutable)]),
                             observed_species=specnumber(level),
                             shannon=vegan::diversity(level, index="shannon"),
                             simpson=vegan::diversity(level, index="simpson"),
                             invsimpson=vegan::diversity(level, index="invsimpson"),
                             Pielou_evenness=vegan::diversity(level,
                                                              index="shannon")/log(specnumber(level)))

# Export diversity table
write.table(level_diversity, file=paste(opts$output, ".txt", sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
