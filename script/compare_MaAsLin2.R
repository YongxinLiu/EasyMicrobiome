#!/usr/bin/env Rscript

# Copyright 2024-2026 Defeng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Bai, et al. 2025. EasyMetagenome: A User‐Friendly and Flexible Pipeline for Shotgun Metagenomic Analysis in Microbiome Research. iMeta 4: e70001. https://doi.org/10.1002/imt2.70001

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 更新
# 2024/11/12：增加MaAsLin2方法用于组间差异比较
# 2025/11/27：规范代码


# Clean enviroment object
rm(list=ls()) 

# 1.1 程序功能描述和主要步骤

# 程序功能：运用MaAsLin2方法计算组间差异
# Functions: Difference analysis using MaAsLin2

options(warn = -1) # Turn off warning


## 设置输入输出文件和参数

# 修改下面`default=`后面的文件和参数。
#
# 输入文件为物种相对丰度表(Species.txt)+分组信息(metadata.txt)
#
# 输入文件"-i", "--input"，metaphlan4/Species.txt; 物种相对丰度表格；
#
# 实验设计"-d", "--design"，默认`metadata.txt`，可手动修改文件位置；
#
# 丰度筛选"-a", "--min_abundance, 默认为'0';
#
# 流行率筛选"-p", "--min_prevalence", 默认为'0.1';
#
# 显著性阈值"-s", "--max_significance", 默认为'0.05'
#
# 标准化方法"-n", "--normalization", 默认为'TSS', 可选TSS, CLR, CSS, NONE, TMM；
#
# 数据转换方法"-t", "--transform", 默认为'Log', 可选LOG, LOGIT, AST, NONE；
#
# 分析方法"-m", "--analysis_method", 默认为'LM', 可选LM, CPLM, NEGBIN, ZINB；
#
# 随机效应"-r", "--random_effects", 默认为"none"; 
#
# 固定效应"-f", "--fixed_effects", 默认为"Group";
#
# p值矫正方法"c", "--correction", 默认为"BH"; 
#
# 分组列名"-o", "--output"，默认为输出目录，输出差异分析结果文件;


# 1.2 依赖包安装

site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
a = rownames(installed.packages())

# install CRAN
install_CRAN <- c("optparse", "dplyr", "reshape2",  "readxl", "tibble",  "openxlsx",
                  "foreach", "data.table",  "gridExtra", "scales", "ggplot2",  "ggh4x",
                  "ggfortify", "ggvenn",  "ggrepel", "vegan", "pairwiseCI",  "vcd",
                  "igraph", "sampling", "CVXR", "DescTools", "phyloseq")
for (i in install_CRAN) {
  if (!i %in% a)
    install.packages(i, repos = site)
}

# install bioconductor
install_bioc <- c( "phyloseq", "ANCOMBC", "Maaslin2")
for (i in install_bioc) {
  if (!i %in% a)
    BiocManager::install(i, update = F) # , site_repository=site
  a = rownames(installed.packages())
}

# install github
if (!"amplicon" %in% a){
  devtools::install_github("microbiota/amplicon")
}


# 1.3 解析命令行
# 解析参数-h显示帮助信息
suppressMessages(
  library("optparse")
)
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="metaphlan4/Species.txt",
                help="Unfiltered OTU table [default %default]"),
    make_option(c("-g", "--metadata"), type="character", default="metadata.txt",
                help="metadata file or metadata [default %default]"),
    make_option(c("-a", "--min_abundance"), type="numeric", default="0",
                help="Minimum abundance used for analysis [default %default]"),
    make_option(c("-p", "--min_prevalence"), type="numeric", default="0.1",
                help="Minimum prevalence used for analysis [default %default]"),
    make_option(c("-s", "--max_significance"), type="numeric", default="0.05",
                help="Significance threshold used for analysis [default %default]"),
    make_option(c("-n", "--normalization"), type="character", default="TSS",
                help="Normalization method used for analysis, TSS, CLR, CSS, NONE, TMM [default %default]"),
    make_option(c("-t", "--transform"), type="character", default="Log",
                help="Transformation method used for analysis, LOG, LOGIT, AST, NONE [default %default]"),
    make_option(c("-m", "--analysis_method"), type="character", default="LM",
                help="Analysis method used for analysis, LM, CPLM, NEGBIN, ZINB [default %default]"),
    make_option(c("-r", "--random_effects"), type="character", default="none",
                help="Random effects [default %default]"),
    make_option(c("-f", "--fixed_effects"), type="character", default="Group",
                help="Fixed effects [default %default]"),
    make_option(c("-c", "--correction"), type="character", default="BH",
                help="Fixed effects [default %default]"),
    make_option(c("-o", "--output"), type="character", default="metaphlan4/",
                help="Output quantile value for filter feature table [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
# print("You are using the following parameters:")
# print(opts)


# 2. 依赖关系检查、安装和加载

suppress <- function(x){invisible(capture.output(suppressMessages(suppressWarnings(x))))}

# 依赖包列表
package_list <- c(
  "optparse","dplyr", "reshape2",  "readxl", "tibble",  "openxlsx",
  "foreach", "data.table",  "gridExtra", "scales", "ggplot2",  "ggh4x",
  "ggfortify", "ggvenn",  "ggrepel", "vegan", "pairwiseCI",  "vcd",
  "igraph", "sampling", "CVXR", "DescTools", "phyloseq"
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


# # Install related packages
# # 基于CRAN安装R包，检测没有则安装 Installing R packages based on CRAN and installing them if they are not detected
# p_list = c("dplyr", "reshape2",  "readxl", "tibble",  "openxlsx",
#            "foreach", "data.table",  "gridExtra", "scales", "ggplot2",  "ggh4x",
#            "ggfortify", "ggvenn",  "ggrepel", "vegan", "pairwiseCI",  "vcd",
#             "igraph", "sampling", "CVXR", "DescTools")
# for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
#   library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}
# 
# #### LOAD REQUIRED R PACKAGES ####
# 
# # options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager", repos = site)
# a = rownames(installed.packages())
# 
# # install CRAN
# install_CRAN <- c("ggplot2", "reshape2", "readxl", "tibble","openxlsx", "foreach", 
#                   "data.table", "gridExtra","scales", "ggh4x", "ggfortify", "ggvenn",
#                   "ggrepel", "vegan", "pairwiseCI", "vcd", "igraph")
# for (i in install_CRAN) {
#   if (!i %in% a)
#   install.packages(i, repos = site)
# }
# 
# # install bioconductor
# install_bioc <- c( "phyloseq", "ANCOMBC", "Maaslin2")
# for (i in install_bioc) {
#   if (!i %in% a)
#     BiocManager::install(i, update = F) # , site_repository=site
#     a = rownames(installed.packages())
# }
# 
# # install github
# if (!"amplicon" %in% a){
#   devtools::install_github("microbiota/amplicon")
# }

# suppress <- function(x){invisible(capture.output(suppressMessages(suppressWarnings(x))))}
# suppress(library(dplyr))
# suppress(library(reshape2))
# suppress(library(readxl))
# suppress(library(phyloseq))
# suppress(library(tibble))
# suppress(library(openxlsx))
# suppress(library(foreach))
# suppress(library(data.table))
# suppress(library(gridExtra))
# suppress(library(scales))
# suppress(library(ggplot2))
# suppress(library(ggh4x))
# suppress(library(ggfortify))
# suppress(library(ggvenn))
# suppress(library(ggrepel))
# suppress(library(vegan))
# suppress(library(pairwiseCI))
# suppress(library(vcd))
# # library(CVXR)
# suppress(library(ANCOMBC))
# suppress(library(Maaslin2))
# suppress(library(igraph))


# 3. 读取输入文件

# 读取实验设计
# PREPARE RELATIVE ABUNDANCE AND COUNT DATA 
metadata <- read.delim(opts$metadata, header=T, sep="\t") 
rownames(metadata) <- metadata$SampleID

# 读取物种相对丰度表格
# read in tables that were previously generated by taxonomic profiling
ra <- read.delim(opts$input, header=T, sep="\t")
ra <- ra[,c('Taxonomy',metadata$SampleID)]

# make table sample x feature
rownames(ra) <- ra$Taxonomy
ra <- data.frame(t(ra[,-1]), check.names=FALSE)
#ra <- ra[, c(-1,-2)]

# compile relative abundance data into phyloseq objects for species
ra.sub <- ra
ra.ps.s <- phyloseq(otu_table(as.matrix(ra.sub), taxa_are_rows=FALSE),
                    sample_data(metadata),
                    tax_table(as.matrix(
                      data.frame(Taxonomy = colnames(ra.sub),
                                 check.names=FALSE, row.names=colnames(ra.sub)))))

# 获取 Group 列中不重复的分组
groups <- unique(sample_data(ra.ps.s)$Group)

# 使用 recode 将分组自动映射为 1, 2, ...
sample_data(ra.ps.s)$Group <- dplyr::recode(
  sample_data(ra.ps.s)$Group,
  !!!setNames(seq_along(groups), groups)
)

# lm.s.npc
library(phyloseq)
library(Maaslin2)
ci <- function(coef, se){
  lower.ci <- coef - 1.96*se
  upper.ci <- coef + 1.96*se
  return(c(lower.ci=lower.ci,upper.ci=upper.ci))
}

ps = phyloseq(otu_table(ra.ps.s)/100, sample_data(ra.ps.s))
#ps = phyloseq(otu_table(ra.ps.s), sample_data(ra.ps.s))
# run MaAsLin2
input_data <- data.frame(otu_table(ps))
library(sampling)
sam_data = as.data.frame(sample_data(ps))
common_columns = colnames(sam_data)[colnames(sam_data) %in% colnames(metadata)]
input_metadata <- data.frame(sam_data[, common_columns])
# capt<- capture.output(fits <- suppressWarnings(Maaslin2(input_data, input_metadata, 
#                                                         output='temp_directory', 
#                                                         min_prevalence=0.05, 
#                                                         #min_variance,
#                                                         #normalization='NONE', 
#                                                         normalization='CLR',
#                                                         #normalization= 'TSS',
#                                                         #transform='LOG',
#                                                         transform = "NONE",
#                                                         correction = "BH",
#                                                         analysis_method='LM',
#                                                         max_significance=0.05,
#                                                         fixed_effects = c('Group'),
#                                                         #standardize=FALSE,
#                                                         standardize=TRUE,
#                                                         #standardize=FALSE,
#                                                         plot_heatmap=TRUE, 
#                                                         plot_scatter=TRUE
# )))

capt<- capture.output(fits <- suppressWarnings(Maaslin2(input_data, input_metadata, 
                                                        output='temp_directory', 
                                                        min_prevalence=opts$min_prevalence, 
                                                        #min_variance,
                                                        #normalization='NONE', 
                                                        normalization=opts$normalization,
                                                        #normalization= 'TSS',
                                                        #transform='LOG',
                                                        transform = opts$transform,
                                                        correction = opts$correction,
                                                        analysis_method = opts$analysis_method,
                                                        max_significance= opts$max_significance,
                                                        fixed_effects = c(opts$fixed_effects),
                                                        #standardize=FALSE,
                                                        standardize=TRUE,
                                                        #standardize=FALSE,
                                                        plot_heatmap=TRUE, 
                                                        plot_scatter=TRUE
)))

# put back original feature names
for (feat in seq_along(fits$results$feature)){
  fits$results$feature[feat] <- taxa_names(ps)[make.names(taxa_names(ps)) == 
                                                 fits$results$feature[feat]]
}

input_metadata = input_metadata[!is.na(input_metadata$Case_status),]
input_metadata = input_metadata[!is.na(input_metadata$Sex),]
input_metadata = input_metadata[!is.na(input_metadata$Age),]
input_metadata = input_metadata[!is.na(input_metadata$smoke),]
input_metadata = input_metadata[!is.na(input_metadata$drink),]
input_metadata = input_metadata[!is.na(input_metadata$BMI),]

ps = phyloseq(otu_table(ps),sample_data(ra.ps.s))
sample_data(ps) = data.frame(sample_data(ps))

res <- data.frame()
for (var in seq_along(unique(fits$results$metadata))){
  # get variable name
  var.name <- unique(fits$results$metadata)[var]
  if (length(table(sample_data(ps)[,var.name])) == 2){
    group.1.index <- sample_data(ps)[,var.name] == 
      names(table(sample_data(ps)[,var.name]))[2]
    group.1.index[is.na(group.1.index)] <- FALSE
    group.2.index <- sample_data(ps)[,var.name] == 
      names(table(sample_data(ps)[,var.name]))[1]
    group.2.index[is.na(group.2.index)] <- FALSE
    n1 <- colSums(otu_table(ps)[group.1.index,] > 0)
    n2 <- colSums(otu_table(ps)[group.2.index,] > 0)
    mean1 <- colMeans(otu_table(ps)[group.1.index,])
    mean2 <- colMeans(otu_table(ps)[group.2.index,])
  }else{
    n1 <- rep(sum(table(sample_data(ps)[,var.name])), ntaxa(ps))
    names(n1) <- taxa_names(ps)
    n2 <- rep(NA, ntaxa(ps))
    names(n2) <- taxa_names(ps)
    mean1 <- colMeans(otu_table(ps))
    mean2 <- rep(NA, ntaxa(ps))
    names(mean2) <- taxa_names(ps)
  }
  # calculate fold change and confidence interval of fold change
  if(length(table(sample_data(ps)[,var.name])) == 2){
    FC <- 2^(fits$results$coef[fits$results$metadata == var.name])
    FC.lower <- c()
    FC.upper <- c()
    for (coef in seq_along(fits$results$coef[fits$results$metadata == var.name])){
      FC.lower <- c(FC.lower, 2^(ci(fits$results$coef[fits$results$metadata == 
                                                        var.name][coef],
                                    fits$results$stderr[fits$results$metadata == 
                                                          var.name][coef])['lower.ci']))
      FC.upper <- c(FC.upper, 2^(ci(fits$results$coef[fits$results$metadata == 
                                                        var.name][coef],
                                    fits$results$stderr[fits$results$metadata ==
                                                          var.name][coef])['upper.ci']))
    }
  }else{
    FC <- NA
    FC.lower <- NA
    FC.upper <- NA
  }
  # summarize results for variable
  correction = "BH"
  rvar <- data.frame(Variable=var.name,
                     Feature=fits$results$feature[fits$results$metadata == var.name],
                     N1=n1[fits$results$feature[fits$results$metadata == var.name]],
                     N2=n2[fits$results$feature[fits$results$metadata == var.name]],
                     Mean1=mean1[fits$results$feature[fits$results$metadata == var.name]],
                     Mean2=mean2[fits$results$feature[fits$results$metadata == var.name]],
                     Beta=fits$results$coef[fits$results$metadata == var.name],
                     SE=fits$results$stderr[fits$results$metadata == var.name],
                     P=fits$results$pval[fits$results$metadata == var.name],
                     FDR=p.adjust(fits$results$pval[fits$results$metadata == var.name], 
                                  method=correction),
                     FC=FC, FC_lower=FC.lower, FC_upper=FC.upper,
                     check.names=FALSE)
  res <- rbind(res, rvar[order(rvar$P),])
  # add untested features if they exist
  if (nrow(rvar) != ntaxa(ps)){
    res <- rbind(res,
                 data.frame(Variable=var.name,
                            Feature=taxa_names(ps)[!(taxa_names(ps) %in% 
                                                       fits$results$feature[fits$results$metadata == var.name])],
                            N1=n1[taxa_names(ps)[!(taxa_names(ps) %in% 
                                                     fits$results$feature[fits$results$metadata == var.name])]],
                            N2=n2[taxa_names(ps)[!(taxa_names(ps) %in% 
                                                     fits$results$feature[fits$results$metadata == var.name])]],
                            Mean1=mean1[taxa_names(ps)[!(taxa_names(ps) %in% 
                                                           fits$results$feature[fits$results$metadata == var.name])]],
                            Mean2=mean2[taxa_names(ps)[!(taxa_names(ps) %in% 
                                                           fits$results$feature[fits$results$metadata == var.name])]],
                            Beta=NA, SE=NA, P=NA, FDR=NA, FC=NA, FC_lower=NA, FC_upper=NA,
                            check.names=FALSE)
    )
  }
}
lm.s.npc = list(result.summary=res, Maaslin2.output=fits)

write.csv(lm.s.npc[["result.summary"]], paste0(opts$output, 'MaAsLin2_overall_difference.csv'))


# Enriched or Depleted
# get FDR q-values ready for plotting
plot.data <- lm.s.npc$result.summary[lm.s.npc$result.summary$Variable == 'Group',
                                     c('Feature','P','FDR','FC')]
plot.data <- plot.data[rowSums(is.na(plot.data)) == 0,]

plot.data$`NPC association` <- ifelse(plot.data[,2] < 0.05,
                                      ifelse(plot.data[,4] > 1, 'depleted',
                                             ifelse(plot.data[,4] < 1, 
                                                    'enriched','opposite directions')),
                                      'not associated')

write.csv(plot.data, paste0(opts$output, 'MaAsLin2_enriched_depleted.csv'))

