#!/usr/bin/env Rscript

# 关闭警告信息
options(warn = -1)

# 设置清华源
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"

# 解析命令行参数
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE)))) {
  install.packages("optparse", repos=site)
  require("optparse", character.only = TRUE)
}

# 设置命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type="character", default="result2/tax/otutab.txt",
              help="OTU composition [default %default]"),
  make_option(c("-g", "--metadata"), type = "character", default = "result2/tax/metadata.txt", 
              help = "Metadata"),
  make_option(c("-t", "--tax"), type = "character", default = "result2/tax/taxonomy.txt", 
              help = "Taxonomy"),
  make_option(c("-p", "--phylo"), type = "character", default = "result2/tax/otus.tree", 
              help = "Taxonomy"),
  make_option(c("-o", "--output"), type="character", default="result2/tax/",
              help="Output directory [default %default]"),
  make_option(c("-w", "--width"), type="numeric", default=120 * 1.5,
              help="Figure width (mm) [default %default]"),
  make_option(c("-e", "--height"), type="numeric", default=80 * 1.5,
              help="Figure height (mm) [default %default]")
)

# 解析命令行
opts = parse_args(OptionParser(option_list=option_list))

if (opts$output == "") { opts$output = opts$input }

# 检查并安装必要的R包
packages_needed <- c("microeco", "ape", "magrittr", "ggplot2")
packages_missing <- packages_needed[!(packages_needed %in% installed.packages()[,"Package"])]
if(length(packages_missing)) {
  install.packages(packages_missing)
}

# microeco一般需要从GitHub安装
if(!require("microeco")){
  if(!require("devtools")) install.packages("devtools")
  devtools::install_github("ChiLiubio/microeco")
}

suppressWarnings(suppressMessages(library(microeco)))
suppressWarnings(suppressMessages(library(ape)))
suppressWarnings(suppressMessages(library(magrittr)))
suppressWarnings(suppressMessages(library(ggplot2)))

# 设定主题
theme_set(theme_bw())

# 设置随机种子，保证可复现
set.seed(123)

# 读入数据
otu_table_16S <- read.delim(opts$input, row.names = 1, header = TRUE, sep = "\t")
otu_table_16S <- otu_table_16S[, 1:12]
sample_info_16S <- read.delim(opts$metadata, row.names = 1, header = TRUE, sep = "\t")
taxonomy_table_16S <- read.delim(opts$tax, row.names = 1, header = TRUE, sep = "\t")
phylo_tree_16S <- read.tree(opts$phylo)

# 处理分类表：统一前缀
taxonomy_table_16S %<>% tidy_taxonomy

# 创建microtable对象
mt <- microtable$new(sample_table = sample_info_16S,
                     otu_table = otu_table_16S,
                     tax_table = taxonomy_table_16S,
                     phylo_tree = phylo_tree_16S)
print(mt)

# 删除非细菌和古菌的OTU
mt$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
print(mt)

# 修剪数据集
mt$tidy_dataset()
print(mt)

# 过滤线粒体、叶绿体污染
mt$filter_pollution(taxa = c("mitochondria", "chloroplast"))
mt$tidy_dataset()
print(mt)

# 查看样本测序深度范围
print(mt$sample_sums() %>% range)

# 样本稀释标准化
mt$rarefy_samples(sample.size = 98)
print(mt)

# 保存处理后的基础数据表
mt$save_table(dirpath = "basic_files", sep = ",")

# 将环境因子加入样本表,因为metadata文件里有环境因子故省略
# mt$sample_table <- data.frame(mt$sample_table, env_data_16S[rownames(mt$sample_table), ])
# print(mt)

# 创建环境因子分析对象
t1 <- trans_env$new(dataset = mt, env_cols = 4:7)

# dbRDA分析，使用bray-curtis距离
t1$cal_ordination(method = "dbRDA", use_measure = "bray")

# 可视化（第一次）
t1$trans_ordination()
t1$plot_ordination(plot_color = "Group")

# 调整箭头长度，更美观
t1$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1.5)

# 重新绘制并保存图片
p11 <- t1$plot_ordination(plot_color = "Group")
print(p11)
ggsave(paste0(opts$output, "RDA_p11.pdf"), p11, width = 6, height = 4)

# 保存R环境
# save(mt, t1, file = "result2/tax/microeco_analysis_result.RData")

# 在 Genus 水平做 RDA 分析
t1$cal_ordination(method = "RDA", taxa_level = "Genus")

# 选择贡献度最高的10个分类单元（OTUs），并且调整箭头长度
t1$trans_ordination(
  show_taxa = 10,
  adjust_arrow_length = TRUE,
  max_perc_env = 1.5,
  max_perc_tax = 1.5,
  min_perc_env = 0.2,
  min_perc_tax = 0.2
)

# t1$res_rda_trans 是转换后的结果，可以直接用于画图
p12 <- t1$plot_ordination(plot_color = "Group")
p12

# 保存图片
ggsave(paste0(opts$output, "RDA_p12.pdf"), p12, width = 6, height = 4)

