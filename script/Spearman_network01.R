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
  make_option(c("-i", "--input"), type="character", default="result/tax/otutab_amplicon.txt",
              help="OTU/ASV table [default %default]"),
  make_option(c("-g", "--group"), type = "character", default = "result/tax/sample_amplicon.txt", 
              help = "Group information [default %default]"),
  make_option(c("-t", "--tax"), type = "character", default = "result/tax/taxonomy_amplicon.txt", 
              help = "Taxonomy information [default %default]"),
  make_option(c("-o", "--output"), type="character", default="result/tax/",
              help="Output directory [default %default]"),
  make_option(c("-w", "--width"), type="numeric", default=120 * 1.5,
              help="Figure width (mm) [default %default]"),
  make_option(c("-e", "--height"), type="numeric", default=80 * 1.5,
              help="Figure height (mm) [default %default]")
)

# 解析命令行
opts = parse_args(OptionParser(option_list=option_list))

if (opts$output == "") { opts$output = opts$input }

#install.packages("BiocManager")
#library(BiocManager)
#install("remotes")
#install("tidyverse")
#install("tidyfst")
#install("igraph")
#install("sna")
#install("phyloseq")
#install("ggalluvial")
#install("ggraph")
#install("WGCNA")
#install("ggnewscale")
#install("pulsar")
#install("patchwork")
#remotes::install_github("taowenmicro/EasyStat")
#remotes::install_github("taowenmicro/ggClusterNet")

# 基于CRAN安装R包，检测没有则安装
p_list = c("BiocManager","MicrobiotaProcess","dplyr","ggplot2","phyloseq","ggtree",
           "ggtreeExtra","ggstar","forcats","conflicted")
for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
  library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

# 加载R包 Load the package
suppressWarnings(suppressMessages(library("BiocManager")))
suppressWarnings(suppressMessages(library("MicrobiotaProcess")))
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("ggplot2")))
suppressWarnings(suppressMessages(library("phyloseq")))
suppressWarnings(suppressMessages(library("ggtree")))
suppressWarnings(suppressMessages(library("ggtreeExtra")))
suppressWarnings(suppressMessages(library("ggstar")))
suppressWarnings(suppressMessages(library("forcats")))
suppressWarnings(suppressMessages(library("conflicted")))
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_scout()
conflicted::conflicts_prefer(phyloseq::tax_table)

# 导入数据
# load data
#OTU<- read.table("Figure3/data/otutab_amplicon.txt",check.names = F, row.names = 1, header = 1, sep = "\t")
OTU<- read.table(opts$input, check.names = F, row.names = 1, header = 1, sep = "\t")
#sample <- read.table("Figure3/data/sample_amplicon.txt",check.names = F, row.names = 1, header = 1, sep = "\t")
sample <- read.table(opts$group,check.names = F, row.names = 1, header = 1, sep = "\t")
#Tax <- read.table("Figure3/data/taxonomy_amplicon.txt",check.names = F, row.names = 1, header = 1, sep = "\t")
Tax <- read.table(opts$tax, check.names = F, row.names = 1, header = 1, sep = "\t")

# 利用phyloseq包重新构造可转换为分析的数据格式
# Reconstructing data formats that can be converted into analysis using the photoseq package
ps <- phyloseq(sample_data(sample),
               otu_table(as.matrix(OTU), taxa_are_rows=TRUE), 
               tax_table(as.matrix(Tax)))

library(ggClusterNet)
library(phyloseq)
library(tidyverse)
#data(ps)
otupath = "./result/"
#  新网络分析-2023年末更新#--------

netpath = paste(otupath,"/network.new/",sep = "")
dir.create(netpath)

rank.names(ps)
library(ggrepel)
library(igraph)
#detach("package:MicrobiotaProcess")

# 8.1 网络分析主函数#--------
tab.r = network.pip(
  ps = ps,
  N = 500,
  # ra = 0.05,
  big = TRUE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.6,
  p.threshold = 0.05,
  maxnode = 2,
  method = "spearman",
  label = TRUE,
  lab = "elements",
  group = "group",
  fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R=10,
  ncpus = 4
)

#  建议保存一下输出结果为R对象，方便之后不进行相关矩阵的运算，节约时间
saveRDS(tab.r,paste0(netpath,"network.pip.sparcc.rds"))
tab.r = readRDS(paste0(netpath,"network.pip.sparcc.rds"))

dat = tab.r[[2]]
cortab = dat$net.cor.matrix$cortab
# 大型相关矩阵跑出来不容易，建议保存，方便各种网络性质的计算
saveRDS(cortab,paste0(netpath,"cor.matrix.all.group.rds"))
cor = readRDS(paste0(netpath,"cor.matrix.all.group.rds"))

#-提取全部图片的存储对象
plot = tab.r[[1]]
# 提取网络图可视化结果
p0 = plot[[1]]
ggsave(paste0(opts$output, "Spearman_network01.pdf"), plot = p0, width = 10, height = 3, dpi = 600)

# 可以对图进行调整，使其变得更美观
dat = tab.r[[2]]
node = dat$net.cor.matrix$node
edge = dat$net.cor.matrix$edge
head(edge)
head(node)

node2  = add.id.facet(node,"group")
head(node2)

p <- ggplot() + geom_segment(
  aes(
    x = X1,
    y = Y1,
    xend = X2,
    yend = Y2,
    color = cor
  ),
  data = edge,
  size = 0.2,
  alpha = 0.7
) +
  geom_point(
    aes(X1, X2, fill = Phylum, size = igraph.degree),
    pch = 21,
    data = node,
    #color = "gray40"
    color = "white"
    #color = NA
  ) +
  facet_wrap(. ~ label, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = c("#66c2a4", "#BC80BD")) +
  # scale_fill_hue()+
  scale_size(range = c(0.6, 3.5)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white", colour = NA)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
p

ggsave(paste0(opts$output, "Spearman_network02.pdf"), plot = p, width = 14, height = 3.5, dpi = 600)

