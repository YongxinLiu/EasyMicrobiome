#!/usr/bin/env Rscript

# Copyright 2026 Defeng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, et al. 2021. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12: 315-330. https://doi.org/10.1007/s13238-020-00724-8

# 手动运行脚本，使用 Ctrl+Shift+H 或 Session 需要设置工作目录
# Set Work Directory - Choose Directory / To Source File Location

# Clean enviroment object
rm(list=ls()) 

# 1.1 简介 Introduction #----

# 程序功能：SparCC可视化
# Functions: Sparcc Visualization


options(warn = -1) # Turn off warning

# 1.2 参数 Parameters #----
# 设置清华源加速下载
# (Optional) Set up Tsinghua Mirror to speed up download
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析包是否安装，没安装则安装，然后加载
# Determine whether the command line parsing package is installed, install it if it is not installed, then load
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-R", "--Correlation"), type="character", default="metaphlan4/R_Centenarians.txt",
                help="Metaphlan4 species table"),
    make_option(c("-P", "--Pvalue"), type="character", default="metaphlan4/P_Centenarians.txt",
                help="Unfiltered OTU table [default %default]"),
    make_option("--Group", type="character"),
    make_option(c("-r", "--output"), type="character", default="metaphlan4/",
                help="Output file for SparCC analysis in different groups  [default %default]") 
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
print("You are using the following parameters:")
print(opts)

options(repos = c(CRAN = "https://cloud.r-project.org"))

packages <- c(
  "reshape2",
  "ggplot2",
  "ggprism",
  "dplyr",
  "plyr",
  "igraph",
  "tidyverse",
  "ggraph",
  "magrittr"
)

installed <- rownames(installed.packages())
missing_pkgs <- packages[!(packages %in% installed)]

if(length(missing_pkgs) > 0){
  install.packages(missing_pkgs, dependencies = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c(
  "clusterProfiler",
  "DOSE",
  "enrichplot",
  "ComplexHeatmap",
  "pathview",
  "KEGGREST"
)

for(pkg in pkgs){
  if(!requireNamespace(pkg, quietly = TRUE)){
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

# if (!requireNamespace("clusterProfiler", quietly = TRUE))
#   BiocManager::install("clusterProfiler", update = FALSE, ask = FALSE)

#library(clusterProfiler)

# load related packages
suppressWarnings(suppressMessages(library("reshape2")))
suppressWarnings(suppressMessages(library("ggplot2")))
suppressWarnings(suppressMessages(library("ggprism")))
# suppressWarnings(suppressMessages(library("vegan")))
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("plyr")))
suppressWarnings(suppressMessages(library("igraph")))
suppressWarnings(suppressMessages(library("tidyverse")))
suppressWarnings(suppressMessages(library("ggraph")))
suppressWarnings(suppressMessages(library("magrittr")))
suppressWarnings(suppressMessages(library("clusterProfiler")))


#r.cor <- read.table("data3/r.cor_count_Y.txt", sep="\t", header=T, check.names=F,row.names = 1)
r.cor <- read.table(opts$Correlation, sep="\t", header=T, check.names=F,row.names = 1)
#p.cor <- read.table("data3/p.cor_count_Y.txt", sep="\t", header=T, check.names=F,row.names = 1)
p.cor <- read.table(opts$Pvalue, sep="\t", header=T, check.names=F,row.names = 1)

r.cor[p.cor>0.05] <- 0

# Build network connection attributes and node attributes
# Convert data to long format for merging and add connection properties
r.cor$from = rownames(r.cor)
p.cor$from = rownames(p.cor)
p_value <-  p.cor %>%
  gather(key = "to", value = "p", -from) %>%
  data.frame() 
p_value$FDR <- p.adjust(p_value$p,"BH")
p_value <- p_value[, -3]

cor.data<- r.cor %>%
  gather(key = "to", value = "r", -from) %>%
  data.frame() %>%
  left_join(p_value, by=c("from","to"))
cor.data <- as.data.frame(cor.data)
cor.data <- cor.data[cor.data$FDR <= 0.05 & cor.data$from != cor.data$to, ]
cor.data <- cor.data[abs(cor.data$r) >= 0.6 & cor.data$from != cor.data$to, ]
cor.data <- cor.data %>%
  plyr::mutate(
    linecolor = ifelse(r > 0,"positive","negative"),
    linesize = abs(r)
  )

# Set node properties
vertices <- c(as.character(cor.data$from),as.character(cor.data$to)) %>%
  as_tibble() %>%
  group_by(value) %>%
  clusterProfiler::summarise()
colnames(vertices) <- "name"

# Build graph data structure and add network basic attributes, save data
# Building a graph data structure
graph <- graph_from_data_frame(cor.data, vertices = vertices, directed = FALSE)
E(graph)$weight <- abs(E(graph)$r)
V(graph)$label <- V(graph)$name

# 利用“louvain”算法进行进行聚类群划分
# calculate community membership and modularity of networks
patients.clusters <- cluster_louvain(graph)
V(graph)$Cluster <- patients.clusters$membership

# save data
#write_graph(graph, file = paste(opts$output, "Centenarians_01.graphml", sep=""), format="graphml")
write_graph(
  graph,
  file = file.path(
    opts$output,
    paste0(opts$Group, "_01.graphml")
  ),
  format = "graphml"
)
