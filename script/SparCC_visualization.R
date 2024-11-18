#!/usr/bin/env Rscript

# Copyright 2024 De-feng Bai <baidefeng@caas.cn>

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
    make_option(c("-i", "--Correlation"), type="character", default="result12/metaphlan4/sxtr_cov_mat_Cancer.tsv",
                help="Metaphlan4 species table"),
    make_option(c("-P", "--Pvalue"), type="character", default="result12/metaphlan4/sxtr_pvals_cancer.two_sided.tsv",
                help="Unfiltered OTU table [default %default]"),
    make_option(c("-r", "--output"), type="character", default="result12/metaphlan4/",
                help="Output file for SparCC analysis in different groups  [default %default]") 
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
print("You are using the following parameters:")
print(opts)


# Install related packages
if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("reshape2","ggplot2","ggprism","dplyr","plyr",
             "igraph","tidyverse","ggraph")) # ,"vegan"
}
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
write_graph(graph, file = paste(opts$output, "Cancer_01.graphml", sep=""), format="graphml")

# 可视化方式1：基于Gephi软件进行可视化 https://gephi.org/
# Visualized in Gephi software
# The same procedure for healthy group

# healthy组和patients组相同
# healthy.clusters <- cluster_louvain(healthy.igraph.s)
# V(healthy.igraph.s)$Cluster <- healthy.clusters$membership

# 可视化方式2：利用igraph进行可视化
g <- graph
# 准备网络图布局数据
# Preparing network diagram layout data。
layout1 <- layout_in_circle(g) 
layout5 <- layout_with_graphopt(g)

# 设置绘图颜色
# Setting the drawing color
#color <- c("#879b56","#ce77ad","#5ea6c2")

color = c("#d2da93","#5196d5","#00ceff","#ff630d","#35978b",
          "#e5acd7","#77aecd","#ec8181","#dfc6a5","#e50719",
          "#d27e43","#8a4984","#fe5094","#8d342e","#f94e54",
          "#ffad00","#36999d","#00fc8d","#b64aa0","#9b82e1")

names(color) <- unique(V(g)$Cluster) 
V(g)$point.col <- color[match(V(g)$Cluster,names(color))] 

# 边颜色按照相关性正负设置
# The edge color is set according to the positive or negative correlation
E(g)$color <- ifelse(E(g)$linecolor == "positive","#ff878c","#5ea6c2")

pdf(file=paste(opts$output, "network_group_Cancer_cluster6.pdf", sep=""), width=10, height=12)
par(mar=c(5,2,1,2))
plot.igraph(g, layout=layout5,
            vertex.color=V(g)$point.col,
            vertex.border=V(g)$point.col,
            vertex.size=6,
            vertex.frame.color="white",
            vertex.label=g$name,
            vertex.label.cex=0.8,
            vertex.label.dist=0, 
            vertex.label.degree = pi/2,
            vertex.label.col="black",
            edge.arrow.size=0.5,
            edge.width=abs(E(g)$r)*15,
            edge.curved = FALSE
)

# 设置图例
legend(
  title = "Cluster",
  list(x = min(layout1[,1])-0.05,
       y = min(layout1[,2])-0.05), 
  legend = c(unique(V(g)$Cluster)),
  fill = color,
  #pch=1
)

legend(
  title = "|r-value|",
  list(x = min(layout1[,1])+0.6,
       y = min(layout1[,2])-0.05),
  legend = c(0.2,0.4,0.6,0.8,1.0),
  col = "black",
  lty=1,
  lwd=c(0.2,0.4,0.6,0.8,1.0)*4,
)

legend(
  title = "Correlation (±)",
  list(x = min(layout1[,1])+1.0,
       y = min(layout1[,2])-0.05),
  legend = c("positive","negative"),
  col = c("#ff878c",rgb(0,147,0,maxColorValue = 255)),
  lty=1,
  lwd=1
)
dev.off()

