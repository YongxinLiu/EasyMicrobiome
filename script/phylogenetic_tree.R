#!/usr/bin/env Rscript

# Copyright 2024 De-feng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, et al. 2021. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12: 315-330. https://doi.org/10.1007/s13238-020-00724-8

# 手动运行脚本，使用 Ctrl+Shift+H 或 Session 需要设置工作目录
# Set Work Directory - Choose Directory / To Source File Location

# Clean enviroment object
rm(list=ls()) 

# 1.1 简介 Introduction #----

# 程序功能：微生物物种系统发育树
# Functions: Phylogenetic analysis for microbiota species


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
    make_option(c("-p", "--input"), type="character", default="result/gtdb_95/tax.unrooted.tree",
                help="MAGs unrooted tree [default %default]"),
    make_option(c("-a", "--annotation"), type="character", default="result/gtdb_95/annotation2.txt",
                help="MAGs annotations [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/gtdb_95/",
                help="Output Phylogenetic tree for microbiota species [default %default]") 
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
print("You are using the following parameters:")
print(opts)


# Install related packages
if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggtreeExtra","ggtree","treeio","tidytree","ggstar","ggplot2",
             "ggnewscale","TDbook"))
}
# load related packages
suppressWarnings(suppressMessages(library("ggtreeExtra")))
suppressWarnings(suppressMessages(library("ggtree")))
suppressWarnings(suppressMessages(library("treeio")))
suppressWarnings(suppressMessages(library("tidytree")))
suppressWarnings(suppressMessages(library("ggstar")))
suppressWarnings(suppressMessages(library("ggplot2")))
suppressWarnings(suppressMessages(library("ggnewscale")))
suppressWarnings(suppressMessages(library("TDbook")))



# Load tree data
#tree <- read.tree("tax.unrooted.tree")
tree <- read.tree(opts$input)

#data_anno <- read.table("annotation2.txt", header = TRUE, sep = "\t")
data_anno <- read.table(opts$annotation, header = TRUE, sep = "\t")

data1 <- subset(data_anno, select = c(ID, Phylum))
data1 <- data1[order(data1$Phylum), ]
color_palette <-c("#bf812d", "#c7eae5", "#80cdc1", "#01665e", "#de77ae", "#d73027", "#4393c3","#e0e0e0")
unique_phyla <- unique(data1$Phylum)
num_unique_phyla <- length(unique_phyla)
if (num_unique_phyla > length(color_palette)) {
  color_mapping <- setNames(rep(color_palette, length.out = num_unique_phyla), unique_phyla)
} else {
  color_mapping <- setNames(color_palette[1:num_unique_phyla], unique_phyla)
}
data1$color <- color_mapping[data1$Phylum]
data1 <- data1[order(data1$ID), ]
data2 <- subset(data_anno, select = c(ID, MY1, MY2, MY3, FY1, FY2, FY3))
data2$Prevalence<- rowMeans(data2[, c("MY1", "MY2", "MY3","FY1", "FY2", "FY3")])
data2$Group <- "Young"
data3 <- subset(data_anno, select = c(ID, ME1, ME2, ME3, FE1, FE2, FE3 ))
data3$Prevalence<- rowMeans(data3[, c("ME1", "ME2", "ME3", "FE1","FE2","FE3")])
data3$Group <- "Elderly"
data4 <- subset(data_anno, select = c(ID, MC1, MC2, MC3, FC1, FC2, FC3))
data4$Prevalence<- rowMeans(data4[, c("MC1", "MC2", "MC3", "FC1","FC2","FC3")])
data4$Group <- "Centenarian"
#如果数据相差太大，可以做Min-Max Normalization

min_max_normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
data2$Prevalence <- min_max_normalize(data2$Prevalence)
data3$Prevalence <- min_max_normalize(data3$Prevalence)
data4$Prevalence <- min_max_normalize(data4$Prevalence)
# Prepare the tree
data1$label <- data1$ID
tree2 <- full_join(tree, data1, by = "label")
# Create the circular layout tree
#p <- ggtree(tree2, layout = "fan", size = 0.15, open.angle = 0)
#p

p2<- ggtree(tree2,
            aes(color=Phylum),#支长颜色按照分组进行着色
            layout="fan",#进化树类型
            open.angle=0,#开口角度
            linewidth=0.75,#分支线条粗细
            show.legend = F)+
  #geom_tiplab(aes(color = label %in% df_label),#设定标签颜色根据筛选条件突出显示特定标签
  #            size=3.5,#字体大小
  #            align = T,#使用虚线连接标签与分支
  #            linetype = 3,linewidth = 0.4,offset = 12.5,show.legend = F)
  new_scale_fill() +
  geom_fruit(
    data=data2,#数据
    geom = geom_col,#绘图类型
    mapping = aes(y=ID, x= Prevalence, fill = Group),
    offset = 0.1,
    pwidth = 0.1,
    width=0.5,
    #show.legend = FALSE
  )+
  scale_fill_manual(
    #values=c("#4285f4", "#34a853", "#fbbc05","#ea4335"),
    values=c("#6a67ce","#ffb900","#fc636b","#aeb6b8","#e53238","lightblue"),
    #values=c("red","black"),
    guide=guide_legend(keywidth=1, keyheight=1, order=2),
    name="Cooling")+
  new_scale_fill() +
  geom_fruit(
    data=data3,#数据
    geom = geom_col,#绘图类型
    mapping = aes(y=ID, x= Prevalence, fill = Group),
    offset = 0.01,
    pwidth = 0.1,
    width=0.5,
    #show.legend = FALSE
  )+
  scale_fill_manual(
    #values=c("#4285f4", "#34a853", "#fbbc05","#ea4335"),
    values=c("#ffb900","#fc636b","#aeb6b8","#e53238","lightblue"),
    #values=c("red","black"),
    guide=guide_legend(keywidth=1, keyheight=1, order=2),
    name="Mature")+
  new_scale_fill() +
  geom_fruit(
    data=data4,#数据
    geom = geom_col,#绘图类型
    mapping = aes(y=ID, x= Prevalence, fill = Group),
    offset = 0.01,
    pwidth = 0.1,
    width=0.5,
    #show.legend = FALSE
  )+
  #scale_color_manual(values=c("black","#1aafd0","#6a67ce","#ffb900","#fc636b","#aeb6b8","#e53238"))+
  scale_color_manual(values=c("#bf812d", "#c7eae5", "#80cdc1", "#01665e", "#de77ae", "#d73027", "#4393c3","#e0e0e0"))+
  new_scale_fill() +
  geom_fruit(
    data=data1,
    geom=geom_tile,
    mapping=aes(y=ID, fill=data1$Phylum),
    #color="grey10",
    width=0.05,
    offset=0.1#,
    #show.legend = FALSE
  )+
  scale_fill_manual(
    #values=c("#4285f4", "#34a853", "#fbbc05","#ea4335"),
    values=c("#bf812d", "#c7eae5", "#80cdc1", "#01665e", "#de77ae", "#d73027", "#4393c3","#e0e0e0"),
    #values=c("red","black"),
    guide=guide_legend(keywidth=1, keyheight=1, order=2),
    name="Taxonomy")
pdf(file=paste(opts$output, "MAGs_phylogenetic_tree01.pdf", sep=""), height = 5.2, width = 11)
p2
dev.off()

