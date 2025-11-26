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

#各个物种占比注释

data3 <- data_anno[, c("ID", "Phylum")]

# 创建一个数据框存储门的名称和相应的百分比#后续可手动将比例尺加到图上
phylum_counts <- table(data3$Phylum)
phylum_proportion <- prop.table(phylum_counts)
phylum_percentage <- phylum_proportion * 100
phylum_percentage_df <- data.frame(
  Phylum = names(phylum_percentage), 
  Percentage = as.numeric(phylum_percentage)
)

phylum_percentage_sorted <- phylum_percentage_df[order(-phylum_percentage_df$Percentage), ]
# 提取占比多的Phylum
top_phylum <- phylum_percentage_sorted$Phylum[1:15]
# 创建新的列，将其他替换为 'others'
data4 <- data3
data4$Phylum <- ifelse(data4$Phylum %in% top_phylum, data4$Phylum, "Others")
# 检查结果
head(data4)

data1 <- subset(data4, select = c(ID, Phylum))
data1 <- data1[order(data1$Phylum), ]
color_palette <- c("#fff7f3", "#deebf7","#9ecae1", "#4292c6", "#2171b5", "#08519c","#bc80bd", "#ffffb3", 
                   "#bebada", "#fb8072", "#fdb462","#b3de69","#fccde5","#8dd3c7","#b35806","#737373")
unique_phyla <- unique(data1$Phylum)
num_unique_phyla <- length(unique_phyla)
if (num_unique_phyla > length(color_palette)) {
  color_mapping <- setNames(rep(color_palette, length.out = num_unique_phyla), unique_phyla)
} else {
  color_mapping <- setNames(color_palette[1:num_unique_phyla], unique_phyla)
}
data1$color <- color_mapping[data1$Phylum]
data1 <- data1[order(data1$ID), ]
data2 <- subset(data_anno, select = c(ID, CO1, CO2, CO3, CO4,MA1, MA2, MA3, MA4,ME1, ME2, ME3, ME4,TH1, TH2, TH3, TH4))
data2$Prevalence<- rowMeans(data2[, c("CO1", "CO2", "CO3", "CO4","MA1", "MA2", "MA3", "MA4","ME1", "ME2", "ME3", "ME4","TH1", "TH2", "TH3", "TH4")])
data2$Group <- "ALL"

# Prepare the tree
data1$label <- data1$ID
tree2 <- full_join(tree, data1, by = "label")
# Create the circular layout tree
p <- ggtree(tree2, layout = "fan", size = 0.15, open.angle = 0)

p
top_phylum_sorted <- sort(top_phylum)
data1$Phylum <- factor(data1$Phylum, levels = c(top_phylum_sorted, "Others"))
p2<- ggtree(tree2,
            aes(color=Phylum),#支长颜色按照分组进行着色
            layout="fan",#进化树类型
            open.angle=0,#开口角度
            linewidth=0.3,#分支线条粗细
            show.legend = F)+
  #geom_tiplab(aes(color = label %in% df_label),#设定标签颜色根据筛选条件突出显示特定标签
  #            size=3.5,#字体大小
  #            align = T,#使用虚线连接标签与分支
  #            linetype = 3,linewidth = 0.4,offset = 12.5,show.legend = F)+
  #scale_color_manual(values=c("black","#1aafd0","#6a67ce","#ffb900","#fc636b","#aeb6b8","#e53238"))+
  scale_color_manual(values=c("#fff7f3", "#deebf7","#9ecae1", "#4292c6", "#2171b5", "#08519c","#bc80bd", "#ffffb3",  
                              "#bebada", "#fb8072", "#fdb462","#b3de69","#fccde5","#8dd3c7","#b35806", "#737373"))+
  new_scale_fill() +
  geom_fruit(
    data=data2,#数据
    geom = geom_col,#绘图类型
    mapping = aes(y=ID, x= Prevalence, fill = Group),
    offset = 0.03,
    pwidth = 0.5,
    width=0.6,
    #show.legend = FALSE
  )+
  scale_fill_manual(
    #values=c("#4285f4", "#34a853", "#fbbc05","#ea4335"),
    values=c("black","#6a67ce","#ffb900","#fc636b","#aeb6b8","#e53238","lightblue"),
    #values=c("red","black"),
    guide=guide_legend(keywidth=1, keyheight=1, order=2),
    name="Cooling")+
  new_scale_fill() +
  geom_fruit(
    data=data1,
    geom=geom_tile,
    mapping=aes(y=ID, fill=data1$Phylum),
    #color="grey10",
    width=0.03,
    offset=-0.15#,
    #show.legend = FALSE
  )+
  scale_fill_manual(
    #values=c("#4285f4", "#34a853", "#fbbc05","#ea4335"),
    values=c("#fff7f3", "#deebf7","#9ecae1", "#4292c6", "#2171b5", "#08519c","#bc80bd", "#ffffb3", 
             "#bebada", "#fb8072", "#fdb462","#b3de69","#fccde5","#8dd3c7","#b35806","#737373"),
    #values=c("red","black"),
    guide=guide_legend(keywidth=1, keyheight=1, order=2),
    name="Host taxonomy")
#p2
pdf(file=paste(opts$output, "MAGs_phylogenetic_tree02.pdf", sep=""), height = 5.2, width = 11)
p2
dev.off()

