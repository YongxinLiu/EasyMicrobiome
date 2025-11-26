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
  make_option(c("-i", "--input"), type="character", default="result2/tax/otus.nwk",
              help="Phylogenetic tree [default %default]"),
  make_option(c("-n", "--anno"), type = "character", default = "result2/tax/annotation2.txt", 
              help = "Annotations"),
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

# 基于CRAN安装R包，检测没有则安装
p_list = c("ggtreeExtra","ggtree","treeio","tidytree","ggstar","ggplot2",
           "ggnewscale","TDbook")
for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
  library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

# 加载R包 Load the package
suppressWarnings(suppressMessages(library("ggtreeExtra")))
suppressWarnings(suppressMessages(library("ggtree")))
suppressWarnings(suppressMessages(library("treeio")))
suppressWarnings(suppressMessages(library("tidytree")))
suppressWarnings(suppressMessages(library("ggstar")))
suppressWarnings(suppressMessages(library("ggplot2")))
suppressWarnings(suppressMessages(library("ggnewscale")))
suppressWarnings(suppressMessages(library("TDbook")))

# 树文件和注释文件是通过EasyAmplicon分析流程获取的，EasyAmplicon网址：https://github.com/YongxinLiu/EasyAmplicon

# Load tree data
# 导入树文件
tree <- read.tree(opts$input)

# Load annotation data
# 导入注释文件
#data_anno <- read.table("Figure2/data/annotation2.txt", header = TRUE, sep = "\t")
data_anno <- read.table(opts$anno, header = TRUE, sep = "\t")

data1 <- subset(data_anno, select = c(OTUID, Phylum))
data1 <- data1[order(data1$Phylum), ]
color_palette <-c("#4393c3","#fec44f","#fa9fb5","#a1d99b",
                  "#dd1c77","#bcbddc","#1c9099","#bf812d",
                  "#c7eae5", "#80cdc1", "#01665e", "#de77ae", 
                  "#d73027", "#e0e0e0")

unique_phyla <- unique(data1$Phylum)
num_unique_phyla <- length(unique_phyla)
if (num_unique_phyla > length(color_palette)) {
  color_mapping <- setNames(rep(color_palette, length.out = num_unique_phyla), unique_phyla)
} else {
  color_mapping <- setNames(color_palette[1:num_unique_phyla], unique_phyla)
}

data1$color <- color_mapping[data1$Phylum]
data1 <- data1[order(data1$OTUID), ]

data2 <- subset(data_anno, select = c(OTUID, feces1, feces2))
data2$Prevalence<- rowMeans(data2[, c("feces1","feces2")])
data2$Group <- "Feces"

data3 <- subset(data_anno, select = c(OTUID, saliva1, saliva2))
data3$Prevalence<- rowMeans(data3[, c("saliva1", "saliva2")])
data3$Group <- "Saliva"

data4 <- subset(data_anno, select = c(OTUID, plaque1, plaque2))
data4$Prevalence<- rowMeans(data4[, c("plaque1", "plaque2")])
data4$Group <- "Plaque"

#如果数据相差太大，可以做Min-Max Normalization
min_max_normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
data2$Prevalence <- min_max_normalize(data2$Prevalence)
data3$Prevalence <- min_max_normalize(data3$Prevalence)
data4$Prevalence <- min_max_normalize(data4$Prevalence)
# Prepare the tree
data1$label <- data1$OTUID
tree2 <- full_join(tree, data1, by = "label")

# Create the circular layout tree
# 绘制系统发育树圈图
# mycol = c("#e6550d","#fd8d3c","#fdae6b","#fdd0a2","#31a354",
#           "#74c476","#a1d99b","#c7e9c0","#756bb1","#9e9ac8",
#           "#bcbddc","#dadaeb","#636363","#969696","#bdbdbd",
#           "#d9d9d9","#3182bd","#6baed6","#9ecae1","#c6dbef")
mycol = c("#bf812d", "#4393c3","#01665e","#de77ae","#d73027", "#80cdc1","#c7eae5", "#e0e0e0")
p1<- ggtree(tree2,
            aes(color=Phylum),#支长颜色按照分组进行着色
            #layout="fan",#进化树类型
            layout="circular",#进化树类型
            open.angle=0,#开口角度
            linewidth=0.45,#分支线条粗细
            branch.length = 'none',
            show.legend = F)+
  scale_color_manual(values = mycol, na.value = "#000000",
                     #guide="none"
                     guide=guide_legend(keywidth=0.5, keyheight=0.5, order=2),
                     name="Phylum")+
  new_scale_fill() +
  geom_fruit(
    data=data2,#数据
    geom = geom_tile,#绘图类型
    #mapping = aes(y=OTUID, x= Prevalence, fill = Group),
    mapping = aes(y=OTUID, fill= Prevalence),
    offset = 0.3,
    pwidth = 0.5,
    width=1.5,
    color = "black"
  )+
  scale_fill_gradient2(low = "#e0e6e6", high = "#d73027",midpoint = 0.005,
                       name = "\u0394 Feces")+
  geom_tiplab(size=2.4,
              fontface=2, offset=0.5)+
  theme(legend.text = element_text(colour="black", size = 8),
        legend.title = element_text(color = "black", size = 8, face = "bold"))+
  new_scale_fill() +
  geom_fruit(
    data=data3,#数据
    geom = geom_tile,#绘图类型
    #mapping = aes(y=OTUID, x= Prevalence, fill = Group),
    mapping = aes(y=OTUID, fill= Prevalence),
    offset = 0.09,
    pwidth = 0.5,
    width=1.5,
    color = "black"
  )+
  scale_fill_gradient2(low = "#e0e6e6",  high = "#01665e",midpoint = 0.005,
                       name = "\u0394 Saliva")+
  theme(legend.text = element_text(colour="black", size = 8),
        legend.title = element_text(color = "black", size = 8, face = "bold"))+
  new_scale_fill() +
  geom_fruit(
    data=data4,#数据
    geom = geom_tile,#绘图类型
    #mapping = aes(y=OTUID, x= Prevalence, fill = Group),
    mapping = aes(y=OTUID, fill= Prevalence),
    offset = 0.09,
    pwidth = 0.5,
    width=1.5,
    color = "black"
  )+
  scale_fill_gradient2(low = "#e0e6e6",  high = "#4393c3",midpoint = 0.000000005,
                       name = "\u0394 Plaque")+
  theme(legend.text = element_text(colour="black", size = 8),
        legend.title = element_text(color = "black", size = 8, face = "bold"))

pdf(paste0(opts$output, "phylogenetic_tree_without_distance.pdf"), width = 10, height = 8)  # Adjust width and height as needed
p1
dev.off()
