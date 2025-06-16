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
  make_option(c("-n", "--anno"), type = "character", default = "result2/tax/annotation3.txt", 
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

# 树文件和注释文件是通过EasyMetagenome分析流程获取的，EasyMetagenome网址：https://github.com/YongxinLiu/EasyMetagenome

# Load tree data
# 导入树文件
tree <- read.tree(opts$input)

# Load annotation data
# 导入注释文件
data_anno <- read.table(opts$anno, header = TRUE, sep = "\t")

data1 <- subset(data_anno, select = c(OTUID, Phylum, Class))
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
data1$color2 <- color_mapping[data1$Class]
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

data5 <- subset(data_anno, select = c(OTUID, All1, All2))
data5$Prevalence<- rowMeans(data5[, c("All1", "All2")])
data5$Group <- "All"

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
p1<- ggtree(tree2,
            aes(color=Phylum),#支长颜色按照分组进行着色
            layout="fan",#进化树类型
            left = TRUE,
            #layout="slanted",#进化树类型
            #layout="circular",#进化树类型
            open.angle=30,#开口角度
            linewidth=0.25,#分支线条粗细
            #branch.length = 'none',
            show.legend = F)+
  geom_aline(aes(color=Phylum),linetype = "solid",size=0.05)+
  scale_color_manual(values=c("#bf812d", "#c7eae5", "#a6bddb", "#01665e", "#de77ae", "#d73027", "#4393c3","#e0e0e0"))

p2 <- rotate_tree(p1, angle = 170)+
  new_scale_fill() +
  geom_fruit(
    data=data1,
    geom=geom_tile,
    mapping=aes(y=OTUID, fill=data1$Class),
    width=0.03,
    offset=0.05
  )+
  scale_fill_manual(
    values=c("#d2da93","#5196d5","#00ceff","#ff630d","#35978b",
             "#e5acd7","#77aecd","#ec8181","#dfc6a5","#e50719",
             "#d27e43","#8a4984","#fe5094","#8d342e","#f94e54",
             "#CAB2D6", "#FFFF99", "#8DD3C7", "#FFED6F", "#BC80BD",
             "#ffad00","#36999d","#00fc8d","#b64aa0","#9b82e1",
             "#66c2a4","#8c96c6","#fdbb84","#a6bddb","#fa9fb5"),
    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=2),
    name="Class")+
  #geom_tiplab(size=5,
  #            fontface=2, offset=4.5)+
  theme(legend.text = element_text(colour="black", size = 20),
        legend.title = element_text(color = "black", size = 20, face = "bold"))+
  new_scale_fill() +
  geom_fruit(
    data=data2,#数据
    geom = geom_col,#绘图类型
    mapping = aes(y=OTUID, x= Prevalence, fill = Group),
    offset = 0.04,
    pwidth = 0.1,
    width=0.5,
  )+
  scale_fill_manual(
    values=c("#4393c3","#fec44f","#fa9fb5","#a1d99b",
             "#dd1c77","#bcbddc"),
    guide=guide_legend(keywidth=1, keyheight=1, order=2),
    name="Sites1")+
  new_scale_fill() +
  geom_fruit(
    data=data3,#数据
    geom = geom_col,#绘图类型
    mapping = aes(y=OTUID, x= Prevalence, fill = Group),
    offset = 0.04,
    pwidth = 0.1,
    width=0.5,
  )+
  scale_fill_manual(
    values=c("#fec44f","#4393c3","#fa9fb5","#a1d99b",
             "#dd1c77","#bcbddc"),
    guide=guide_legend(keywidth=1, keyheight=1, order=2),
    name="Sites2")+
  new_scale_fill() +
  geom_fruit(
    data=data4,#数据
    geom = geom_col,#绘图类型
    mapping = aes(y=OTUID, x= Prevalence, fill = Group),
    offset = 0.04,
    pwidth = 0.1,
    width=0.5,
  )+
  scale_fill_manual(
    values=c("#fa9fb5","#fec44f","#4393c3","#a1d99b",
             "#dd1c77","#bcbddc"),
    guide=guide_legend(keywidth=1, keyheight=1, order=2),
    name="Sites3")+
  new_scale_fill() +
  geom_fruit(
    data=data5,#数据
    geom = geom_col,#绘图类型
    mapping = aes(y=OTUID, x= Prevalence, fill = Group),
    offset = 0.04,
    pwidth = 0.1,
    width=0.5,
  )+
  scale_fill_manual(
    values=c("#bcbddc","#fa9fb5","#fec44f","#4393c3","#a1d99b",
             "#dd1c77"),
    guide=guide_legend(keywidth=1, keyheight=1, order=2),
    name="Sites4")
#pdf(file=paste("Figure3/results/phylogenetic_tree01_04.pdf", sep=""), height = 10.2, width = 14)
pdf(paste0(opts$output, "phylogenetic_tree_with_distance.pdf"), width = 10.2, height = 14)  # Adjust width and height as needed
p2
dev.off()

