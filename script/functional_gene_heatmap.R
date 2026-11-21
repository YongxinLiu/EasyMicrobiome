#!/usr/bin/env Rscript

# Copyright 2016-2022 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, et al. 2021. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12: 315-330. https://doi.org/10.1007/s13238-020-00724-8

# 手动运行脚本，使用 Ctrl+Shift+H 或 Session 需要设置工作目录
# Set Work Directory - Choose Directory / To Source File Location

# Clean enviroment object
rm(list=ls()) 

# 1.1 简介 Introduction #----

# 程序功能：展示不同门水平功能基因的相对丰度
# Functions: Display the relative abundance of functional genes at different phylum levels



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
    make_option(c("-i", "--input"), type="character", default="result/dbcan3/ARG_F_new3.txt",
                help="MAGs taxonomy tables [default %default]"),
    make_option(c("-g", "--group"), type="character", default="result/dbcan3/group_ARG_new2.txt",
                help="MAGs taxonomy tables [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/dbcan3/",
                help="Display the relative abundance of functional genes at different phylum levels [default %default]") 
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
print("You are using the following parameters:")
print(opts)

# 设置输出文件缺省值，如果为空，则为输入+pcoa.pdf
if(opts$output == ""){
  opts$output = paste0(opts$input, ".heatmap.pdf")}


# Install related packages
if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggplot2","scales","dplyr","tidyr","ComplexHeatmap"))
}
# load related packages
suppressWarnings(suppressMessages(library("ggplot2")))
suppressWarnings(suppressMessages(library("scales")))
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("tidyr")))
suppressWarnings(suppressMessages(library("ComplexHeatmap")))


# 读取数据
#data3 <- read.table("data/ARG_F_new3.txt", header = TRUE, sep = "\t")
data3 <- read.table(opts$input, header = TRUE, sep = "\t")
data4 <- data3 %>% 
  select(-ORF_ID, -Antibiotic)

data5 <- data4 %>%
  group_by(AMR_Gene_Family) %>%
  summarise(across(everything(), sum, na.rm = TRUE))

data5 <- as.data.frame(data5)
rownames(data5) <- data5$AMR_Gene_Family
data5 <- data5[,-1]
data_p <- data5

# Decreased sort by abundance
mean_sort = data_p[(order(-rowSums(data_p))), ]
mean_sort = as.data.frame(mean_sort)
mean_sort2 = t(mean_sort)
mean_sort2 = mean_sort2[order(-mean_sort2[,1]),]
mean_sort3 = t(mean_sort2)
mean_sort3 = as.data.frame(mean_sort3)

# Filter Top 10
other = colSums(mean_sort3[11:dim(mean_sort3)[1], ])
mean_sort3 = mean_sort3[(11 - 1):1, ]
mean_sort3 = rbind(other,mean_sort3)
rownames(mean_sort3)[1] = c("others")
mean_sort3 = as.data.frame(mean_sort3)
mean_sort4 = mean_sort3[(order(-rowSums(mean_sort3))), ]
mean_sort4 <- mean_sort4[c(1:5, 7:11, 6), ]
data <- mean_sort4
data_row_select2 <- mean_sort4

# 绘图测试
pheatmap(data_row_select2, color = colorRampPalette(c("lightgrey", "firebrick3"))(100)) #出图

# 分组注释
#annotation_col<-read.table("data/group_ARG_new2.txt",header=T,sep="\t",row.names=1)
annotation_col<-read.table(opts$group,header=T,sep="\t",row.names=1)
annotation_col <- as.data.frame(annotation_col)

df1 <- annotation_col

# 计算 Phylum 的比例
phylum_percentage_df <- df1 %>%
  count(Group) %>%
  mutate(Percentage = n / sum(n) * 100) %>%
  arrange(desc(Percentage))

# 提取前15个 Phylum，并按字母顺序排列
top_phylum <- sort(phylum_percentage_df$Group[1:15])

# 用于替换其他的 Phylum
df2 <- df1 %>%
  mutate(Group = ifelse(Group %in% top_phylum, Group, "Others"))

annotation_col2 <- as.data.frame(df2)

data_row_select32 <- mean_sort4

group = annotation_col2
group$Group=as.factor(group$Group)

topanno=HeatmapAnnotation(df=group,#列注释
                          border = F,
                          show_annotation_name = F,
                          col = list(Group=c(
                            'Actinomycetota'="#fff7f3",
                            'Bacillota'="#c7eae5",
                            'Bacillota_A'="#deebf7",
                            'Bacillota_E'="#4292c6",
                            'Bacillota_G'="#2171b5",
                            'Bacillota_I'="#08519c",
                            'Bacteroidota'="#bc80bd",
                            'Bdellovibrionota'="#ffffb3",
                            'Chloroflexota'="#bebada",
                            'Gemmatimonadota'="#fb8072",
                            'Myxococcota'="#fdb462",
                            'Patescibacteria'="#b3de69",
                            'Planctomycetota'="#fccde5",
                            'Pseudomonadota'="#8dd3c7",
                            'Verrucomicrobiota'="#b35806",
                            'Others'= "#737373"
                          )),
                          gp =gpar(lwd = 0.0),
                          height = unit(0.01, "mm"),
                          annotation_height = unit(0.02, "mm"),
                          annotation_width = unit(0.02, "mm"),
                          simple_anno_size = unit(3.0, "mm"),
                          annotation_name_gp = gpar(fontsize = 10),
                          annotation_legend_param = list(title_gp = gpar(fontsize = 12),
                                                         labels_gp = gpar(fontsize = 12)
                          ),
                          labels = c(),
                          show_legend = T)

# 保存为pdf文件
pdf(file=paste(opts$output, "MAGs_functional_genes_heatmap.pdf", sep=""), width = 12, height = 3.5)
p2 = Heatmap(data_row_select32, row_names_side = "right", cluster_rows = FALSE, cluster_columns = TRUE, 
             row_names_gp = gpar(fontsize = 8),
             rect_gp = gpar(col = "white", lwd = 0.05),
             col = colorRampPalette(colors = c("lightgrey","#1E90FF"))(50),
             show_column_names = FALSE, top_annotation = topanno,
             column_split = group, 
             column_names_gp = gpar(fontsize = 0)
)
#p2
# 绘制热图
draw(p2)
#关闭设备并保存图像
dev.off()




