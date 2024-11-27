#!/usr/bin/env Rscript

# Copyright 2024 De-feng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, Yang Bai. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein Cell 2021(12) 5:315-330 doi: 10.1007/s13238-020-00724-8

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 更新
# 2021/5/31: 更新引文，添加数据表转置和标准化的选项

# 1.1 程序功能描述和主要步骤

# 程序功能：Alpha多样性箱线图+统计
# Functions: Alpha boxplot

options(warn = -1) # Turn off warning


## 设置输入输出文件和参数

# 修改下面`default=`后面的文件和参数。
#
# 输入文件为原始alpha多样性(vegan.txt)+分组信息(metadata.txt)
#
# 输入文件"-i", "--input"，result/alpha/vegan.txt; alpha多样性表格
#
# 实验设计"-d", "--design"，默认`metadata.txt`，可手动修改文件位置；
#
# 分组列名"-n", "--group"，默认将metadata.txt中的Group列作为分组信息，可修改为任意列名；
#
# 分组列名"-o", "--output"，默认为输出目录，图片文件名为alpha_boxplot_+多样性指数名+.pdf；统计文本位于代码运行目录中alpha_boxplot_TukeyHSD.txt；
#
# 图片宽"-w", "--width"，默认89 mm，根据图像布局可适当增大或缩小
#
# 图片高"-e", "--height"，默认59 mm，根据图像布局可适当增大或缩小


# 1.2 解析命令行
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T)
}


site = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"

options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = site)

a = rownames(installed.packages())

install_bioc <- c("ggplot2", "multcompView")

for (i in install_bioc) {
  if (!i %in% a)
    BiocManager::install(i, update = F, site_repository=site)
  a = rownames(installed.packages())
}

if (!"amplicon" %in% a){
  devtools::install_github("microbiota/amplicon")
}


# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="result12/metaphlan4/alpha.txt",
                help="Alpha diversity matrix [default %default]"),
    make_option(c("-a", "--alpha_index"), type="character", default="shannon",
                help="Group name [default %default]"),
    make_option(c("-d", "--design"), type="character", default="result12/metadata.txt",
                help="Design file or metadata [default %default]"),
    make_option(c("-t", "--transpose"), type="logical", default=FALSE,
                help="Design file or metadata [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Group",
                help="Group name [default %default]"),
    make_option(c("-s", "--scale"), type="logical", default=FALSE,
                help="Normalize to 100 [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result12/metaphlan4/",
                help="Output pdf directory, with prefix alpha_boxplot_; Stat in alpha_boxplot_TukeyHSD.txt [default %default]"),
    make_option(c("-x", "--xlabAngle"), type="logical", default=FALSE,
                help="X lab set in angle [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=89,
                help="Figure width in mm [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=59,
                help="Figure heigth in mm [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
suppressWarnings(dir.create(dirname(opts$output), showWarnings = F))


index = opts$alpha_index
# Source and edited from package amplicon
alpha_boxplot2 <- function(alpha_div, metadata, index = "shannon", groupID = "group", levels = c(), outlier = FALSE){
  #groupID = opts$group
  #levels = c("Cancer","Normal")
  #index = opts$alpha_index
  p_list = c("ggplot2", "dplyr", "multcompView")
  for (p in p_list) {
    if (!requireNamespace(p)) {
      install.packages(p)
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE, 
                                           quietly = TRUE, warn.conflicts = FALSE))
  }
  idx = rownames(metadata) %in% rownames(alpha_div)
  metadata = metadata[idx, , drop = F]
  alpha_div = alpha_div
  idx = rownames(metadata) %in% rownames(alpha_div)
  metadata = metadata[idx, , drop = F]
  alpha_div = alpha_div[rownames(metadata), ]
  sampFile = as.data.frame(metadata[, groupID], row.names = row.names(metadata))
  df = cbind(alpha_div[rownames(sampFile), index], sampFile)
  colnames(df) = c(index, "group")
  max = max(df[, c(index)])
  min = min(df[, index])
  x = df[, c("group", index)]
  y = x %>% group_by(group) %>% summarise_(Max = paste("max(", index, ")", sep = ""))
  y = as.data.frame(y)
  rownames(y) = y$group
  df$y = y[as.character(df$group), ]$Max + (max - min) * 0.05
  levels(as.factor(df$group))
  df = df %>%
    mutate(group = ordered(df$group,levels=levels))
  df$class = index
  compaired = list(c(levels[1], levels[2]))
  wt = wilcox.test(df[[index]] ~ df$group, alternative=c("two.sided"))
  FDR = p.adjust(wt$p.value, method = "BH")
  p1 = ggplot(df, aes(x = group, y = .data[[index]])) +
    geom_jitter(aes(color=group),position = position_jitter(0.15), size = 1.3, alpha = 1) +
    geom_boxplot(position=position_dodge(width =0.4),width=0.5, size = 0.4,
                 fill = "transparent", 
                 outlier.shape = NA,
                 linetype = "dashed", color = "black") +
    stat_boxplot(aes(ymin=..lower..,ymax=..upper..,
                     fill=group
    ),
    color="black",
    fill = "transparent",position=position_dodge(width =0.4),width=0.5, size = 0.4,outlier.shape = NA)+
    stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),
                 width=0.18,color="black",size = 0.4)+
    stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),
                 width=0.18,color="black",size = 0.4)+
    labs(x = NULL, y = NULL, color = groupID) + 
    scale_y_continuous(labels = label_number(accuracy = 0.1)) +
    scale_fill_manual(values = c("#00a9ba","#fe9390"))+
    scale_color_manual(values = c("#00a9ba","#fe9390"))+
    geom_signif(comparisons = compaired,
                step_increase = 0.3,
                map_signif_level = F,
                test = wilcox.test,
                color = "black",
                size = 0.2,
                textsize = 3
    )+
    mytheme+
    facet_grid(.~class)
  p1
}


# 2. 依赖关系检查、安装和加载

suppressWarnings(suppressMessages(library(amplicon)))
suppressWarnings(suppressMessages(library(scales)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(ggtern)))
suppressWarnings(suppressMessages(library(ggplot2)))

mytheme = theme_bw() + theme(text = element_text(family = "sans", size = 8))+
  theme(legend.position="none",
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=10, colour="black", family = "sans", angle = 0),
        axis.text.x = element_text(size=10, colour="black", family = "sans", angle = 0, hjust = 0),
        axis.title= element_text(size=10, family = "sans"),
        strip.text.x = element_text(size=10, angle = 0),
        strip.text.y = element_text(size=10, angle = 0),
        plot.title = element_text(size=10, angle = 0),
        strip.background.x = element_rect(fill = "#E5E4E2", colour = "black", size = 0.2))+
  theme(axis.text.x=element_text(angle=0,vjust=1, hjust=0.6))+
  theme(axis.line = element_line(size = 0.1, colour = "black"))

# 3. 读取输入文件

# 读取OTU表
alpha_div = read.table(opts$input, header=T, row.names=1, sep="\t", comment.char="")


# 条件判断是否转置
if (opts$transpose){
  alpha_div = as.data.frame(t(alpha_div))
}

# 条件判断是否标准化
if (opts$scale){
  alpha_div = alpha_div/rowSums(alpha_div,na=T)*100
}

# 读取实验设计
metadata = read.table(opts$design, header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)

p = alpha_boxplot2(alpha_div, metadata, index = opts$alpha_index, groupID = opts$group, levels = c("Cancer","Normal"))
if (opts$xlabAngle){
  p = p + theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
}
# Saving figure
# 保存图片，大家可以修改图片名称和位置，长宽单位为毫米
ggsave(paste0(opts$output,"p_boxplot_",opts$alpha_index,".pdf"), p, width = opts$width, height = opts$height, units = "mm")


