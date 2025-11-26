#!/usr/bin/env Rscript

# Copyright 2016-2021 Yong-Xin Liu <yxliu@genetics.ac.cn / metagenome@126.com>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, Yang Bai. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein Cell 2021(12) 5:315-330 doi: 10.1007/s13238-020-00724-8

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

#----1. 参数 Parameters#----

#----1.1 功能描述 Function description#----

# 程序功能：物种组成堆叠柱状图
# Functions: Taxonomy stackplot

options(warn = -1) # Turn off warning


## 设置输入输出文件和参数

# 修改下面`default=`后面的文件和参数。
#
# 输入文件为物种相对丰度矩阵(sum_p/c/o/f/g.txt)+分组信息(metadata.txt)
#
# 输入文件"-i", "--input"，tax/sum_p.txt; 物种组成表
#
# 实验设计"-d", "--design"，默认`metadata.txt`，可手动修改文件位置；
#
# 分组列名"-n", "--group"，默认将metadata.txt中的Group列作为分组信息，可修改为任意列名；
#
# 图片宽"-w", "--width"，默认89 mm，根据图像布局可适当增大或缩小
#
# 图片高"-e", "--height"，默认59 mm，根据图像布局可适当增大或缩小


#----1.2 参数缺省值 Default values#----
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="metaphlan4/Phylum.txt",
                help="Taxonomy composition [default %default]"),
    make_option(c("-d", "--design"), type="character", default="metadata.txt",
                help="Design file or metadata [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Group",
                help="Group name [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="Output directory; name according to input [default %default]"),
    make_option(c("-l", "--legend"), type="numeric", default=6,
                help="Legend number [default %default]"),
    make_option(c("-c", "--color"), type="character", default="Paired",
                help="color ggplot, manual1, Paired or Set3 [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=181,
                help="Figure width in mm [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=118,
                help="Figure heidth in mm [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  # suppressWarnings(dir.create(opts$output))
}
# 设置输出文件缺省值，如果为空，则为输入+pcoa.pdf
if(opts$output==""){opts$output=opts$input}


#----1.3. 加载包 Load packages#----

suppressWarnings(suppressMessages(library(amplicon)))
suppressWarnings(suppressMessages(library(patchwork)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(magrittr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggprism)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(plyr)))

#----2. 读取文件 Read files#----

#----2.1 实验设计 Metadata#----
metadata = read.table(opts$design, header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)

#----2.2 物种组成矩阵Taxonomy matrix#----
#taxonomy = read.table(opts$input, header=T, row.names=1, sep="\t", comment.char="", quote = "")
taxonomy = read.table(opts$input, header=T, sep="\t", comment.char="", quote = "")


# sum
# 计算微生物相对丰度之和，避免有重复统计
data = taxonomy
data <- aggregate(.~ Taxonomy,data=data,sum)
rownames(data) = data$Taxonomy
data = data[, -1]

# 计算相对丰度
# Calculate relative abundance
data = apply(data , 2, function(x) x/sum(x))
data = as.data.frame(data)

# Decreased sort by abundance
# 相对丰度按降序排列
mean_sort = data[(order(-rowSums(data))), ]
mean_sort = as.data.frame(mean_sort)
mean_sort2 = t(mean_sort)
mean_sort2 = mean_sort2[order(-mean_sort2[,1]),]
mean_sort3 = t(mean_sort2)
mean_sort3 = as.data.frame(mean_sort3)

# Phylum水平展示前5个
# Top 5
other = colSums(mean_sort3[opts$legend:dim(mean_sort3)[1], ])
mean_sort3 = mean_sort3[(opts$legend - 1):1, ]
mean_sort3 = rbind(other,mean_sort3)
rownames(mean_sort3)[1] = c("others")
mean_sort3 = as.data.frame(mean_sort3)

# Add taxonomy
# 加入微生物分类信息
mean_sort3$tax = rownames(mean_sort3)
data_all = as.data.frame(reshape2::melt(mean_sort3, id.vars = c("tax")))
data_all$group = data_all$variable
data_all$group = as.character(data_all$group)
data_all$group = gsub("[0-9]","", data_all$group)

# 从 metadata 建立映射表：首字母（大写） -> 完整组名（取第一个匹配）
map <- metadata %>%
  mutate(prefix = toupper(substr(Group, 1, 1))) %>%
  # 若同一首字母有多个完整组，默认保留第一个；如需其它策略可改这里
  group_by(prefix) %>%
  slice(1) %>%
  ungroup() %>%
  select(prefix, full_group = Group)

# 将 data_all 的 group 首字母转换并左连接映射表，然后用完整名替换（若找不到则保留原值）
data_all <- data_all %>%
  mutate(prefix = toupper(substr(group, 1, 1))) %>%
  left_join(map, by = "prefix") %>%
  mutate(group = if_else(!is.na(full_group), full_group, group)) %>%
  select(-prefix, -full_group)

# 给分组排序
# Sort for different groups
levels(as.factor(data_all$group))

data_all2 <- data_all %>%
  mutate(group = factor(group, levels = unique(group)))


# 根据样本数量确定每个分面的宽度，图例在顶部
# Determine the width of each facet based on the number of samples, the legend is at the top
plots <- lapply(split(data_all2, data_all2$group), function(df) {
  group_name <- unique(df$group)
  ggplot(df, aes(x = factor(variable, levels = unique(df$variable)),
                 y = value, fill = factor(tax, levels = unique(df$tax)))) +
    geom_bar(stat = "identity", position = "stack", width = 1) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    theme_classic() +
    labs(x = group_name, y = NULL) +
    scale_fill_manual(values = c("#e5acd7",  "#00ceff", "#ff630d", "#35978b","#d2da93",
                                 "#5196d5", "#77aecd", "#ec8181", "#dfc6a5", "#e50719",
                                 "#d27e43", "#8a4984", "#fe5094", "#8d342e", "#f94e54",
                                 "#ffad00", "#36999d", "#00fc8d", "#b64aa0", "#9b82e1")) +
    guides(fill = guide_legend(title = "Phylum"))
})

# 移除后三个分面的所有 y 轴元素和图例
# Remove all y-axis elements and legends for the last three facets
for (i in 2:2) {
  plots[[i]] <- plots[[i]] + theme(axis.text.y = element_blank(),
                                   axis.text.x = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.ticks.x = element_blank(),
                                   axis.title.y = element_blank(),
                                   axis.line.y = element_blank(),
                                   legend.position = "none")
}

# 为第一个分面保留 y 轴标签和图例
# Keep y-axis label and legend for the first facet
plots[[1]] <- plots[[1]] + 
  ylab("Percentage (%)") + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "left",
        legend.justification = c("left", "top"))

# 每个分面的宽度由样本数量决定
# The width of each facet is determined by the number of samples
sample_counts <- table(data_all2$group)
relative_widths <- sample_counts / sum(sample_counts)

# 使用 patchwork 组合图形，设置每个分面的宽度并统一图例
# Use patchwork to combine graphics, set the width of each facet and unify the legend
p <- wrap_plots(plots) +
  plot_layout(widths = relative_widths, guides = "collect") &
  theme(legend.position = "top", 
        legend.justification = "center",
        legend.direction = "horizontal", 
        legend.key.size = unit(0.3, "cm"), 
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.1, "cm"),
        axis.title.y = element_text(size = 10),
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) 

# 保存图像
# Save plot
#ggsave("results/Phylum_fungi_top5_3.pdf", p03, width = 139 * 1.5, height = 80 * 1.5, unit = 'mm')
#---3.2 保存 Saving#----
# 大家可以修改图片名称和位置，长宽单位为毫米
ggsave(paste0(opts$output,".sample.order.pdf"), p, width = opts$width, height = opts$height, units = "mm")
