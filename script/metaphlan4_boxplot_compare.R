#!/usr/bin/env Rscript
# 
# Copyright 2024-2026 Defeng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Bai, et al. 2025. EasyMetagenome: A User‐Friendly and Flexible Pipeline for Shotgun Metagenomic Analysis in Microbiome Research. iMeta 4: e70001. https://doi.org/10.1002/imt2.70001

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 更新
# 2024/11/12：增加利用Metaphlan4相对丰度表组间比较箱线图加统计显著性
# 2025/11/25：规范代码


# 1.1 程序功能描述和主要步骤

# 程序功能：物种相对丰度组间比较箱线图
# Functions: Box plot comparing relative abundance of species between groups

# 主要步骤: 
# - 读取Metaphlan2物种组表 result/metaphlan2/taxonomy.spf(己标准化为100)
# - 物种组成按丰度均值排序
# - 筛选Top N个物种绘制热图

options(warn = -1) # Turn off warning

## 设置输入输出文件和参数

# 修改下面`default=`后面的文件和参数。
#
# 输入文件为原始物种相对丰度(metaphlan4/taxonomy.spf)+分组信息(metadata.txt)
#
# 输入文件"-i", "--input"，metaphlan4/taxonomy.tsv; 微生物相对丰度表格
#
# 分组列名"-t", "--taxonomy"，默认'Phylum'，代表门级别；
#
# 输入文件"-n", "--TopN"，默认'25'; 相对丰度top 25;
#
# 实验设计"-g", "--metadata"，默认`metadata.txt`，可手动修改文件位置；
#
# 图片宽"-w", "--width"，默认89 mm，根据图像布局可适当增大或缩小
#
# 图片高"-e", "--height"，默认59 mm，根据图像布局可适当增大或缩小
# 
# 分组列名"-o", "--output"，默认为输出目录；


# 1.2 依赖包安装

site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
a = rownames(installed.packages())

# install CRAN
install_CRAN <- c("ggplot2","dplyr","reshape2","ggpubr", "optparse", "tidyr")
for (i in install_CRAN) {
  if (!i %in% a)
    install.packages(i, repos = site)
  require(i,character.only=T)
  a = rownames(installed.packages())
}

# install bioconductor
install_bioc <- c("ggplot2", "multcompView")
for (i in install_bioc) {
  if (!i %in% a)
    BiocManager::install(i, update = F)
  a = rownames(installed.packages())
}

# install github
if (!"amplicon" %in% a){
  devtools::install_github("microbiota/amplicon")
}


# 1.3 解析命令行
# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="metaphlan4/taxonomy.spf", help="Metaphlan2 [default %default]"),
    make_option(c("-t", "--taxonomy"), type="character", default="Phylum", help="Taxonomy level [default %default]"),
    make_option(c("-n", "--TopN"), type="numeric", default="25", help="Number of taxonomy showing [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=8,
                help="Width of figure [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=5,
                help="Height of figure [default %default]"),
    make_option(c("-o", "--output"), type="character", default="", help="Output boxplot filename [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  
  # 调置如果无调设置输出，根据其它参数设置默认输出
  prefix = gsub("taxonomy.spf$", "", opts$input, perl = T)
  if (opts$output==""){opts$output=paste0(prefix, "boxplot", opts$taxonomy, opts$TopN)}
  
  # 显示输入输出确认是否正确
  # Metaphlan2物种组成表
  print(paste("The input file: ", opts$input,  sep = ""))
  # 绘制的分类级别, 默认为种
  print(paste("Taxonomy level: ", opts$taxonomy, ". Default if Species", sep = ""))
  # 选择绘制高丰度物种数量，默认30，最大为物种级别非冗余条目数量
  print(paste("Number of taxonomy showing: ", opts$TopN,  sep = ""))
  # 输出文件名，不填则为输入目录+boxplot+taxonomy
  print(paste("Output boxplot filename: ", opts$output, sep = ""))
}


# 2. 依赖关系检查、安装和加载

# 依赖包列表
package_list <- c(
  "ggplot2","dplyr","reshape2","ggpubr", "optparse"
)

# 批量安装和加载
for (p in package_list) {
  # 如果未安装，则安装
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  }
  
  # 批量加载，抑制警告和消息
  suppressWarnings(
    suppressMessages(
      library(p, character.only = TRUE)
    )
  )
}


# 设置绘图主题
mytheme = theme_bw() + theme(text = element_text(family = "sans", size = 6))+
  theme(#legend.position="none",
    legend.text = element_text(size=12),
    legend.title = element_blank(), 
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(size=12, colour="black", family = "sans", angle = 0), 
    axis.text.x = element_text(size=12, colour="black", family = "sans", angle = 0, hjust = 0),
    axis.title= element_text(size=12),
    strip.text.x = element_text(size=12, angle = 0),
    strip.text.y = element_text(size=12, angle = 0),
    plot.title = element_text(size=12, angle = 0),
    strip.background.x = element_rect(fill = "#E5E4E2", colour = "black", size = 0.2))+
  theme(axis.text.x=element_text(angle=0,vjust=1, hjust=0.6))+
  theme(axis.line = element_line(size = 0.1, colour = "black"))


# 3. 读取输入文件

# 读取metaphlan2文件
# 默认的quote会跳过2/3的数据，导致行减少产生NA，改默认值为空
taxonomy = read.table(opts$input, header=T, sep="\t", quote = "", row.names=NULL, comment.char="")
print(paste0("All taxonomy annotations are ", dim(taxonomy)[1], " lines!"))

# 去除NA，否则无法计算
taxonomy = na.omit(taxonomy)

# 显示样本总数据，有冗余
# colSums(taxonomy)


# 4. 差异比较并绘图

# 按指定组合并
grp = taxonomy[, opts$taxonomy, drop=F]
abu = taxonomy[,9:dim(taxonomy)[2]]
merge = cbind(abu, grp)

# group_by传变量，前面加".dots="
mergeTax = merge %>% group_by(.dots=opts$taxonomy) %>% summarise_all(sum)

# 合并后表格转换为数据框
mergeTax = as.data.frame(mergeTax)

# 按丰度排序
idx = order(rowMeans(mergeTax[,2:dim(mergeTax)[2]]), decreasing = T)
mergeTax = mergeTax[idx,]

# 添加行名
rownames(mergeTax)=mergeTax[,1]

# 筛选TopN绘图
# remove rownames line
Top = mergeTax[,-1]
# normalization to percentage 100
#Top = as.data.frame(t(t(Top)/colSums(Top,na=T) * 100))
Top = as.data.frame(t(t(Top)/colSums(Top,na=T)))

# Select TopN line for plotting
Top = head(Top, n=opts$TopN)
Top = cbind(row.names(Top), Top)
colnames(Top)[1] = opts$taxonomy

# melt有问题
# data_all = as.data.frame(melt(Top, id.vars = opts$taxonomy))
# 
# # set taxonomy order by abundance, default by alphabet
# data_all[[opts$taxonomy]]  = factor(as.character(data_all[[opts$taxonomy]]), levels=rownames(Top))
# 
# data_all$group = data_all$variable
# data_all$group = gsub("[0-9]","", data_all$group)
# data_all = data_all[, -2]


# pivot_longer 替代 melt
data_all <- Top %>%
  pivot_longer(
    cols = -!!sym(opts$taxonomy),        # 除了 taxonomy 列之外的所有列
    names_to = "variable",               # 原来的 melt variable.name
    values_to = "value"                  # 原来的 melt value.name
  )

# set taxonomy order by abundance, default by rownames
data_all[[opts$taxonomy]] <- factor(as.character(data_all[[opts$taxonomy]]), levels = Top[[opts$taxonomy]])

# 提取 group 信息
data_all$group <- gsub("[0-9]", "", data_all$variable)

# 只保留 taxonomy 列和 value 列
data_all <- data_all[, c(opts$taxonomy, "value", "group")]

# 绘图
p <- ggplot(data_all, aes_string(x = opts$taxonomy, y = "value", fill = "group")) +
  stat_boxplot(geom = "errorbar", width = 0.4, position = position_dodge(0.8)) +
  geom_boxplot(width = 0.6, alpha = 1, position = position_dodge(0.8), outlier.shape = NA) +
  mytheme +
  theme(legend.position = "top") +
  stat_compare_means(aes(group = group), method = "wilcox.test", label = "p.signif") +
  labs(x = opts$taxonomy, y = "Percentage (%)") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0.1)) +
  geom_jitter(aes(color = group),
              shape = 21, size = 1.6, alpha = 0.5,
              fill = "transparent",
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)
  ) +
  scale_fill_manual(values = c("#00a9ba", "#ec8883")) +
  scale_color_manual(values = c("#00a9ba", "#ec8883")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# 保存
ggsave(paste0(opts$output, "_compare.pdf"), p, width = opts$width, height = opts$height)
