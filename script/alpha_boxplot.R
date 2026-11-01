#!/usr/bin/env Rscript

# Copyright 2016-2026 Yong-Xin Liu <liuyongxin@caas.cn / metagenome@126.com>

# If used this script, please cited:
# Bai, et al. 2025. EasyMetagenome: A User‐Friendly and Flexible Pipeline for Shotgun Metagenomic Analysis in Microbiome Research. iMeta 4: e70001. https://doi.org/10.1002/imt2.70001

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 更新
# 2021/5/31: 更新引文，添加数据表转置和标准化的选项
# 2025/11/25: 更新引文，更新模板

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


# 1.2 依赖包安装

site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
a = rownames(installed.packages())

# install CRAN
install_CRAN <- c("ggplot2", "BiocManager", "optparse")
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



install_bioc <- c("ggplot2", "multcompView")

for (i in install_bioc) {
  if (!i %in% a)
    BiocManager::install(i, update = F, site_repository=site)
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
    make_option(c("-i", "--input"), type="character", default="result/alpha/vegan.txt",
                help="Alpha diversity matrix [default %default]"),
    make_option(c("-a", "--alpha_index"), type="character", default="richness",
                help="Group name [default %default]"),
    make_option(c("-d", "--design"), type="character", default="result/metadata.txt",
                help="Design file or metadata [default %default]"),
    make_option(c("-t", "--transpose"), type="logical", default=FALSE,
                help="Design file or metadata [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Group",
                help="Group name [default %default]"),
    make_option(c("-s", "--scale"), type="logical", default=FALSE,
                help="Normalize to 100 [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/alpha/",
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


# 2. 依赖关系检查、安装和加载
#install.packages("remotes")
#remotes::install_github("microbiota/amplicon")
suppressWarnings(suppressMessages(library(amplicon)))


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

alpha_boxplot <- function(alpha_div, metadata, index = "richness", groupID = "Group",
                          outlier = TRUE
) {
  # 依赖关系检测与安装
  p_list = c("ggplot2", "dplyr", "multcompView") # "agricolae"
  for(p in p_list){
    if (!requireNamespace(p)){
      install.packages(p)}
    suppressPackageStartupMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
  }

  # 测试默认参数
  # library(amplicon)
  # index = "richness"
  # groupID = "Group"
  # metadata = read.table("metadata2.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
  # metadata = subset(metadata, Group %in% c("KO","OE"))

  # 交叉筛选
  idx = rownames(metadata) %in% rownames(alpha_div)
  metadata = metadata[idx,,drop=F]
  alpha_div = alpha_div[rownames(metadata),]

  # 提取样品组信息,默认为group可指定
  sampFile = as.data.frame(metadata[, groupID],row.names = row.names(metadata))
  # colnames(sampFile)[1] = "group"

  # 合并alpha_div和metadata
  df = cbind(alpha_div[rownames(sampFile),index], sampFile)
  colnames(df) = c(index,"group")

  # 统计各种显著性
  model = aov(df[[index]] ~ group, data=df)
  # 计算Tukey显著性差异检验
  Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
  # 提取比较结果
  Tukey_HSD_table = as.data.frame(Tukey_HSD$group)

  # 保存统计结果
  # 保存一个制表符，解决存在行名时，列名无法对齐的问题
  write.table(paste(date(), "\nGroup\t", groupID, "\n\t", sep=""), file=paste("alpha_boxplot_TukeyHSD.txt",sep=""),append = T, quote = F, eol = "", row.names = F, col.names = F)
  # 保存统计结果，有warning正常
  suppressWarnings(write.table(Tukey_HSD_table, file=paste("alpha_boxplot_TukeyHSD.txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

  # 函数：将Tukey检验结果P值转换为显著字母分组
  # 输入文件为图基检验结果和分组
  generate_label_df = function(TUKEY, variable){
    # library(multcompView)
    # 转换P值为字母分组
    ## 提取图基检验中分组子表的第4列P adjust值
    Tukey.levels = TUKEY[[variable]][,4]
    # 方法1.multcompLetters函数将两两p值转换为字母，data.frame并生成列名为Letters的数据框
    Tukey.labels = data.frame(multcompLetters(Tukey.levels)['Letters'])
    # 方法2. 解决字母顺序相反的问题
    # library(multcomp)
    # tuk <- cld(glht(model, alternative = 'two.sided', linfct = mcp(group = 'Tukey')), sig = p, decreasing = TRUE)
    # Tukey.labels <- data.frame(Letters=tuk$mcletters$Letters, stringsAsFactors = FALSE)

    # 按分组名字母顺序
    ## 提取字母分组行名为group组名
    Tukey.labels$group = rownames(Tukey.labels)
    # 按组名的字母顺序排列，默认的Levels
    Tukey.labels=Tukey.labels[order(Tukey.labels$group), ]
    return(Tukey.labels)
  }

  # 当只有两组时，用LSD标注字母
  if (length(unique(df$group)) == 2){
    # LSD检验，添加差异组字母
    library(agricolae)
    out = LSD.test(model, "group", p.adj="none")
    stat = out$groups
    # 分组结果添入Index
    df$stat=stat[as.character(df$group),]$groups
    # 当大于两组时，用multcompView标注字母
  }else{
    # library(multcompView)
    LABELS = generate_label_df(Tukey_HSD , "group")
    df$stat=LABELS[as.character(df$group),]$Letters
  }

  # 设置分组位置为各组y最大值+高的5%
  max=max(df[,c(index)])
  min=min(df[,index])
  x = df[,c("group",index)]
  #y = x %>% group_by(group) %>% summarise(Max=paste('max(',index,')',sep=""))
  y = x %>% group_by(group) %>% summarise(Max = max(.data[[index]], na.rm = TRUE))
  y=as.data.frame(y)
  rownames(y)=y$group
  df$y=y[as.character(df$group),]$Max + (max-min)*0.05


  if (outlier) {
    # 绘图 plotting
    p = ggplot(df, aes(x=group, y=.data[[index]], color=group)) +
      geom_boxplot(alpha=1,
                   # outlier.shape = NA,
                   # outlier.size=0,
                   size=0.7,
                   width=0.5, fill="transparent") +
      labs(x="Groups", y=paste(index, "index"), color=groupID) + theme_classic() +
      geom_text(data=df, aes(x=group, y=y, color=group, label=stat)) +
      geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)+
      theme(text=element_text(family="sans", size=7))
    p
  } else{
    # 绘图 plotting
    p = ggplot(df, aes(x=group, y=.data[[index]], color=group)) +
      geom_boxplot(alpha=1,
                   outlier.shape = NA,
                   outlier.size=0,
                   size=0.7,
                   width=0.5, fill="transparent") +
      labs(x="Groups", y=paste(index, "index"), color=groupID) + theme_classic() +
      geom_text(data=df, aes(x=group, y=y, color=group, label=stat)) +
      geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)+
      theme(text=element_text(family="sans", size=7))
    p
  }


}

p = alpha_boxplot(alpha_div, index = opts$alpha_index, metadata, groupID = opts$group)
if (opts$xlabAngle){
  p = p + theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
}
# Saving figure
# 保存图片，大家可以修改图片名称和位置，长宽单位为毫米
ggsave(paste0(opts$output,"boxplot_",opts$alpha_index,".pdf"), p, width = opts$width, height = opts$height, units = "mm")
