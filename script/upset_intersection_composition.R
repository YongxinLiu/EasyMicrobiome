#!/usr/bin/env Rscript

# Copyright 2026 De-feng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, et al. 2021. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12: 315-330. https://doi.org/10.1007/s13238-020-00724-8

# 手动运行脚本，使用 Ctrl+Shift+H 或 Session 需要设置工作目录
# Set Work Directory - Choose Directory / To Source File Location

# Clean enviroment object
rm(list=ls()) 

# 1.1 简介 Introduction #----

# 程序功能：功能通路注释结果分组比较UpSet图
# Functions: UpSet diagram for group comparison of functional pathway annotation results


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
    make_option(c("-i", "--input"), type="character", default="eggnog/data_venn2.txt",
                help="Unfiltered OTU table [default %default]"),
    make_option(c("-g", "--group"), type="character", default="Centenarians,Young",
                help="Group names separated by comma"),
    make_option(c("-o", "--output"), type="character", default="eggnog/",
                help="Output MAGs UpSet plot for different databases  [default %default]") 
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
print("You are using the following parameters:")
print(opts)



# Install related packages
if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggplot2","colorspace","RColorBrewer","UpSetR","VennDiagram",
             "formattable")) # ,"vegan"
}
# load related packages
suppressWarnings(suppressMessages(library("ggplot2")))
suppressWarnings(suppressMessages(library("colorspace")))
suppressWarnings(suppressMessages(library("RColorBrewer")))
suppressWarnings(suppressMessages(library("UpSetR")))
suppressWarnings(suppressMessages(library("VennDiagram")))
suppressWarnings(suppressMessages(library("formattable")))

# 提取交集信息函数
get_intersect_members <- function(x, ...) {
  
  library(dplyr)
  library(tibble)
  library(rlang)
  
  x <- x[, sapply(x, is.numeric)]
  
  x <- rownames_to_column(x, var = "ID")
  
  target_groups <- c(...)
  
  cond_list <- lapply(names(x)[-1], function(item){
    
    if(item %in% target_groups){
      expr(!!sym(item) >= 1)
    } else {
      expr(!!sym(item) == 0)
    }
    
  })
  
  x %>%
    filter(!!!cond_list) %>%
    column_to_rownames(var = "ID")
}


# 载入数据
Data4Pic = read.table(opts$input, header=T, row.names=1)
Data4Pic[Data4Pic>0] = 1

# 解析分组
target_groups <- unlist(
  strsplit(opts$group, ",")
)

print("Target groups:")
print(target_groups)

# 求交集
Intersect <- do.call(
  get_intersect_members,
  c(list(Data4Pic), as.list(target_groups))
)

# 提取交集行名
Intersect_ID <- rownames(Intersect)

# 从原始数据提取
data <- Data4Pic[Intersect_ID, ]

# 组间交集提取
write.csv(Intersect, file=paste(opts$output, "/Cogs_01.csv", sep=""))

# 提取COG分类（行名）
cog <- rownames(data)

# 单字母COG
single_cog <- cog[nchar(cog) == 1]

# 标准COG顺序
cog_levels <- c(
  "E","D","C","T","G","M",
  "Z","B","K","N","O","P",
  "A","J","L","Q"
)

# 统计单字母频数
tab <- table(factor(single_cog, levels = cog_levels))

# Other
other_n <- sum(nchar(cog) > 1)

# 合并
value <- c(as.numeric(tab), other_n)

# 分组名称
group1 <- c(cog_levels, "Other")

# 构建结果表
result <- data.frame(
  COG = group1,
  value = value
)

# 计算比例
result$Rel <- result$value / sum(result$value)

# 百分比
result$Per <- sprintf("%.2f%%", result$Rel * 100)

# group1列
result$group1 <- result$COG

# 输出
print(result)

data1 <- result

# COGs饼图
# data1 <- read.table(opts$composition, header = 1, check.names = F, sep = "\t")
data1 <- as.data.frame(data1)

data1$group1 <- factor(data1$COG, levels = c("E","D","C","T","G","M",
                                             "Z","B","K","N","O","P",
                                             "A","J","L","Q","Other"
))

#确定标签位置
#Determine label position
data1$ymax<-cumsum(data1$Rel)
data1$ymin<-c(0,head(data1$ymax,n=-1))
data1$labelposition<-(data1$ymax + data1$ymin)/2

#绘制圆环图
#Draw a circular diagram
p1 <- ggplot(data1,aes(ymax=ymax,ymin=ymin,
                       xmax=3,xmin=2))+
  geom_rect(aes(fill=group1))+
  geom_text(x=2.5,aes(y=labelposition,label=paste0(group1,"\n(",Per,")")),size=4, color = "black")+
  xlim(1,3)+
  coord_polar(theta="y")+
  theme_void()+
  #guides(color=guide_legend(nrow=4,byrow=TRUE))+
  #theme(legend.position = "none")+
  scale_fill_manual(values = c("#f94e54","#5196d5","#00ceff","#ff630d","#35978b",
                               "#e5acd7","#77aecd","#ec8181","#dfc6a5","#e50719",
                               "#d27e43","#8a4984","#fe5094","#8d342e","#d2da93",
                               "#ffad00","lightgrey","#00fc8d","#b64aa0","#9b82e1"))
p1
ggsave(file=paste(opts$output, "Cogs_pie_example.pdf", sep=""), plot = p1, width = 6, height = 6)