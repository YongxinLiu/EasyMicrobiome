#!/usr/bin/env Rscript

# Copyright 2016-2022 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Tong Chen, Yong-Xin Liu, Luqi Huang. 2022. ImageGP: An easy-to-use data visualization web server for scientific researchers. iMeta 1: e5. https://doi.org/10.1002/imt2.5
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, Yang Bai. 2021. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein & Cell 12: 315-330. https://doi.org/10.1007/s13238-020-00724-8


# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

#----1. 参数 Parameters#----

#----1.1 功能描述 Function description#----

# 程序功能：物种组成弦图
# Functions: Taxonomy circlize

options(warn = -1) # Turn off warning


## 设置输入输出文件和参数

# 修改下面`default=`后面的文件和参数。
#
# 输入文件为特征表(otutab.txt)+分组信息(metadata.tsv)
#
# 输入文件"-i", "--input"，otutab.txt; 特征表
#
# 实验设计"-d", "--design"，默认`metadata.tsv`，可手动修改文件位置；
#
# 分组列名"-n", "--group"，默认将metadata.tsv中的Group列作为分组信息，可修改为任意列名；
#
# 分组列名"-c", "--compare_pair"，默认将比较metadata.tsv中的Group列的前两个值，建议手动设定；
#
# 分组列名"-t", "--threshold"，丰度筛选阈值，默认千分之1
#
# 图片宽"-w", "--width"，默认89 mm，根据图像布局可适当增大或缩小
#
# 图片高"-e", "--height"，默认59 mm，根据图像布局可适当增大或缩小


#----1.2 参数缺少值 Default values#----
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
if (TRUE){
    option_list=list(
        # 差异比较结果
        make_option(c("-i", "--input"), type="character", default="result/picrust/KO-WT.txt", 
                    help="OTU table in counts;  [default %default]"),
        make_option(c("-d", "--data"), type="character", default="MeanA",
                    help="data name [default %default]"),
        # 注释文件
        make_option(c("-a", "--annotation"), type="character", default="result/picrust/pathway2.anno.txt",
                    help="metadata file;  [default %default]"),
        # holm, hochberg, hommel, bonferroni, BH, BY, fdr, none
        make_option(c("-o", "--output"), type="character", default="result/picrust/KO-WT.bar.pdf",
                    help="Output prefix; [default %default]"),
        # 图片宽mm
        make_option(c("-w", "--width"), type="numeric", default=181,
                    help="Figure width;  [default %default]"),
        # 图片高mm
        make_option(c("-e", "--height"), type="numeric", default=247,
                    help="Figure heidth;  [default %default]")
    )
    opts=parse_args(OptionParser(option_list=option_list))
    
    # 调置如果无调设置输出，根据其它参数设置默认输出
    if (opts$output==""){
        opts$output=paste(opts$input,".bar.pdf", sep="")}
    
    # 显示输入输出参数，用户确认是否正确
    print("Parameters are as follows. Please check it!")
    print(paste("The input data matrix file is ", opts$input,  sep=""))
    print(paste("Output figure width ", opts$width,  sep=""))
    print(paste("Output figure height ", opts$height,  sep=""))
    print(paste("The output file is ", opts$output, sep=""))
}
# 设置输出文件缺省值，如果为空，则为输入+pcoa.pdf
suppressWarnings(dir.create(dirname(opts$output), showWarnings = F))


#----1.3. 加载包 Load packages#----
# suppressWarnings(suppressMessages(library(amplicon)))
# 依赖包列表：差异分析、绘图、热图、数据变换和开发者工具
package_list=c("reshape2", "ggplot2")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
    if(!suppressWarnings(suppressMessages(require(p, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)))){
        install.packages(p, repos=site)
        suppressWarnings(suppressMessages(library(p, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)))
    }
}


#----2. 读取文件 Read files#----

# 读取特征表
dat=read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char="", quote="") 
# 读取注释
anno=read.table(opts$annotation, header=T, row.names= 1, sep="\t", comment.char="", quote="") 
# 添加注释
KEGG = cbind(dat, anno[rownames(dat),])

# KEGG第一层分类中名字特别长，需要自动换行(替换空格为换行\n)
swr = function(string, nwrap = 12){
    paste(strwrap(string,width = nwrap),collapse = "\n")
}
swr = Vectorize(swr)
KEGG$Category <- swr(KEGG$Category)
KEGG$KEGG_Pathways = rownames(KEGG)

# 绘制L2(KEGG_Pathways)级和丰度的柱状图，按L1(Category)着色并分面
opts$data = sym(opts$data)
p <- ggplot(KEGG,aes(!!opts$data, KEGG_Pathways)) +
    geom_bar(aes(fill = Category),stat = "identity",width = 0.6) +
    xlab("Relative abundance (%)") + 
    ylab("KEGG Pathway") +
    theme(panel.background = element_rect(fill = "white",colour='black'),
          panel.grid.major = element_line(color = "grey",linetype = "dotted",size = 0.3),
          panel.grid.minor = element_line(color = "grey",linetype = "dotted",size = 0.3),
          axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_text(colour='black', size=8,face = "bold"),
          axis.title.y=element_text(colour='black', size=8),
          axis.text.x=element_text(colour='black',size=8),
          axis.text.y = element_text(color = "black",size = 8),
          legend.position = "none",
          strip.text.y = element_text(angle = 0, size = 10, face = "bold")) +
    facet_grid(Category~.,space = "free_y",scales = "free_y")
# 预览
p

## 保存图像
ggsave(opts$output, p, width=opts$width, height=opts$height, units="mm")
