#!/usr/bin/env Rscript

# Copyright 2016-2022 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo & Yang Bai. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein Cell 41, 1-16, doi:10.1007/s13238-020-00724-8 (2020).

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
# 分组列名"-t", "--threhold"，丰度筛选阈值，默认千分之1
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
        # 原始OTU表counts值
        make_option(c("-i", "--input"), type="character", default="result/otutab.txt", # otutab.txt stamp/tax_6Genus.txt
                    help="OTU table in counts;  [default %default]"),
        # 元数据/实验设计文件
        make_option(c("-d", "--metadata"), type="character", default="result/metadata.txt",
                    help="metadata file;  [default %default]"),
        # 分组列名     
        make_option(c("-n", "--group"), type="character", default="Group",
                    help="Group name;  [default %default]"),
        # 组间比较
        make_option(c("-c", "--compare"), type="character", default="KO-OE",
                    help="Groups comparison;  [default %default]"),
        # 组间比较方法
        make_option(c("-m", "--method"), type="character", default="t.test",
                    help="Compare method, default t.test, alternative wilcox [default %default]"),
        # 显著性阈值
        make_option(c("-p", "--pvalue"), type="numeric", default=0.05,
                    help="Threshold of P-value, [default %default]"),
        # adjust method
        make_option(c("-f", "--fdr"), type="character", default="none",
                    help="adjust methods: holm, hochberg, hommel, bonferroni, BH, BY, fdr, none [default %default]"),
        # 相对丰度，默认千一
        make_option(c("-t", "--threshold"), type="numeric", default=0.1,
                    help="Relative abundance,  [default %default]"),
        # holm, hochberg, hommel, bonferroni, BH, BY, fdr, none
        make_option(c("-o", "--output"), type="character", default="result/stamp/",
                    help="Output prefix; [default %default]"),
        # 图片宽mm
        make_option(c("-w", "--width"), type="numeric", default=89,
                    help="Figure width;  [default %default]"),
        # 图片高mm
        make_option(c("-e", "--height"), type="numeric", default=59,
                    help="Figure heidth;  [default %default]")
    )
    opts=parse_args(OptionParser(option_list=option_list))
    
    # 调置如果无调设置输出，根据其它参数设置默认输出
    if (opts$output==""){
        opts$output=paste("stamp/",opts$compare, sep="")}
    
    # 显示输入输出参数，用户确认是否正确
    print("Parameters are as follows. Please check it!")
    print(paste("The input data matrix file is ", opts$input,  sep=""))
    print(paste("The metadata file is ", opts$metadata,  sep=""))
    print(paste("Group name is ", opts$group,  sep=""))
    print(paste("Group compare is ", opts$compare,  sep=""))
    print(paste("Threshold of P-value is ", opts$pvalue,  sep=""))
    print(paste("Threshold of FDR is ", opts$fdr,  sep=""))
    print(paste("Threshold of relative abundance is ", opts$threshold,  sep=""))
    print(paste("Output figure width ", opts$width,  sep=""))
    print(paste("Output figure height ", opts$height,  sep=""))
    print(paste("The output file is ", opts$output, sep=""))
}
# 设置输出文件缺省值，如果为空，则为输入+pcoa.pdf
suppressWarnings(dir.create(dirname(opts$output), showWarnings = F))


#----1.3. 加载包 Load packages#----
# suppressWarnings(suppressMessages(library(amplicon)))
# 依赖包列表：差异分析、绘图、热图、数据变换和开发者工具
package_list=c("tidyverse", "ggplot2","BiocManager","pheatmap","dplyr","devtools")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
    if(!suppressWarnings(suppressMessages(require(p, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)))){
        install.packages(p, repos=site)
        suppressWarnings(suppressMessages(library(p, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)))
    }
}


#----2. 读取文件 Read files#----

# 读取OTU表
dat=read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char="") 

# 读取实验设计
metadata=read.table(opts$metadata, header=T, row.names= 1, sep="\t", comment.char="") 

# 将选定的分组列统一命名为group
metadata$group=metadata[,opts$group]

# 标准化和按丰度筛选
idx=rownames(metadata) %in% colnames(dat)
table(idx)
metadata=metadata[idx,,drop=F]
dat=dat[, rownames(metadata)]

#----丰度过滤#----
# 标准化为百分比
if (TRUE){ # normalize
    norm=t(t(dat)/colSums(dat,na=T)*100)
}else{
    norm=as.matrix(data)
}
# 按丰度筛选标准化特征表和原始值
idx=rowMeans(norm) > opts$threshold
norm=norm[idx, ]
colSums(norm)
dat=dat[idx, ]

#----差异比较组筛选#----
group_list=strsplit(opts$compare,'-')[[1]]
idx=metadata$group %in% group_list
table(idx)
sub_metadata=metadata[idx,,drop=F]
sub_dat=as.matrix(dat[, rownames(sub_metadata)])


#----3. 统计保存 Stat and saving#----

#----3.1 统计 Stat#----
data <- t(sub_dat)
data1 <- data.frame(data, sub_metadata$group)
colnames(data1) <- c(colnames(data),"Group")
data1$Group <- as.factor(data1$Group)

## t-test

if (opts$method == "t.test"){
diff <- data1 %>% 
    select_if(is.numeric) %>%
    map_df(~ broom::tidy(t.test(. ~ Group,data=data1)), .id='var')
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
diff$p.value <- p.adjust(diff$p.value, opts$fdr)
diff <- diff %>% filter(p.value < opts$pvalue)
}

## wilcox
if (opts$method == "wilcox"){
    diff1 <- data1 %>%
    select_if(is.numeric) %>%
    map_df(~ broom::tidy(wilcox.test(. ~ Group,data=data1)), .id='var')

diff1$p.value <- p.adjust(diff1$p.value, opts$fdr)
diff1 <- diff %>% filter(p.value < opts$pvalue)
}

## 绘图数据构建
## 左侧条形图
abun.bar <- data1[,c(diff$var,"Group")] %>% 
    gather(variable,value,-Group) %>% 
    group_by(variable,Group) %>% 
    summarise(Mean=mean(value))

## 右侧散点图
diff.mean <- diff[,c("var","estimate","conf.low","conf.high","p.value")]
diff.mean$Group <- c(ifelse(diff.mean$estimate >0,levels(data1$Group)[1],
                            levels(data1$Group)[2]))
diff.mean <- diff.mean[order(diff.mean$estimate,decreasing=TRUE),]

## 左侧条形图
library(ggplot2)
cbbPalette <- c("#E69F00", "#56B4E9")
abun.bar$variable <- factor(abun.bar$variable,levels=rev(diff.mean$var))
p1 <- ggplot(abun.bar,aes(variable,Mean,fill=Group)) +
    scale_x_discrete(limits=levels(diff.mean$var)) +
    coord_flip() +
    xlab("") +
    ylab("Mean proportion (%)") +
    theme(panel.background=element_rect(fill='transparent'),
          panel.grid=element_blank(),
          axis.ticks.length=unit(0.4,"lines"), 
          axis.ticks=element_line(color='black'),
          axis.line=element_line(colour="black"),
          axis.title.x=element_text(colour='black', size=12,face="bold"),
          axis.text=element_text(colour='black',size=10,face="bold"),
          legend.title=element_blank(),
          legend.text=element_text(size=12,face="bold",colour="black",
                                   margin=margin(r=20)),
          legend.position=c(-1,-0.1),
          legend.direction="horizontal",
          legend.key.width=unit(0.8,"cm"),
          legend.key.height=unit(0.5,"cm"))


for (i in 1:(nrow(diff.mean) - 1)) 
    p1 <- p1 + annotate('rect', xmin=i+0.5, xmax=i+1.5, ymin=-Inf, ymax=Inf, 
                        fill=ifelse(i %% 2 == 0, 'white', 'gray95'))

p1 <- p1 + 
    geom_bar(stat="identity",position="dodge",width=0.7,colour="black") +
    scale_fill_manual(values=cbbPalette)


## 右侧散点图
diff.mean$var <- factor(diff.mean$var,levels=levels(abun.bar$variable))
diff.mean$p.value <- signif(diff.mean$p.value,3)
diff.mean$p.value <- as.character(diff.mean$p.value)
p2 <- ggplot(diff.mean,aes(var,estimate,fill=Group)) +
    theme(panel.background=element_rect(fill='transparent'),
          panel.grid=element_blank(),
          axis.ticks.length=unit(0.4,"lines"), 
          axis.ticks=element_line(color='black'),
          axis.line=element_line(colour="black"),
          axis.title.x=element_text(colour='black', size=12,face="bold"),
          axis.text=element_text(colour='black',size=10,face="bold"),
          axis.text.y=element_blank(),
          legend.position="none",
          axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.title=element_text(size=15,face="bold",colour="black",hjust=0.5)) +
    scale_x_discrete(limits=levels(diff.mean$var)) +
    coord_flip() +
    xlab("") +
    ylab("Difference in mean proportions (%)") +
    labs(title="95% confidence intervals") 

for (i in 1:(nrow(diff.mean) - 1)) 
    p2 <- p2 + annotate('rect', xmin=i+0.5, xmax=i+1.5, ymin=-Inf, ymax=Inf, 
                        fill=ifelse(i %% 2 == 0, 'white', 'gray95'))

p2 <- p2 +
    geom_errorbar(aes(ymin=conf.low, ymax=conf.high), 
                  position=position_dodge(0.8), width=0.5, size=0.5) +
    geom_point(shape=21,size=3) +
    scale_fill_manual(values=cbbPalette) +
    geom_hline(aes(yintercept=0), linetype='dashed', color='black')


p3 <- ggplot(diff.mean,aes(var,estimate,fill=Group)) +
    geom_text(aes(y=0,x=var),label=diff.mean$p.value,
              hjust=0,fontface="bold",inherit.aes=FALSE,size=3) +
    geom_text(aes(x=nrow(diff.mean)/2 +0.5,y=0.85),label="P-value (corrected)",
              srt=90,fontface="bold",size=5) +
    coord_flip() +
    ylim(c(0,1)) +
    theme(panel.background=element_blank(),
          panel.grid=element_blank(),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          axis.text=element_blank(),
          axis.title=element_blank())

#----3.2 保存表格 Saving#----
## 图像拼接
library(patchwork)
(p <- p1 + p2 + p3 + plot_layout(widths=c(4,6,2)))

## 保存图像
ggsave(paste0(opts$output, "_", gsub(".txt", "", basename(opts$input)), ".pdf"), p, width=opts$width*2, height=opts$height*dim(diff)[1]/10, units="mm")
#----3.2 保存表格 Saving#----
filename=paste0(opts$output, "_", basename(opts$input))
write.table(diff, file=filename, append=F, quote=F, sep='\t', row.names=F, col.names=T)
