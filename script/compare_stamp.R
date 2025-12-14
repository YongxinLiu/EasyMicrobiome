#!/usr/bin/env Rscript

# Copyright 2024-2026 Defeng Bai <baidefeng@caas.cn> Yong-Xin Liu <liuyongxin@caas.cn / metagenome@126.com>

# If used this script, please cited:
# Bai, et al. 2025. EasyMetagenome: A User‐Friendly and Flexible Pipeline for Shotgun Metagenomic Analysis in Microbiome Research. iMeta 4: e70001. https://doi.org/10.1002/imt2.70001

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 更新
# 2024/11/12: 增加利用R语言实现的stamp分析
# 2025/11/26: 规范格式，p-value矫正方法实现自定义输入

# 1.1 程序功能描述和主要步骤

# 程序功能：物种组成比较STAMP分析
# Functions: STAMP analysis

options(warn = -1) # Turn off warning


## 设置输入输出文件和参数

# 修改下面`default=`后面的文件和参数。
#
# 输入文件为特征表(Genus.txt)+分组信息(metadata.txt)
#
# 输入文件"-i", "--input"，默认'metaphlan4/Genus.txt'; 属水平特征表
#
# 实验设计"-d", "--metadata"，默认`metadata.txt`，可手动修改文件位置；
#
# 分组列名"-n", "--group"，默认将metadata.txt中的Group列作为分组信息，可修改为任意列名；
#
# 分组列名"-c", "--compare"，默认将比较metadata.txt中的Group列的前两个值，建议手动设定；
#
# 分组列名"-m", "--method"，默认'wilcox.test'；
#
# 分组列名"-p", "--pvalue"，默认'0.05',可手动设定；
# 
# 分组列名"-f", "--fdr"，默认'none',可手动设定pvalue矫正的方法(holm, hochberg, hommel, bonferroni, BH, BY, fdr, none)；
#
# 分组列名"-t", "--threhold"，丰度筛选阈值，默认千分之1；
#
# 分组列名"-o", "--output"，结果输出文件夹，默认'metaphlan4/'；
#
# 图片宽"-w", "--width"，默认89 mm，根据图像布局可适当增大或缩小
#
# 图片高"-e", "--height"，默认59 mm，根据图像布局可适当增大或缩小


# 1.2 依赖包安装

site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
a = rownames(installed.packages())

# install CRAN
install_CRAN <- c("tidyverse", "ggplot2", "BiocManager", "pheatmap", "dplyr",
                  "patchwork", "purrr", "broom", "tidyr", "cli", "devtools")
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
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T)
}

# 解析参数-h显示帮助信息
if (TRUE){
    option_list=list(
        # 原始OTU表counts值
        make_option(c("-i", "--input"), type="character", default="metaphlan4/Genus.txt",
                    help="OTU table in counts;  [default %default]"),
        # 元数据/实验设计文件
        make_option(c("-d", "--metadata"), type="character", default="metadata.txt",
                    help="metadata file;  [default %default]"),
        # 分组列名     
        make_option(c("-n", "--group"), type="character", default="Group",
                    help="Group name;  [default %default]"),
        # 组间比较
        make_option(c("-c", "--compare"), type="character", default="Centenarians-Young",
                    help="Groups comparison;  [default %default]"),
        # 组间比较方法
        make_option(c("-m", "--method"), type="character", default="wilcox.test",
                    help="Compare method, default t.test, alternative wilcox [default %default]"),
        # 显著性阈值
        make_option(c("-p", "--pvalue"), type="numeric", default=0.05,
                    help="Threshold of P-value, [default %default]"),
        # adjust method
        make_option(c("-f", "--fdr"), type="character", default="none",
                    help="adjust methods: holm, hochberg, hommel, bonferroni, BH, BY, fdr, none [default %default]"),
        # 相对丰度，默认千一
        make_option(c("-t", "--threshold"), type="numeric", default=0.001,
                    help="Relative abundance,  [default %default]"),
        make_option(c("-o", "--output"), type="character", default="metaphlan4/",
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
suppressWarnings(dir.create(dirname(opts$output), showWarnings = F))


# 2. 依赖关系检查、安装和加载

# 依赖包列表
package_list <- c(
  "tidyverse", "ggplot2", "BiocManager", "pheatmap", "dplyr",
  "patchwork", "purrr", "broom", "tidyr", "cli", "devtools"
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


# 3.读取输入文件

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

# 差异比较组筛选
group_list=strsplit(opts$compare,'-')[[1]]
idx=metadata$group %in% group_list
table(idx)
sub_metadata=metadata[idx,,drop=F]
sub_dat=as.matrix(dat[, rownames(sub_metadata)])


# 4.统计分析

# 数据准备
data <- t(sub_dat)
data1 <- data.frame(data, sub_metadata$group)
colnames(data1) <- c(colnames(data),"Group")
data1$Group <- as.factor(data1$Group)

# t-test
# if (opts$method == "t.test"){
# diff <- data1 %>% 
#     select_if(is.numeric) %>%
#     map_df(~ broom::tidy(t.test(. ~ Group,data=data1)), .id='var')
# # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
# diff$p.value <- p.adjust(diff$p.value, opts$fdr)
# diff <- diff %>% filter(p.value < opts$pvalue)
# }

# diff1 <- data1 %>%
#   select_if(is.numeric) %>%
#   map_df(~ broom::tidy(wilcox.test(. ~ Group,data=data1)), .id='var')
# 
# diff1$p.value <- p.adjust(diff1$p.value, opts$fdr)
# diff1 <- diff1 %>% dplyr::filter(p.value < opts$pvalue)

# diff1 <- data1 %>%
#   select_if(is.numeric) %>%
#   map_df(~ broom::tidy(wilcox.test(. ~ Group, data = data1, conf.int = TRUE, conf.level = 0.95)),
#          .id = "var")

# 批量差异分析-wilcox方法
diff1 <- data1 %>%
  select_if(is.numeric) %>%
  map_df(function(x) {
    tryCatch(
      broom::tidy(wilcox.test(x ~ data1$Group, conf.int = TRUE, conf.level = 0.95)),
      error = function(e) tibble(
        statistic = NA,
        p.value = NA,
        estimate = NA,
        conf.low = NA,
        conf.high = NA,
        method = NA,
        alternative = NA
      )
    )
  }, .id = "var")

# 多重检验校正
diff1$p.value <- p.adjust(diff1$p.value, method = opts$fdr)

# 筛选显著结果
diff1 <- diff1 %>% filter(p.value < opts$pvalue)


# 5. 绘图

# 左侧条形图
abun.bar <- data1[,c(diff1$var,"Group")] %>% 
    gather(variable,value,-Group) %>% 
    group_by(variable,Group) %>% 
    dplyr::summarise(Mean=mean(value))

# 右侧散点图
diff.mean <- diff1[,c("var","estimate","conf.low","conf.high","p.value")]
diff.mean$Group <- c(ifelse(diff.mean$estimate >0,levels(data1$Group)[1],
                            levels(data1$Group)[2]))
diff.mean <- diff.mean[order(diff.mean$estimate,decreasing=TRUE),]

# 左侧条形图
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


# 右侧散点图
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

# 图像拼接
(p <- p1 + p2 + p3 + plot_layout(widths=c(4,6,2)))

# Saving figure
# 保存图片，大家可以修改图片名称和位置，长宽单位为毫米
ggsave(
  filename = paste0(opts$output, "_STAMP.pdf"),
  plot = p,                                   
  width = opts$width,                        
  height = opts$height,                      
  units = "mm",                            
  dpi = 300,                                
  device = cairo_pdf                        
)

