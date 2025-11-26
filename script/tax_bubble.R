#!/usr/bin/env Rscript

# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# 如果使用本脚本，请参考：
# 提供文献引用等

# 设置命令行参数
library(optparse)
library(ggplot2)
library(reshape2)

# 解析命令行参数
option_list = list(
  make_option(c("-i", "--input"), type="character", default="./result2/tax/sum_g.txt",
              help="输入物种丰度表文件 [default %default]"),
  make_option(c("-g", "--group"), type="character", default="../metadata.txt",
              help="分组信息文件 [default %default]"),
  make_option(c("-c", "--group_col"), type="character", default="Group", 
              help="在metadata文件中指定用于分组的列名 [default %default]"),
  make_option(c("-o", "--output"), type="character", default="top20_bubble.pdf",
              help="输出PDF图片文件名 [default %default]"),
  make_option(c("-w", "--width"), type="numeric", default=10,
              help="输出图片宽度 [default %default]"),
  make_option(c("-e", "--height"), type="numeric", default=3.0,
              help="输出图片高度 [default %default]"),
  make_option(c("-s", "--size_range"), type="numeric", default=15,
              help="气泡大小范围的最大值 [default %default]"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE, 
              help="显示版本信息"),
  make_option(c("-n", "--top_n"), type="integer", default=20,
              help="指定丰度排名前N的物种数量 [default %default]")
)

# 解析命令行参数
opts = parse_args(OptionParser(option_list=option_list))

# 处理版本信息
if (opts$version) {
  cat("Version: 1.0.0\n")
  quit(status=0)
}

# 设置工作路径
#setwd(dirname(opts$input))

# 加载种水平物种丰度表，并设置列名和分隔符
data <- read.table(opts$input, header=TRUE, sep="\t", row.names=1)

# 加载分组表group
group <- read.delim(opts$group, sep='\t', stringsAsFactors=FALSE)

# 计算相对丰度（按列归一化，即每个样本中物种的百分比）
rel_abund <- data / colSums(data) * 100

# 按物种的全局总丰度排序（行和降序）
data <- data[order(rowSums(data), decreasing = TRUE), ]
top_species <- rel_abund[1:opts$top_n, ]


# 计算Others的丰度（每个样本中剩余物种的总和）
others <- colSums(rel_abund) - colSums(top_species)
others <- pmax(others, 0)  # 避免负值（理论上不应出现）
others_df <- t(data.frame(Others = others))

# 合并前N物种和Others
data2 <- rbind(top_species, others_df)
# 删除最后一列all列
data2 <- data2[, -ncol(data2)]
# 构建数据框
group_subset <- group[, c("SampleID", "Group")]  # 只保留 SampleID 和 Group
data3 <- cbind(group_subset, t(data2))          # 合并时只带入关键列
data3 <- as.data.frame(data3)

# 数据宽转长
data4 <- melt(data3)

# 根据用户指定的列名获取分组信息
data4$group <- factor(data4[[opts$group_col]], levels=unique(group[[opts$group_col]]))

# 绘图
p <- ggplot(data4, aes(x=SampleID, y=variable, size=value, colour=Group)) +
  geom_point(alpha = 0.7) +
  scale_colour_manual(values=c("#c74732", "#1EB5B8", "#A07DB7", "#8DD3C7")) +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(size="Relative Abundance(%)", x='Samples', y='Taxonomy') +
  scale_size(range=c(0, opts$size_range), breaks=c(0, 10, 20, 40), limits=c(0, 400)) +
  guides(colour=guide_legend(override.aes=list(size=8))) +
  theme(axis.text=element_text(colour='black', size=7)) +
  facet_wrap(~group, ncol=3, scales="free_x")

# 保存绘图为PDF格式
ggsave(opts$output, p, width=opts$width, height=opts$height, device="pdf")
