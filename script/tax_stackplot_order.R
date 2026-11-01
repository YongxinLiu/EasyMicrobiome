#!/usr/bin/env Rscript

# Copyright 2024-2026 Defeng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Bai, et al. 2025. EasyMetagenome: A User‐Friendly and Flexible Pipeline for Shotgun Metagenomic Analysis in Microbiome Research. iMeta 4: e70001. https://doi.org/10.1002/imt2.70001

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 更新
# 2024/11/12：增加按丰度排序堆叠柱状图功能
# 2025/11/27：规范代码
# 2026/5/22：解决代码报错，规范图例显示


# 1.1 程序功能描述和主要步骤

# 程序功能：排序分面堆叠柱状图
# Functions: Sort faceted stacked bar chart

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
# 图例设置"-l", "--legend"，默认`6`，同一列最多显示6个；
#
# 颜色设置"-l", "--color"，默认`Paired`，设置ggplot2绘图颜色；
#
# 图片宽"-w", "--width"，默认89 mm，根据图像布局可适当增大或缩小
#
# 图片高"-e", "--height"，默认59 mm，根据图像布局可适当增大或缩小


# 1.2 依赖包安装

site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
a = rownames(installed.packages())

# install CRAN
install_CRAN <- c("ggplot2", "BiocManager", "optparse","patchwork","reshape2","magrittr",
                  "ggprism", "dplyr", "plyr")
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
}
# 设置输出文件缺省值
if(opts$output==""){opts$output=opts$input}


# 2. 依赖关系检查、安装和加载
# 依赖包列表
package_list <- c(
  "ggplot2", "BiocManager", "optparse","patchwork","reshape2","magrittr",
  "ggprism", "dplyr", "plyr","scales","grid"
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

# =========================================
# 1. 读取输入文件
# =========================================

# 实验设计 Metadata
metadata <- read.table(
  opts$design,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  comment.char = "",
  stringsAsFactors = FALSE
)

# 物种组成矩阵 Taxonomy matrix
taxonomy <- read.table(
  opts$input,
  header = TRUE,
  sep = "\t",
  comment.char = "",
  quote = "",
  stringsAsFactors = FALSE
)

# =========================================
# 2. 数据处理
# =========================================

# 按 Taxonomy 合并重复物种
data <- aggregate(. ~ Taxonomy, data = taxonomy, sum)

# 设置行名
rownames(data) <- data$Taxonomy
data <- data[, -1]

# 相对丰度标准化
data <- apply(data, 2, function(x) x / sum(x))
data <- as.data.frame(data)

# =========================================
# 3. 丰度排序
# =========================================

# 按总丰度降序排列
mean_sort <- data[order(-rowSums(data)), ]
mean_sort <- as.data.frame(mean_sort)

# 样本按第一物种丰度排序
mean_sort2 <- t(mean_sort)
mean_sort2 <- mean_sort2[order(-mean_sort2[,1]), ]
mean_sort3 <- t(mean_sort2)
mean_sort3 <- as.data.frame(mean_sort3)

# =========================================
# 4. Top N物种处理
# =========================================

# Top N + others
other <- colSums(mean_sort3[opts$legend:nrow(mean_sort3), ])

mean_sort3 <- mean_sort3[(opts$legend - 1):1, ]
mean_sort3 <- rbind(other, mean_sort3)

rownames(mean_sort3)[1] <- "others"
mean_sort3 <- as.data.frame(mean_sort3)

# =========================================
# 5. 转长表
# =========================================

# 添加 taxonomy 名称
mean_sort3$tax <- rownames(mean_sort3)

data_all <- reshape2::melt(mean_sort3, id.vars = "tax")
data_all <- as.data.frame(data_all)

# group 提取
data_all$group <- as.character(data_all$variable)
data_all$group <- gsub("[0-9]", "", data_all$group)

# =========================================
# 6. Metadata分组映射
# =========================================

# 从 metadata 建立映射表：首字母 -> 完整组名
map <- metadata %>%
  mutate(prefix = toupper(substr(Group, 1, 1))) %>%
  group_by(prefix) %>%
  slice(1) %>%
  ungroup() %>%
  select(prefix, full_group = Group)

# 替换 group 名称
data_all <- data_all %>%
  mutate(prefix = toupper(substr(group, 1, 1))) %>%
  left_join(map, by = "prefix") %>%
  mutate(group = if_else(!is.na(full_group), full_group, group)) %>%
  select(-prefix, -full_group)

# 保持 group 顺序
data_all2 <- data_all %>%
  mutate(group = factor(group, levels = unique(group)))

# =========================================
# 7. 分面作图
# =========================================

plots <- lapply(split(data_all2, data_all2$group), function(df) {
  
  group_name <- unique(df$group)
  
  ggplot(
    df,
    aes(
      x = factor(variable, levels = unique(df$variable)),
      y = value,
      fill = factor(tax, levels = unique(df$tax))
    )
  ) +
    geom_bar(stat = "identity", position = "stack", width = 1) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    theme_classic() +
    labs(x = group_name, y = NULL) +
    scale_fill_manual(
      values = c(
        "#e5acd7", "#00ceff", "#ff630d", "#35978b", "#d2da93",
        "#5196d5", "#77aecd", "#ec8181", "#dfc6a5", "#e50719",
        "#d27e43", "#8a4984", "#fe5094", "#8d342e", "#f94e54",
        "#ffad00", "#36999d", "#00fc8d", "#b64aa0", "#9b82e1"
      ),
    ) 
})

# =========================================
# 8. 去掉后续分面的坐标轴（不处理legend）
# =========================================

if (length(plots) > 1) {
  for (i in 2:length(plots)) {
    plots[[i]] <- plots[[i]] +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none"
      )
  }
}

# 第一个分面也去掉 legend
plots[[1]] <- plots[[1]] +
  ylab("Percentage (%)") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

# =========================================
# 9. patchwork组合（legend自动宽度，不重叠）
# =========================================

sample_counts <- table(data_all2$group)
relative_widths <- as.numeric(sample_counts) / sum(sample_counts)

# 主图（无legend）
main_plot <- wrap_plots(
  plots,
  widths = relative_widths
)

# 先把 tax 转 factor
data_all2$tax <- factor(data_all2$tax, levels = unique(data_all2$tax))

# 单独legend图
legend_plot <- ggplot(
  data_all2,
  aes(
    x = variable,
    y = value,
    fill = factor(tax, levels = unique(data_all2$tax))
  )
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c(
      "#e5acd7", "#00ceff", "#ff630d", "#35978b", "#d2da93",
      "#5196d5", "#77aecd", "#ec8181", "#dfc6a5", "#e50719",
      "#d27e43", "#8a4984", "#fe5094", "#8d342e", "#f94e54",
      "#ffad00", "#36999d", "#00fc8d", "#b64aa0", "#9b82e1"
    ),
    name = strsplit(basename(opts$output), "[._]")[[1]][1]
  ) +
  theme_classic()+
  theme(
    legend.position = "right",
    legend.box.margin = margin(0, 10, 0, 10),
    # 图例标题字体
    legend.title = element_text(size = 8),
    # 图例标签字体
    legend.text = element_text(size = 6),
    # 图例方块大小
    legend.key.size = unit(0.4, "cm"),
    # 图例行间距
    legend.spacing.y = unit(0.1, "cm")
  )


# 提取legend
legend <- cowplot::get_legend(legend_plot)

# 计算legend真实宽度（单位：inch）
legend_width <- sum(convertWidth(legend$widths, "in", valueOnly = TRUE))

# 主图给固定宽度，legend按真实宽度占比
main_width <- 4   # 主图基准宽度（可调整）

p <- cowplot::plot_grid(
  main_plot,
  legend,
  nrow = 1,
  rel_widths = c(main_width, legend_width),
  align = "h"
)

# 保存 Saving 
# 大家可以修改图片名称和位置，长宽单位为毫米
ggsave(paste0(opts$output,".sample.order.pdf"), p, width = opts$width, height = opts$height, units = "mm")
