#!/usr/bin/env Rscript

# Copyright 2024-2026 Defeng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Bai, et al. 2025. EasyMetagenome: A User‐Friendly and Flexible Pipeline for Shotgun Metagenomic Analysis in Microbiome Research. iMeta 4: e70001. https://doi.org/10.1002/imt2.70001

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 更新
# 2024/11/12：增加按丰度排序堆叠柱状图功能
# 2025/11/27：规范代码
# 2026/5/22：解决代码报错，规范图例显示
# 2026/11/16：合并 unclassified 和 others，并将 others 固定放在堆叠柱底部


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
                  "ggprism", "dplyr", "plyr", "cowplot")
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
  "ggprism", "dplyr", "plyr", "cowplot", "scales", "grid"
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

# 将 Taxonomy 转换为字符型，避免因子类型影响后续名称替换
 taxonomy$Taxonomy <- as.character(taxonomy$Taxonomy)

# 识别原始分类表中名称包含 unclassified 或 others 的分类
# 不区分大小写，例如：s__unclassified、g__unclassified、s__others、others
 is_other <- !is.na(taxonomy$Taxonomy) &
   grepl("unclassified|others", taxonomy$Taxonomy, ignore.case = TRUE)

# 将上述分类统一重命名为 others
 taxonomy$Taxonomy[is_other] <- "others"

# 按 Taxonomy 合并重复物种，同时将原有 unclassified 和 others 的丰度相加
 data <- aggregate(. ~ Taxonomy, data = taxonomy, sum)

# 设置行名，并删除 Taxonomy 列
 rownames(data) <- data$Taxonomy
 data <- data[, -1, drop = FALSE]

# 相对丰度标准化：每个样本单独除以该样本的丰度总和
# 对总丰度为 0 的样本保留 0，避免产生 NaN 或 Inf
 data_matrix <- as.matrix(data)
 sample_sum <- colSums(data_matrix, na.rm = TRUE)
 data_matrix <- sweep(data_matrix, 2, sample_sum, FUN = "/")
 data_matrix[, sample_sum == 0] <- 0
 data <- as.data.frame(data_matrix, check.names = FALSE)

# =========================================
# 3. 丰度排序及样本排序
# =========================================

# 按各物种在所有样本中的总相对丰度降序排列
# 已经合并的 others 不参与 Top N 物种排名
all_taxa <- rownames(data)
known_others <- "others" %in% all_taxa
candidate_taxa <- setdiff(all_taxa, "others")

if (length(candidate_taxa) > 0) {
  candidate_taxa <- candidate_taxa[
    order(rowSums(data[candidate_taxa, , drop = FALSE]), decreasing = TRUE)
  ]
}

# =========================================
# 4. Top N 物种处理
# =========================================

# opts$legend 表示图例中最多显示的分类数，包含 others
# 因此保留 opts$legend - 1 个丰度最高的物种，其余全部并入 others
# 例如 opts$legend = 6 时，显示 Top 5 物种 + others
 top_n <- max(as.integer(opts$legend) - 1L, 0L)
 top_n <- min(top_n, length(candidate_taxa))

if (top_n > 0) {
  top_taxa <- candidate_taxa[seq_len(top_n)]
} else {
  top_taxa <- character(0)
}

# 未进入 Top N 的物种将与原有 others 合并
remaining_taxa <- setdiff(candidate_taxa, top_taxa)
need_others <- known_others || length(remaining_taxa) > 0

if (need_others) {
  # 初始化 others 在每个样本中的丰度
  other_values <- setNames(rep(0, ncol(data)), colnames(data))

  # 加入输入文件中原本就存在的 others（其中已包含 unclassified）
  if (known_others) {
    other_values <- other_values + as.numeric(data["others", , drop = TRUE])
    names(other_values) <- colnames(data)
  }

  # 加入未进入 Top N 的其余物种丰度
  if (length(remaining_taxa) > 0) {
    other_values <- other_values + colSums(
      data[remaining_taxa, , drop = FALSE]
    )
  }

  # 构建 others 行，并放在数据矩阵第一行
  others_df <- as.data.frame(as.list(other_values), check.names = FALSE)
  rownames(others_df) <- "others"

  if (length(top_taxa) > 0) {
    data_top <- rbind(others_df, data[top_taxa, , drop = FALSE])
  } else {
    data_top <- others_df
  }
} else {
  # 当没有 others 且所有物种均进入 Top N 时，不额外添加 others
  data_top <- data[top_taxa, , drop = FALSE]
}

# 按丰度最高的物种对样本排序，保持原脚本的样本排序逻辑
if (length(top_taxa) > 0) {
  sample_order <- colnames(data_top)[
    order(-as.numeric(data_top[top_taxa[1], ]))
  ]
} else {
  sample_order <- colnames(data_top)[
    order(-as.numeric(data_top["others", ]))
  ]
}
data_top <- data_top[, sample_order, drop = FALSE]

# 明确设置分类顺序：others 第一，其余物种按总丰度从高到低排列
# 后续绘图将配合 position_stack(reverse = TRUE)，使 others 位于柱子底部
tax_levels <- rownames(data_top)

# =========================================
# 5. 转长表
# =========================================

# 添加 taxonomy 名称，并将矩阵转换为长表
mean_sort3 <- data_top
mean_sort3$tax <- rownames(mean_sort3)

data_all <- reshape2::melt(mean_sort3, id.vars = "tax")
data_all <- as.data.frame(data_all)

# 固定分类因子顺序：others 放在第一位，后续物种按总丰度排序
# 配合 position_stack(reverse = TRUE)，确保 others 位于堆叠柱底部
data_all$tax <- factor(data_all$tax, levels = tax_levels)

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
      fill = tax
    )
  ) +
    geom_bar(
      stat = "identity",
      position = position_stack(reverse = TRUE),
      width = 1
    ) +
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
    ) +
    # 图例从上到下按堆叠柱由顶部到底部排列，others 放在图例最下面
    guides(fill = guide_legend(reverse = TRUE))
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

# 再次明确 tax 的因子顺序，确保主图和图例顺序一致
data_all2$tax <- factor(as.character(data_all2$tax), levels = tax_levels)

# 单独绘制 legend；使用与主图相同的堆叠顺序
legend_plot <- ggplot(
  data_all2,
  aes(
    x = variable,
    y = value,
    fill = tax
  )
) +
  geom_bar(
    stat = "identity",
    position = position_stack(reverse = TRUE)
  ) +
  scale_fill_manual(
    values = c(
      "#e5acd7", "#00ceff", "#ff630d", "#35978b", "#d2da93",
      "#5196d5", "#77aecd", "#ec8181", "#dfc6a5", "#e50719",
      "#d27e43", "#8a4984", "#fe5094", "#8d342e", "#f94e54",
      "#ffad00", "#36999d", "#00fc8d", "#b64aa0", "#9b82e1"
    ),
    name = strsplit(basename(opts$output), "[._]")[[1]][1]
  ) +
  # 图例从上到下按堆叠柱由顶部到底部排列，others 放在图例最下面
  guides(fill = guide_legend(reverse = TRUE)) +
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
