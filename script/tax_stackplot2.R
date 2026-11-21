#!/usr/bin/env Rscript

# 关闭警告信息
options(warn = -1)

# 设置清华源
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"

# 解析命令行参数
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE)))) {
  install.packages("optparse", repos=site)
  require("optparse", character.only = TRUE)
}

option_list = list(
  make_option(c("-i", "--input"), type="character", default="kraken2/bracken.S.txt",
              help="Taxonomy composition [default %default]"),
  make_option(c("-d", "--design"), type="character", default="metadata.txt",
              help="Design file [default %default]"),
  make_option(c("-n", "--group"), type="character", default="Group",
              help="Group name [default %default]"),
  make_option(c("-o", "--output"), type="character", default="",
              help="Output directory [default %default]"),
  make_option(c("-l", "--legend"), type="numeric", default=12,
              help="Legend number [default %default]"),
  make_option(c("-c", "--color"), type="character", default="Paired",
              help="Color palette (manual1, Paired, Set3) [default %default]"),
  make_option(c("-w", "--width"), type="numeric", default=181,
              help="Figure width (mm) [default %default]"),
  make_option(c("-e", "--height"), type="numeric", default=118,
              help="Figure height (mm) [default %default]"),
  make_option(c("-t", "--top_n"), type="integer", default=30,
              help="前N个物种 [default %default]"),
  make_option(c("-s", "--sort_order"), type="character", default=NULL, 
              help="手动指定分组顺序（逗号分隔），如 'acbdef'")
)

opts = parse_args(OptionParser(option_list=option_list))

if (opts$output == "") { opts$output = opts$input }

suppressWarnings(suppressMessages(library(amplicon)))
suppressWarnings(suppressMessages(library(ggplot2)))

# 读取元数据
metadata = read.table(opts$design, header=TRUE, row.names=1, sep="\t", comment.char="", stringsAsFactors = FALSE)

# 读取物种组成矩阵
taxonomy = read.table(opts$input, header=TRUE, row.names=1, sep="\t", comment.char="", quote = "")

# 计算物种的总丰度
taxonomy$sum_abundance <- rowSums(taxonomy)

# 按丰度降序排序
taxonomy_sorted <- taxonomy[order(taxonomy$sum_abundance, decreasing = TRUE), ]

# 选取前 N 物种（去除 sum_abundance 列）
top_n_species <- min(opts$top_n, nrow(taxonomy_sorted))
taxonomy_top_n <- taxonomy_sorted[1:top_n_species, -ncol(taxonomy_sorted)]

# 计算相对丰度（归一化）
taxonomy_rel_abundance <- sweep(taxonomy_top_n, 2, colSums(taxonomy), FUN = "/") * 100
# **自动获取当前数据中的组名**
group_levels <- unique(metadata[[opts$group]])
# **判断用户是否提供了手动排序**
if (!is.null(opts$sort_order)) {
  custom_order <- unlist(strsplit(opts$sort_order, ","))
  missing_groups <- setdiff(group_levels, custom_order)
  
  # **如果有未包含的组，按原顺序附加**
  custom_order <- c(custom_order, missing_groups)
} else {
  # **如果用户未提供顺序，默认按字母顺序**
  custom_order <- sort(group_levels)
}

# **重新设定分组顺序**
metadata[[opts$group]] <- factor(metadata[[opts$group]], levels = custom_order)

# 自定义颜色集
colorset2 = c("#d2da93","#5196d5","#00ceff","#ff630d","#35978b",
              "#e5acd7","#77aecd","#ec8181","#dfc6a5","#e50719",
              "#d27e43","#8a4984","#fe5094","#8d342e","#f94e54",
              "#CAB2D6", "#FFFF99", "#8DD3C7", "#FFED6F", "#BC80BD",
              "#ffad00","#36999d","#00fc8d","#b64aa0","#9b82e1",
              "#66c2a4","#8c96c6","#fdbb84","#a6bddb","#fa9fb5"
              )
# 设置颜色22种
colorset1 = c('#98d66d','#45b1f9','#ffa6f0','#f76fc3','#85c1b8',
              '#a584ff','#ffb444','#c45cc0','#7ebfe5','#cec0c9',
              '#467584','#005ff9','#bc8c38','#bcba6d','#91749e',
              '#b2798d','#fcef5d','#b23301','#235931',"#892e4f",
              '#fabc75','#f75c39')


colorset3 = c('#a3c9dc','#2072a9','#add488','#3a9939','#f19795',
              '#d51f1f','#ffb444','#ef7c1b','#c6b0d1','#643b90',
              '#b2798d','#005ff9','#bc8c38','#bcba6d','#91749e',
              '#b2798d','#fcef5d','#b23301','#235931',"#892e4f",
              '#fabc75','#f75c39')


# 生成样本堆叠图
p_sample = tax_stackplot(taxonomy_rel_abundance, metadata, topN = opts$top_n, groupID = opts$group, style = "sample", sorted = "abundance") +
  theme(legend.text = element_text(size = 6),
        legend.key.size = unit(4, "mm"),
        legend.spacing.y = unit(2, "mm"))+
  scale_y_continuous(labels = scales::percent, expand = c(0,0))

if (opts$color == "manual1") {
  p_sample = p_sample + scale_fill_manual(values = colorset2)
} else if (opts$color == "Paired") {
  p_sample = p_sample + scale_fill_brewer(palette = "Paired")
} else if (opts$color == "Set3") {
  p_sample = p_sample + scale_fill_brewer(palette = "Set3")
}

ggsave(paste0(opts$output, ".sample.pdf"), p_sample, width = opts$width, height = opts$height, units = "mm")

# 生成分组堆叠图
p_group = tax_stackplot(taxonomy_rel_abundance, metadata, topN = opts$top_n, groupID = opts$group, style = "group", sorted = "abundance") +
  theme(legend.text = element_text(size = 6),
        legend.key.size = unit(4, "mm"),
        legend.spacing.y = unit(2, "mm"))+
  scale_y_continuous(labels = scales::percent, expand = c(0,0))

if (opts$color == "manual1") {
  p_group = p_group + scale_fill_manual(values = colorset3)
} else if (opts$color == "Paired") {
  p_group = p_group + scale_fill_brewer(palette = "Paired")
} else if (opts$color == "Set3") {
  p_group = p_group + scale_fill_brewer(palette = "Set3")
}

ggsave(paste0(opts$output, ".group.pdf"), p_group, width = opts$width, height = opts$height, units = "mm")
