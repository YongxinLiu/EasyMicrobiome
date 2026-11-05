#!/usr/bin/env Rscript

# 加载包
suppressPackageStartupMessages({
  library(optparse)
  library(VennDiagram)
  library(grid)
})

# 定义命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "输入文件路径，如：otu_group_exist.txt"),
  make_option(c("-g", "--groups"), type = "character", help = "分组名称，用逗号分隔，例如 feces,plaque,saliva"),
  make_option(c("-o", "--output"), type = "character", default = "venn.pdf", help = "输出 PDF 文件名")
)

# 解析参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 检查必要参数
if (is.null(opt$input) || is.null(opt$groups)) {
  print_help(opt_parser)
  stop("请提供输入文件路径和分组名称 (--input 和 --groups)", call. = FALSE)
}

# 拆分组
group_names <- strsplit(opt$groups, ",")[[1]]

# 检查组数限制
if (length(group_names) < 2 || length(group_names) > 4) {
  stop("目前仅支持 2 到 4 个分组的 Venn 图")
}

# 固定颜色（最多4组）
fill_colors <- c( "#1EB5B8", "#edae11","#A07DB7","#7AA82C")[1:length(group_names)]

# 读取输入文件
venn_data <- read.table(opt$input, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 构造 OTU 列表
venn_list <- list()
for (grp in group_names) {
  venn_list[[grp]] <- unique(venn_data$OTU[venn_data$Group == grp])
}

# 输出 Venn 图
pdf(opt$output, width = 6, height = 6)
venn.plot <- venn.diagram(
  x = venn_list,
  category.names = group_names,
  filename = NULL,
  fill = fill_colors,
  alpha = 0.6,
  cat.col = fill_colors,
  lwd = 0,
  col = "transparent",
  margin = 0.1
)
grid.draw(venn.plot)
dev.off()
