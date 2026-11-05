#rm(list=ls())
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

# 设置命令行参数
option_list <- list(
  #make_option(c("-i", "--input"), type = "character", default = "./result2/compare/saliva-plaque2.txt", help = "组间比较结果表格"),
  make_option(c("-i", "--input"), type="character", default="result2/compare/data5.txt",
              help="Taxonomy composition [default %default]"),
  make_option(c("-o", "--output"), type="character", default="result2/compare/",
              help="Output directory [default %default]"),
  make_option(c("-w", "--width"), type="numeric", default=120 * 1.5,
              help="Figure width (mm) [default %default]"),
  make_option(c("-e", "--height"), type="numeric", default=80 * 1.5,
              help="Figure height (mm) [default %default]")
)

# 解析命令行
opts = parse_args(OptionParser(option_list=option_list))

if (opts$output == "") { opts$output = opts$input }

# 安装和加载包
# install.packages('ggtern')
suppressWarnings(suppressMessages(library(ComplexHeatmap)))
suppressWarnings(suppressMessages(library(circlize)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(readr)))

# 读取数据
data01 <- read.table(opts$input, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

# 提取热图使用的数据列
data_matrix <- as.matrix(data01[, 3:20])

# z-score标准化
data_matrix = apply(data_matrix, 1, function(x){
  return((x-mean(x))/sd(x))
})

data_matrix = as.matrix(t(data_matrix))

# 设置颜色（等价于 pheatmap 的 colorRampPalette）
col_fun <- colorRamp2(c(min(data_matrix), 0, max(data_matrix)), c("#5e7fc2", "white", "#d05771"))

# 设置列注释（annotation_col）
group_info <- data.frame(
  Group = factor(rep(c("Feces", "Plaque", "Saliva"), c(6, 6, 6))),
  row.names = colnames(data_matrix)
)

# 设置颜色注释（列组颜色）
group_color <- list(Group = c("Feces" = "#d73027", "Plaque" = "#80cdc1", "Saliva" = "#BC80BD"))

# 创建列注释对象
column_ha <- HeatmapAnnotation(df = group_info, col = group_color)

# 创建热图对象
ht <- Heatmap(data_matrix,
              name = "Expression",  # 热图图例的名称
              col = col_fun,        # 颜色映射函数
              top_annotation = column_ha,  # 列注释
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              row_km = 3,           # 类似 cutree_rows = 3
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              border = FALSE)

# 保存为 PDF
pdf(paste0(opts$output, "heatmap_amplicon05.pdf"), width = 7.5, height = 7)  # Adjust width and height as needed
draw(ht)
dev.off()

