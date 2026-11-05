#!/usr/bin/env Rscript

# 加载必要的库
suppressPackageStartupMessages({
  library(optparse)
  library(amplicon)
  library(ggtern)
})

# 定义命令行选项
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="输入文件路径（OTU表）", metavar="FILE"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="元数据文件路径", metavar="FILE"),
  make_option(c("-g", "--group"), type="character", default="Group",
              help="元数据中的分组列名 [默认 %default]"),
  make_option(c("-o", "--output"), type="character", default="Ternary.pdf",
              help="输出PDF文件路径 [默认 %default]"),
  make_option(c("-t", "--topn"), type="integer", default=10,
              help="显示的前N个分类单元，其余归为Others [默认 %default]"),
  make_option(c("-W", "--width"), type="numeric", default=9,
              help="图形宽度（英寸） [默认 %default]"),
  make_option(c("-H", "--height"), type="numeric", default=9,
              help="图形高度（英寸） [默认 %default]"),
  make_option(c("-p", "--palette"), type="character", default="Set3",
              help="颜色方案（RColorBrewer调色板） [默认 %default]"),
  make_option(c("-l", "--taxlevel"), type="character", default="Phylum",
              help="分类水平（如Phylum/Class/Order/Family/Genus） [默认 %default]")
)

# 解析参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查必要参数
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("必须提供输入文件路径", call.=FALSE)
}

# 验证分类水平是否有效
valid_taxlevels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
if (!opt$taxlevel %in% valid_taxlevels) {
  stop(paste("无效的分类水平，请选择以下之一:", paste(valid_taxlevels, collapse=", ")), call.=FALSE)
}

# 读取数据
tax_table <- read.table(opt$input, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# 提取分类单元名称（根据用户指定的分类水平修改列名）
colnames(tax_table)[1] <- opt$taxlevel
tax_id <- tax_table[, 1]

# 确保数值列为 numeric 类型
tax_numeric <- tax_table[, -1]
tax_numeric[] <- lapply(tax_numeric, function(x) as.numeric(as.character(x)))

# 计算相对丰度
rela_tax_values <- apply(tax_numeric, 2, function(x) x / sum(x, na.rm=TRUE))

# 重新组合数据
rela_tax <- data.frame(Taxa=tax_id, rela_tax_values, stringsAsFactors=FALSE)
colnames(rela_tax)[1] <- opt$taxlevel
colnames(rela_tax)[-1] <- colnames(tax_table)[-1]

# 如果有元数据文件，按分组计算均值
if (!is.null(opt$metadata)) {
  metadata <- read.table(opt$metadata, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  if (!opt$group %in% colnames(metadata)) {
    stop(paste("分组列", opt$group, "不在元数据文件中"), call.=FALSE)
  }
  
  # 获取各组的样本列
  groups <- unique(metadata[[opt$group]])
  group_cols <- lapply(groups, function(g) {
    which(colnames(rela_tax) %in% metadata$SampleID[metadata[[opt$group]] == g])
  })
  
  # 计算各组的平均相对丰度
  rela_group_tax <- t(apply(rela_tax[, -1], 1, function(x) {
    sapply(group_cols, function(cols) mean(x[cols], na.rm=TRUE))
  }))
  
  colnames(rela_group_tax) <- groups
} else {
  # 如果没有元数据，假设前中后各6个样本为三组
  warning("未提供元数据，将自动按列顺序分组（每组6个样本）", call.=FALSE)
  group_cols <- list(2:7, 8:13, 14:19)
  groups <- c("Group1", "Group2", "Group3")
  rela_group_tax <- t(apply(rela_tax[, -1], 1, function(x) {
    sapply(group_cols, function(cols) mean(x[cols], na.rm=TRUE))
  }))
  colnames(rela_group_tax) <- groups
}

# 转换为数据框
df <- data.frame(Taxa=tax_id, rela_group_tax, stringsAsFactors=FALSE)
colnames(df)[1] <- opt$taxlevel

# 计算总丰度并排序
df$total <- rowSums(df[, -1, drop=FALSE])
sort_df <- df[order(df$total, decreasing=TRUE), ]

# 取前N个分类单元，其余归为Others
top_n <- sort_df[1:min(opt$topn, nrow(sort_df)), ]

# 定义 start_row 和 end_row
start_row <- min(opt$topn, nrow(sort_df)) + 1
end_row <- nrow(sort_df)

# 计算 Others 的丰度（避免索引越界）
others <- if (start_row <= end_row) {
  colSums(sort_df[start_row:end_row, -1, drop = FALSE])
} else {
  setNames(rep(0, ncol(sort_df) - 1), colnames(sort_df)[-1])
}
# 合并数据
df2 <- rbind(top_n[, -1, drop = FALSE], others)

# 计算平均丰度（排除最后一列）
df2$Abundance <- rowMeans(df2[, 1:(ncol(df2) - 1)], na.rm = TRUE)

# 添加分类信息
df2[[opt$taxlevel]] <- factor(
  c(top_n[[opt$taxlevel]], "Others"),
  levels = c(top_n[[opt$taxlevel]], "Others")
)

# 绘制三元图
p <- ggtern(data=df2, aes_string(x=colnames(df2)[1], y=colnames(df2)[2], z=colnames(df2)[3])) +
  geom_point(aes_string(color=opt$taxlevel, size="Abundance")) +
  theme_bw() +
  theme_arrowdefault() +
  theme(
    legend.text=element_text(size=12),
    legend.title=element_text(size=10, face="bold"),
    legend.key=element_rect(fill="transparent"),
    legend.key.size=unit(0.5, "cm"),
    axis.title=element_text(size=10, face="bold"),
    axis.text=element_text(size=10)
  ) +
  guides(color=guide_legend(override.aes=list(size=3))) +
  scale_size_continuous(range=c(4, 10)) +
  scale_color_brewer(palette=opt$palette) +
  theme(legend.position="right") +
  labs(title=paste("Ternary Plot at", opt$taxlevel, "Level"))

# 保存结果
ggsave(opt$output, plot=p, device="pdf", width=opt$width, height=opt$height)

cat(paste0(opt$taxlevel, "水平的三元图已保存到: ", normalizePath(opt$output), "\n"))