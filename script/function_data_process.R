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
  make_option(c("-i", "--input"), type="character", default="result2/tax/KEGG.PathwayL2.raw.txt",
              help="Function annotation table [default %default]"),
  make_option(c("-g", "--group"), type = "character", default = "result2/tax/metadata_amplicon.txt", 
              help = "Group information [default %default]"),
  make_option(c("-o", "--output"), type="character", default="result2/tax/pathway",
              help="Output directory [default %default]"),
  make_option(c("-w", "--width"), type="numeric", default=120 * 1.5,
              help="Figure width (mm) [default %default]"),
  make_option(c("-e", "--height"), type="numeric", default=80 * 1.5,
              help="Figure height (mm) [default %default]")
)

# 解析命令行
opts = parse_args(OptionParser(option_list=option_list))

if (opts$output == "") { opts$output = opts$input }

# 基于CRAN安装R包，检测没有则安装 Installing R packages based on CRAN and installing them if they are not detected
p_list = c("ggplot2", "patchwork", "dplyr", "reshape2", "ggprism", "plyr",
           "magrittr","ggfun","cowplot","DESeq2","edgeR" )
for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
  library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

# 加载R包 Loading R packages
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(patchwork)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(ggprism)))
suppressWarnings(suppressMessages(library(plyr)))
suppressWarnings(suppressMessages(library(magrittr)))
suppressWarnings(suppressMessages(library(ggfun)))
suppressWarnings(suppressMessages(library(cowplot)))


# 读取您的数据
data <- read.table(opts$input, header=TRUE, sep="\t")

# 处理PathwayL2列，将空格替换为下划线
data$PathwayL2 <- gsub(" ", "_", data$PathwayL2)

# 移除KO列（如不需要）
data$KO <- NULL

# 写入新文件
write.table(data, paste0(opts$output, "pathway_count_data.txt"), sep="\t", row.names=FALSE, quote=FALSE)


# 生成分组文件
#design <- read.table("Figure3/data/metadata_amplicon.txt", header=TRUE, sep="\t")
design <- read.table(opts$group, header=TRUE, sep="\t")
simple_design <- data.frame(
  sample = design$SampleID,
  Group = design$Group,
  sample_new = design$SampleID
)

# 分别保存不同组的设计文件
write.table(subset(simple_design, Group=="feces"), paste0(opts$output, "design_feces.txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(subset(simple_design, Group=="plaque"), paste0(opts$output, "design_plaque.txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(subset(simple_design, Group=="saliva"), paste0(opts$output, "design_saliva.txt"), sep="\t", row.names=FALSE, quote=FALSE)

pathways <- data.frame(
  ID = data[[1]],  # 获取 data 第一列的所有值
  Pathway = gsub("_", " ", data[[1]])  # 假设第一列包含路径名称并去掉下划线
)
#write.table(pathways, "Figure3/data/Difference_pathway.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(pathways, paste0(opts$output, "Difference_pathway.txt"), sep="\t", row.names=FALSE, quote=FALSE)

# 1. 差异分析函数
perform_diff_analysis <- function(count_data, design_data, group1_name, group2_name) {
  # 创建DGEList对象
  dge <- DGEList(counts = count_data, 
                 group = design_data$Group)
  
  # 过滤和标准化（保持不变）
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  
  # 创建设计矩阵
  design <- model.matrix(~ 0 + design_data$Group)
  colnames(design) <- c(group1_name, group2_name)
  
  # 对比矩阵
  contrast <- makeContrasts(
    contrasts = paste(group2_name, group1_name, sep = "-"),
    levels = design
  )
  
  # 估计离散度和拟合模型
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  
  # 执行检验
  qlf <- glmQLFTest(fit, contrast = contrast)
  
  # 获取结果
  res <- topTags(qlf, n = Inf)$table
  
  # 添加比较组信息
  res$Pathway <- rownames(res)
  res$Group <- paste(group1_name, "vs", group2_name)
  
  # 重命名和选择列
  final_res <- res %>%
    mutate(B_coef = logFC) %>%
    select(Pathway, Group, FDR, B_coef)
  
  return(final_res)
}


# 2. 读取数据
#pathway_counts <- read.table("Figure3/data/pathway_count_data.txt", header=TRUE, sep="\t")
pathway_counts <- data
#design_feces <- read.table("Figure3/data/design_feces.txt", header=TRUE)
design_feces <- subset(simple_design, Group=="feces")
#design_plaque <- read.table("Figure3/data/design_plaque.txt", header=TRUE)
design_plaque <- subset(simple_design, Group=="plaque")
#design_saliva <- read.table("Figure3/data/design_saliva.txt", header=TRUE)
design_saliva <- subset(simple_design, Group=="saliva")

# 3. 执行差异分析（示例：比较feces vs plaque）
# 合并feces vs plaque两个组的数据
fp_counts <- cbind(
  pathway_counts[, 1, drop = FALSE], # 保留第一列
  pathway_counts[, c(design_feces$sample, design_plaque$sample)] # 选择feces和plaque样本的列
)
# #设置第一列为行名
rownames(fp_counts) <- fp_counts[, 1]  # 将第一列作为行名
fp_counts <- fp_counts[, -1]  # 删除第一列（因为已经作为行名）

# #确保数据是数值型矩阵，适合差异分析
fp_counts <- as.matrix(fp_counts)

fp_design <- rbind(
  data.frame(sample = design_feces$sample, Group = "feces"),
  data.frame(sample = design_plaque$sample, Group = "plaque")
)

# #运行差异分析
# #执行分析时指定比较组名称
diff_fp <- perform_diff_analysis(
  fp_counts,
  fp_design,
  group1_name = "feces",
  group2_name = "plaque"
)

# 合并feces vs saliva两个组的数据
fs_counts <- cbind(
  pathway_counts[, 1, drop = FALSE], # 保留第一列
  pathway_counts[, c(design_feces$sample, design_saliva$sample)] # 选择feces和plaque样本的列
)
# #设置第一列为行名
rownames(fs_counts) <- fs_counts[, 1]  # 将第一列作为行名
fs_counts <- fs_counts[, -1]  # 删除第一列（因为已经作为行名）

# #确保数据是数值型矩阵，适合差异分析
fs_counts <- as.matrix(fs_counts)

fs_design <- rbind(
  data.frame(sample = design_feces$sample, Group = "feces"),
  data.frame(sample = design_saliva$sample, Group = "saliva")
)

# #运行差异分析
# #执行分析时指定比较组名称
diff_fs <- perform_diff_analysis(
  fs_counts,
  fs_design,
  group1_name = "feces",
  group2_name = "saliva"
)

# 合并plaque vs saliva两个组的数据
ps_counts <- cbind(
  pathway_counts[, 1, drop = FALSE], # 保留第一列
  pathway_counts[, c(design_plaque$sample, design_saliva$sample)] # 选择feces和plaque样本的列
)
# #设置第一列为行名
rownames(ps_counts) <- ps_counts[, 1]  # 将第一列作为行名
ps_counts <- ps_counts[, -1]  # 删除第一列（因为已经作为行名）

# #确保数据是数值型矩阵，适合差异分析
ps_counts <- as.matrix(ps_counts)

ps_design <- rbind(
  data.frame(sample = design_plaque$sample, Group = "plaque"),
  data.frame(sample = design_saliva$sample, Group = "saliva")
)

# #运行差异分析
# #执行分析时指定比较组名称
diff_ps <- perform_diff_analysis(
  ps_counts,
  ps_design,
  group1_name = "plaque",
  group2_name = "saliva"
)

# 4. 创建最终结果文件
all_results <- bind_rows(diff_fp, diff_fs, diff_ps) 
difference_pathway <- data.frame(
  ID = 1:nrow(all_results),
  Pathway = all_results$Pathway,
  Group = all_results$Group,
  FDR = all_results$FDR,
  B_coef = all_results$B_coef
)

# 5. 保存结果
#write.table(difference_pathway, "Figure3/data/Difference_pathway2.txt", 
#            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(difference_pathway, paste0(opts$output, "Difference_pathway2.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)


# 计算流行率
# Load data
#data_pathway <- read.table("Figure3/data/pathway_count_data.txt", sep = "\t", header = TRUE, check.names = FALSE)
data_pathway <- data
#design_A <- read.table("Figure3/data/design_feces.txt", sep = "\t", header = TRUE, row.names = 1)
design_A <- design_feces
#design_B <- read.table("Figure3/data/design_plaque.txt", sep = "\t", header = TRUE, row.names = 1)
design_B <- design_plaque
#design_C <- read.table("Figure3/data/design_saliva.txt", sep = "\t", header = TRUE, row.names = 1)
design_C <- design_saliva
#difference_pathway <- read.table("Figure3/data/Difference_pathway.txt", sep = "\t", header = TRUE, row.names = 1)
difference_pathway <- pathways

# Preprocess data
rownames(data_pathway) <- data_pathway$PathwayL2
data_pathway <- data_pathway[, -1] %>% apply(2, function(x) x / sum(x))
data_pathway02 <- as.data.frame(t(data_pathway))

conflicted::conflicts_prefer(dplyr::mutate)
# Function to process each design group
process_design <- function(design, data_pathway02, all_counts_value) {
  data_pathway_group <- data_pathway02[rownames(data_pathway02) %in% rownames(design), ] %>% t() %>% as.data.frame()
  
  zero_counts <- rowSums(data_pathway_group == 0)
  data_Pathway2 <- data_pathway_group %>%
    mutate(zero_counts = zero_counts,
           sample_percent = round(1 - zero_counts / all_counts_value, 6))
  
  data_species3 <- data_Pathway2[rownames(data_Pathway2) %in% rownames(difference_pathway), ]
  return(data_species3)
}

# Process each design group and write to CSV
data_species3_C <- process_design(design_A, data_pathway02, 30)
#write.csv(data_species3_C, "Figure3/data/data_pathway3_prevalence_A.csv")
write.csv(data_species3_C, paste0(opts$output, "data_pathway3_prevalence_feces.csv"))

data_species3_E <- process_design(design_B, data_pathway02, 35)
#write.csv(data_species3_E, "Figure3/data/data_pathway3_prevalence_B.csv")
write.csv(data_species3_E, paste0(opts$output, "data_pathway3_prevalence_plaque.csv"))

data_species3_Y <- process_design(design_C, data_pathway02, 50)
#write.csv(data_species3_Y, "Figure3/data/data_pathway3_prevalence_C.csv")
write.csv(data_species3_Y, paste0(opts$output, "data_pathway3_prevalence_saliva.csv"))

