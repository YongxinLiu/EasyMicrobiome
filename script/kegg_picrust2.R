# 读取您的数据
data <- read.table("KEGG.PathwayL2.raw.txt", header=TRUE, sep="\t")

# 处理PathwayL2列，将空格替换为下划线
data$PathwayL2 <- gsub(" ", "_", data$PathwayL2)

# 移除KO列（如不需要）
data$KO <- NULL

# 写入新文件
write.table(data, "pathway_count_data.txt", sep="\t", row.names=FALSE, quote=FALSE)


design <- read.table("../result/metadata.txt", header=TRUE, sep="\t")
simple_design <- data.frame(
  sample = design$SampleID,
  Group = design$Group,
  sample_new = design$SampleID
)

# 分别保存不同组的设计文件
write.table(subset(simple_design, Group=="feces"), "design_feces.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(subset(simple_design, Group=="plaque"), "design_plaque.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(subset(simple_design, Group=="saliva"), "design_saliva.txt", sep="\t", row.names=FALSE, quote=FALSE)


pathways <- data.frame(
  ID = data[[1]],  # 获取 data 第一列的所有值
  Pathway = gsub("_", " ", data[[1]])  # 假设第一列包含路径名称并去掉下划线
)
write.table(pathways, "Difference_pathway.txt", sep="\t", row.names=FALSE, quote=FALSE)


# 安装差异分析包（如果尚未安装）
if (!requireNamespace("DESeq2")) install.packages("DESeq2")
if (!requireNamespace("edgeR")) install.packages("edgeR")
library(DESeq2)
library(edgeR)


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
pathway_counts <- read.table("pathway_count_data.txt", header=TRUE, sep="\t")
design_feces <- read.table("design_feces.txt", header=TRUE)
design_plaque <- read.table("design_plaque.txt", header=TRUE)
design_saliva <- read.table("design_saliva.txt", header=TRUE)

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
write.table(difference_pathway, "Difference_pathway2.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)


# 基于CRAN安装R包，检测没有则安装 Installing R packages based on CRAN and installing them if they are not detected
p_list = c("ggplot2", "patchwork", "dplyr", "reshape2", "ggprism", "plyr",
           "magrittr","ggfun","cowplot" )
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


# 计算流行率
# Load data
data_pathway <- read.table("pathway_count_data.txt", sep = "\t", header = TRUE, check.names = FALSE)
design_A <- read.table("design_feces.txt", sep = "\t", header = TRUE, row.names = 1)
design_B <- read.table("design_plaque.txt", sep = "\t", header = TRUE, row.names = 1)
design_C <- read.table("design_saliva.txt", sep = "\t", header = TRUE, row.names = 1)
difference_pathway <- read.table("Difference_pathway.txt", sep = "\t", header = TRUE, row.names = 1)

# Preprocess data
rownames(data_pathway) <- data_pathway$PathwayL2
data_pathway <- data_pathway[, -1] %>% apply(2, function(x) x / sum(x))
data_pathway02 <- as.data.frame(t(data_pathway))

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
write.csv(data_species3_C, "data_pathway3_prevalence_A.csv")

data_species3_E <- process_design(design_B, data_pathway02, 35)
write.csv(data_species3_E, "data_pathway3_prevalence_B.csv")

data_species3_Y <- process_design(design_C, data_pathway02, 50)
write.csv(data_species3_Y, "data_pathway3_prevalence_C.csv")


# 载入数据
# Load data
data <- read.table("Difference_pathway2.txt",header = TRUE,row.names = 1,sep = "\t")
data[which(data$FDR<0.05),'sig'] <- '*'

# 排序
# Set order
# 获取 Pathway 列中的所有唯一值作为排序的 levels
pathway_levels <- unique(data$Pathway)

# 将 Pathway 列的值按照 Pathway 列中的唯一值进行排序
data$Pathway2 = factor(data$Pathway, levels = pathway_levels)

# 如果需要将 Pathway2 转换为有序因子
data = data %>%
  mutate(Pathway2 = ordered(Pathway2, levels = pathway_levels))

# 绘图
# Plot
p1 <- ggplot(data, aes(Group,Pathway2)) +  
  geom_tile(aes(fill=B_coef), color="white") +    
  geom_text(aes(label=sig), color="black", size=6,vjust = 0.8) + 
  scale_fill_gradient2(low='#0085c1', high='red',mid = 'white', limit=c(-1.3,1.3),name=paste0("B-coef.")) +
  labs(x=NULL,y=NULL) +  
  theme_classic()+
  theme(axis.text.x = element_text(size=10,angle = 45,hjust = 1,color = "black"),            
        axis.text.y = element_text(size=10,color = "black"), 
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.length = unit(2.0, "mm"),
        panel.background=element_blank(),
        legend.position = "left")
#p1
#ggsave(paste("results/age_fheatmap01",".pdf", sep=""), p1, width=89 * 1.5, height=180 * 1.5, unit='mm')


# 1. 加载必要的包
library(dplyr)

# 1. 读取原始数据（确保字符串不作为因子读取）
pathway_data <- read.table("pathway_count_data.txt", 
                           header = TRUE,
                           sep = "\t",
                           row.names = 1,
                           check.names = FALSE,
                           stringsAsFactors = FALSE)

# 2. 定义分组
group_info <- data.frame(
  sample = colnames(pathway_data),
  group = rep(c("feces", "plaque", "saliva"), each = 6),
  stringsAsFactors = FALSE
)

# 3. 计算组内百分比（不筛选0值）
result <- group_info %>%
  group_by(group) %>%
  group_modify(~ {
    samples <- .x$sample
    group_data <- pathway_data[, samples, drop = FALSE]
    
    # 计算每个pathway在组内的总和
    pathway_sums <- rowSums(group_data)
    
    # 计算全组总和
    group_total <- sum(pathway_sums)
    
    # 计算百分比（乘以100，不四舍五入）
    data.frame(
      Pathway = rownames(group_data),
      `Relative abundance` = 100 * pathway_sums / group_total,  # 直接乘以100
      check.names = FALSE
    )
  }) %>%
  ungroup()

# 4. 验证每组总和是否为100%
result %>%
  group_by(group) %>%
  summarise(total = sum(`Relative abundance`)) %>%
  print()

# 5. 保存结果（保留所有pathway，包括0值）
write.table(result, "pathway_percent_abundance.txt",
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# 查看结果（显示前6行）
head(result)

#p1
#ggsave(paste("results/age_fheatmap01",".pdf", sep=""), p1, width=89 * 1.5, height=180 * 1.5, unit='mm')

# Bar plot prevalence
# 柱状图展示每个物种的流行率
imp_species <- read.table("pathway_percent_abundance.txt", header = TRUE, sep = "\t")
imp_species <- imp_species %>%
  mutate(Pathway = ordered(Pathway, levels = pathway_levels))

# 获取所有唯一的组名
groups <- unique(imp_species$Group)

# 为每个组指定不同的颜色
colors <- c("feces" = "#edae11", "plaque" = "#0f6657", "saliva" = "#c74732")


# 使用 lapply 绘制每个组的图
plots <- lapply(seq_along(groups), function(i) {
  g <- groups[i]
  p <- ggplot(subset(imp_species, Group == g), aes(x = Pathway, y = Abundance)) +
    geom_bar(stat = "identity", fill = colors[g]) +
    coord_flip() +
    theme_classic() +
    theme(text = element_text(family = 'sans', size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_blank()
    ) +
    geom_hline(yintercept = c(5, 10, 15, 20), colour = 'black', lwd = 0.36, linetype = "dashed") +
    facet_set(label = paste(g)) +
    labs(y = "RA(%)")
  
  if (i != 1) {
    p <- p + theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank())
  }
  
  return(p)
})


# 使用 plot_grid 将三个图放在一起，并确保宽度一致
final_plot <- plot_grid(plotlist = plots, ncol = 3, rel_widths = c(1, 1, 1))

# 显示最终图形
#print(final_plot)

library(cowplot)
p2 <- ggdraw() +
  draw_plot(p1, 0, 0.0015, 0.5, 0.933)+
  draw_plot(final_plot, 0.48, 0, 0.5, 0.97)


#p2

pdf("age_fungi_heatmap_ra_bar2.pdf", height = 7.2, width = 8)
p2
dev.off()