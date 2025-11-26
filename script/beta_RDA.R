# 清理工作环境
rm(list = ls())

# 检查并安装必要的R包
packages_needed <- c("microeco", "ape", "magrittr", "ggplot2")
packages_missing <- packages_needed[!(packages_needed %in% installed.packages()[,"Package"])]
if(length(packages_missing)) {
  install.packages(packages_missing)
}

# microeco一般需要从GitHub安装
if(!require("microeco")){
  if(!require("devtools")) install.packages("devtools")
  devtools::install_github("ChiLiubio/microeco")
}
# 加载包
library(microeco)
library(ape)
library(magrittr)
library(ggplot2)

# 设定主题
theme_set(theme_bw())

# 设置随机种子，保证可复现
set.seed(123)

# 读入数据
sample_info_16S <- read.delim("Nnapore_env/result/metadata.txt", row.names = 1, header = TRUE, sep = "\t")
otu_table_16S <- read.delim("Nnapore_env/result/otutab.txt", row.names = 1, header = TRUE, sep = "\t")
taxonomy_table_16S <- read.delim("Nnapore_env/result/taxonomy.txt", row.names = 1, header = TRUE, sep = "\t")
phylo_tree_16S <- read.tree("Nnapore_env/result/otus.tree")
# env_data_16S <- read.delim("Nnapore_env/env_data.txt", row.names = 1, header = TRUE, sep = "\t")

# 处理分类表：统一前缀
taxonomy_table_16S %<>% tidy_taxonomy

# 创建microtable对象
mt <- microtable$new(sample_table = sample_info_16S,
                     otu_table = otu_table_16S,
                     tax_table = taxonomy_table_16S,
                     phylo_tree = phylo_tree_16S)
print(mt)

# 删除非细菌和古菌的OTU
mt$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
print(mt)

# 修剪数据集
mt$tidy_dataset()
print(mt)

# 过滤线粒体、叶绿体污染
mt$filter_pollution(taxa = c("mitochondria", "chloroplast"))
mt$tidy_dataset()
print(mt)

# 查看样本测序深度范围
print(mt$sample_sums() %>% range)

# 样本稀释标准化
mt$rarefy_samples(sample.size = 98)
print(mt)

# 保存处理后的基础数据表
mt$save_table(dirpath = "basic_files", sep = ",")

# 将环境因子加入样本表,因为metadata文件里有环境因子故省略
# mt$sample_table <- data.frame(mt$sample_table, env_data_16S[rownames(mt$sample_table), ])
# print(mt)

# 创建环境因子分析对象
t1 <- trans_env$new(dataset = mt, env_cols = 4:7)

# dbRDA分析，使用bray-curtis距离
t1$cal_ordination(method = "dbRDA", use_measure = "bray")

# 可视化（第一次）
t1$trans_ordination()
t1$plot_ordination(plot_color = "Group")

# 调整箭头长度，更美观
t1$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1.5)

# 重新绘制并保存图片
p11 <- t1$plot_ordination(plot_color = "Group")
print(p11)
ggsave("p11.pdf", p12, width = 6, height = 4)


# 保存R环境
save(mt, t1, file = "microeco_analysis_result.RData")


# 在 Genus 水平做 RDA 分析
t1$cal_ordination(method = "RDA", taxa_level = "Genus")

# 选择贡献度最高的10个分类单元（OTUs），并且调整箭头长度
t1$trans_ordination(
  show_taxa = 10,
  adjust_arrow_length = TRUE,
  max_perc_env = 1.5,
  max_perc_tax = 1.5,
  min_perc_env = 0.2,
  min_perc_tax = 0.2
)

# t1$res_rda_trans 是转换后的结果，可以直接用于画图
p12 <- t1$plot_ordination(plot_color = "Group")
p12

# 保存图片
ggsave("p12.pdf", p12, width = 6, height = 4)



# First create trans_diff object as a demonstration
# Using random forest method to find differentially abundant genera between groups
t2 <- trans_diff$new(
  dataset = mt,
  method = "rf",          # Random forest method
  group = "Group",        # Grouping variable
  taxa_level = "Genus"    # Analysis at genus level
)

# Then create trans_env object for environmental correlation analysis
# Adding environmental data (columns 4-11 from env_data_16S)
t1 <- trans_env$new(
  dataset = mt,
  add_data = env_data_16S[, 4:11]
)

# Calculate correlations using significant taxa from trans_diff results
# Selecting top 40 differentially abundant genera
t1$cal_cor(
  use_data = "other",             # Using custom taxa list
  p_adjust_method = "fdr",        # FDR p-value adjustment
  other_taxa = t2$res_diff$Taxa[1:40]  # Top 40 significant genera
)

# Generate correlation heatmap plot
p14 <- t1$plot_cor()

# Display the plot
p14

# Save high-resolution image
ggsave(
  "p14.png",
  p14,
  width = 8,      # 8 inches wide
  height = 12,    # 12 inches tall
  dpi = 1200      # High resolution
)