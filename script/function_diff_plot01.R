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
  make_option(c("-i", "--input"), type="character", default="result2/tax/pathway/Difference_pathway21.txt",
              help="Data table for plots [default %default]"),
  make_option(c("-p", "--pathway"), type = "character", default = "result2/tax/pathway/pathway_count_data.txt", 
              help = "Group information [default %default]"),
  make_option(c("-s", "--statistic"), type = "character", default = "result2/tax/pathway/pathway_percent_abundance2.txt", 
              help = "Function relative abundance and prevalence information [default %default]"),
  make_option(c("-o", "--output"), type="character", default="result2/tax/pathway/",
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
conflicted::conflicts_prefer(dplyr::mutate)

# 载入数据
# Load data
data <- read.table(opts$input,header = TRUE,row.names = 1,sep = "\t")
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
  scale_fill_gradient2(low='#0085c1', high='red',mid = 'white', limit=c(-1.6,1.55),name=paste0("B-coef.")) +
  labs(x=NULL,y=NULL) +  
  theme_classic()+
  theme(axis.text.x = element_text(size=10,angle = 45,hjust = 1,color = "black"),            
        axis.text.y = element_text(size=10,color = "black"), 
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.length = unit(2.0, "mm"),
        panel.background=element_blank(),
        legend.position = "left")

# 1. 读取原始数据（确保字符串不作为因子读取）
pathway_data <- read.table(opts$pathway,
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

conflicted::conflicts_prefer(dplyr::summarise)
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
write.table(result, paste0(opts$output, "pathway_percent_abundance.txt"),
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
# 查看结果（显示前6行）
head(result)

# Bar plot prevalence
# 柱状图展示每个物种的流行率
imp_species <- read.table(opts$statistic, header = TRUE, sep = "\t")
imp_species <- imp_species %>%
  mutate(Pathway = ordered(Pathway, levels = pathway_levels))

# 获取所有唯一的组名
groups <- unique(imp_species$group)
# 为每个组指定不同的颜色
colors <- c("feces" = "#edae11", "plaque" = "#0f6657", "saliva" = "#c74732")

# 使用 lapply 绘制每个组的图
plots <- lapply(seq_along(groups), function(i) {
  g <- groups[i]
  p <- ggplot(subset(imp_species, group == g), aes(x = Pathway, y = Abundance)) +
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

p2 <- ggdraw() +
  draw_plot(p1, 0, -0.0395, 0.5, 0.983)+
  draw_plot(final_plot, 0.48, 0, 0.5, 0.97)


#pdf("Figure3/results/age_fungi_heatmap_ra_bar2.pdf", height = 10.2, width = 14)
pdf(paste0(opts$output, "amplicon_function_diff_heatmap01.pdf"), width = 14, height = 10.2)  # Adjust width and height as needed
p2
dev.off()


