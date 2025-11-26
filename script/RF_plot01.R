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
  make_option(c("-i", "--input"), type="character", default="result2/tax/RF_model/Species_imp_rf21.txt",
              help="Model results table [default %default]"),
  make_option(c("-t", "--tax"), type = "character", default = "result2/tax/taxonomy_amplicon.txt", 
              help = "Taxonomy information [default %default]"),
  make_option(c("-n", "--optimal"), type = "numeric", default = "", 
              help = "Optimal biomarker number [default %default]"),
  make_option(c("-o", "--output"), type="character", default="result2/tax/RF_model",
              help="Output directory [default %default]"),
  make_option(c("-w", "--width"), type="numeric", default=120 * 1.5,
              help="Figure width (mm) [default %default]"),
  make_option(c("-e", "--height"), type="numeric", default=80 * 1.5,
              help="Figure height (mm) [default %default]")
)

# 解析命令行
opts = parse_args(OptionParser(option_list=option_list))

if (opts$output == "") { opts$output = opts$input }

# 基于CRAN安装R包，检测没有则安装
p_list = c("reshape2","ggplot2","ggprism","dplyr","plyr","cols4all",
           "openxlsx","tidyverse","readr")
for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
  library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

# 加载R包 Load the package
suppressWarnings(suppressMessages(library("reshape2")))
suppressWarnings(suppressMessages(library("ggplot2")))
suppressWarnings(suppressMessages(library("ggprism")))
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("plyr")))
suppressWarnings(suppressMessages(library("cols4all")))
suppressWarnings(suppressMessages(library("openxlsx")))
suppressWarnings(suppressMessages(library("tidyverse")))
suppressWarnings(suppressMessages(library("readr")))


# 读取 taxonomy 文件（假设是 tab 分隔）
#tax <- read_tsv("data/taxonomy_amplicon.txt", col_names = TRUE)
tax <- read_tsv(opts$tax, col_names = TRUE)
# X7 为匹配列，X3 为 Phylum
phylum_map <- tax %>% select(Phylum, Genus) %>% distinct()
# 读取目标表
#species <- read_tsv("results/model_amplicon/Species_imp_rf21.txt")
species <- read_tsv(opts$input)
# 假设第1列是物种名，对应 taxonomy 第7列
colnames(species)[1] <- "Genus"
# 左连接
imp_species <- left_join(species, phylum_map, by = c("Genus" = "Genus"))
# 重命名列
colnames(imp_species)[ncol(imp_species)] <- "Phylum"
# 写入结果
#write_tsv(merged, "results/model_amplicon/Species_imp_phylum_rf.txt")

optimal = opts$optimal

## 绘制模型挑选出的物种的柱状图，横轴为MeanDecreaseAccuracy，纵轴为物种，以Phylum着色
#imp_species = read.table("results/model_RF/Species_imp_phylum_rf2.txt", header=T, row.names= 1, sep="\t")
#imp_species = read.table("results/model_amplicon/Species_imp_rf2_amplicon.txt", header=T, row.names= 1, sep="\t")
imp_species = tail(imp_species, n = optimal)
imp_species$Species = factor(rownames(imp_species), levels = rownames(imp_species))
p03_species = ggplot(imp_species, aes(x = Species, y = MeanDecreaseAccuracy, fill = Phylum)) + 
  geom_bar(stat = "identity") + theme_classic()+
  #scale_fill_manual(values = c("#63B8FF","orange","#4AB3AA", "#D10640"))+
  #  scale_color_manual(values = c("#63B8FF", "orange","#4AB3AA","#D10640"))+
  scale_fill_manual(values = c("#63B8FF","#4AB3AA", "#dfc6a5", "#e5acd7","#D10640"))+
  scale_color_manual(values = c("#63B8FF","#4AB3AA", "#dfc6a5", "#e5acd7","#D10640"))+
  coord_flip() + #main_theme+
  theme(legend.position = c(0.85,0.8))+
  scale_y_continuous(expand = c(0,0))+
  labs(y = "Mean Decrease Accuracy", x = "Genus")
#ggsave(paste("results/model_RF/Species_top_feautre_top6markers_phylum",".pdf", sep=""), p03_species, width=119 * 1.5, height=80 * 1.5, unit='mm')
#ggsave(paste("results/model_amplicon/Species_top_feautre_top24markers_phylum2",".pdf", sep=""), p03_species, width=119 * 1.5, height=110 * 1.5, unit='mm')
#ggsave(paste0(opts$output, "Species_top_feautre_top24markers_phylum2.pdf"), plot = p03_species, width = 119 * 1.5, height = 110 * 1.5, dpi = 600)
ggsave(paste0(opts$output, "Species_top_feautre_topmarkers_phylum2.pdf"), plot = p03_species, width = 6, height = 10, dpi = 600)
p03_species
