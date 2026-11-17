#!/usr/bin/env Rscript

# Copyright 2024-2026 Defeng Bai <baidefeng@caas.cn>

# If used this script, please cited:
# Bai, et al. 2025. EasyMetagenome: A User‐Friendly and Flexible Pipeline for Shotgun Metagenomic Analysis in Microbiome Research. iMeta 4: e70001. https://doi.org/10.1002/imt2.70001

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

# 更新
# 2024/11/12：增加利用MaAsLin2差异分析结果绘制火山图代码
# 2025/11/27：规范代码
# 2026/5/22：更新依赖包安装
# 2026/11/17: 更新输入结果为MaAsLin3差异分析结果

# Copyright 2024-2026 Defeng Bai <baidefeng@caas.cn>
# 火山图展示组间差异（适配 MaAsLin3 输出格式）
options(warn = -1)

## ========================= 1. 依赖包安装 =========================
site <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN"

install_if_missing <- function(pkg, type = "CRAN", repo = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (type == "CRAN") {
      install.packages(pkg, repos = site)
    } else if (type == "BIOC") {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = site)
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else if (type == "GITHUB") {
      if (!requireNamespace("devtools", quietly = TRUE))
        install.packages("devtools", repos = site)
      devtools::install_github(repo, force = TRUE)
    }
  }
}

cran_pkgs <- c("ggplot2", "BiocManager", "optparse", "dplyr", "magrittr", "ggrepel")
for (pkg in cran_pkgs) install_if_missing(pkg, "CRAN")

bioc_pkgs <- c("multcompView")
for (pkg in bioc_pkgs) install_if_missing(pkg, "BIOC")

install_if_missing("amplicon",  "GITHUB", repo = "microbiota/amplicon")
install_if_missing("ggvolcano", "GITHUB", repo = "BioSenior/ggvolcano")

## ========================= 2. 命令行参数 =========================
if (!suppressWarnings(suppressMessages(
  require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos = site)
  require("optparse", character.only = TRUE)
}

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              default = "result/metaphlan4/MaAsLin3_result/MaAsLin3_overall_difference.csv",
              help = "MaAsLin3 overall difference result [default %default]"),
  make_option(c("-d", "--trend"), type = "character",
              default = "result/metaphlan4/MaAsLin3_result/MaAsLin3_enriched_depleted.csv",
              help = "MaAsLin3 enriched/depleted trend table [default %default]"),
  make_option(c("-o", "--output"), type = "character",
              default = "result/metaphlan4/",
              help = "Output directory [default %default]"),
  make_option(c("-w", "--width"), type = "numeric", default = 89,
              help = "Figure width (mm) [default %default]"),
  make_option(c("-e", "--height"), type = "numeric", default = 59,
              help = "Figure height (mm) [default %default]")
)
opts <- parse_args(OptionParser(option_list = option_list))
print("You are using the following parameters:")
print(opts)

## ========================= 3. 加载依赖 =========================
package_list <- c("ggplot2", "dplyr", "magrittr", "ggrepel", "ggvolcano")
for (p in package_list) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = site)
  suppressWarnings(suppressMessages(library(p, character.only = TRUE)))
}

## ========================= 4. 读取输入（修复重复行名报错） =========================

# 趋势表：MaAsLin3_enriched_depleted.csv
# 不用 row.names = 1，避免重复行名报错
data_trend <- read.csv(opts$trend, header = TRUE, check.names = FALSE)

# 如有多组比较，可选定一组（按实际 Group_value 取值修改，默认不过滤）
# data_trend <- data_trend[data_trend$Group_value == "Centenarians", ]

data_trend$species <- data_trend$Feature
data_trend_use <- data_trend[, c("species", "Association")]

# 差异汇总表：MaAsLin3_overall_difference.csv
data_MWAS <- read.csv(opts$input, header = TRUE, check.names = FALSE)

# 如多模型，可选定 abundance（prevalence 的 Mean/FC 常为 NA，默认不过滤）
# data_MWAS <- data_MWAS[data_MWAS$Model == "abundance", ]

data_MWAS$species <- data_MWAS$Feature

# merge（按 species）
data_vol <- merge(data_trend_use, data_MWAS, by = "species")
data_vol <- as.data.frame(data_vol)

# 计算 -log10(FDR)
data_vol$padj2 <- -log10(data_vol$FDR)

# 剔除 NA（prevalence 模型的 FDR/Beta 可能为 NA）
data_vol <- data_vol[!is.na(data_vol$padj2) & !is.na(data_vol$Beta), ]

# 确认 Association 列存在
if (!"Association" %in% colnames(data_vol)) {
  stop("merge 后缺少 Association 列，请检查趋势表列名是否为 'Association'")
}
print(table(data_vol$Association, useNA = "ifany"))

## ========================= 5. 火山图 =========================
assoc_levels <- unique(data_vol$Association)
color_map <- c(
  "enriched"        = "#a60026",
  "depleted"        = "#74add1",
  "not_significant" = "grey"
)
color_map <- color_map[names(color_map) %in% assoc_levels]
missing_lv <- setdiff(assoc_levels, names(color_map))
if (length(missing_lv) > 0)
  color_map <- c(color_map, setNames(rep("grey", length(missing_lv)), missing_lv))

p_volcano1 <- ggplot(data = data_vol, aes(x = -Beta, y = padj2)) +
  geom_point(alpha = 0.4, size = 2.0, aes(color = Association)) +
  scale_color_manual(values = color_map) +
  geom_vline(xintercept = 0, lty = 4, col = "black", lwd = 0.4) +
  geom_hline(yintercept = -log10(0.75), lty = 4, col = "black", lwd = 0.4) +
  labs(x = "Coef. (by MaAsLin3)",
       y = bquote(atop(-Log[10]~italic(FDR)))) +
  theme_bw()

## ========================= 6. 标注 Top 5 =========================
up_data1 <- dplyr::filter(data_vol, Association == "enriched")
up_data2 <- dplyr::arrange(up_data1, dplyr::desc(padj2))
up_data_5 <- head(up_data2, 5)

down_data1 <- dplyr::filter(data_vol, Association == "depleted")
down_data2 <- dplyr::arrange(down_data1, dplyr::desc(padj2))
down_data_5 <- head(down_data2, 5)

p_volcano2 <- p_volcano1 +
  geom_text_repel(data = up_data_5,
                  aes(x = -Beta, y = padj2, label = Feature),
                  size = 2, max.overlaps = 20) +
  geom_text_repel(data = down_data_5,
                  aes(x = -Beta, y = padj2, label = Feature),
                  size = 2, max.overlaps = 20) +
  theme(
    legend.position = c(0.20, 0.75),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    panel.grid = element_blank()
  )

## ========================= 7. 保存 =========================
if (!dir.exists(opts$output)) dir.create(opts$output, recursive = TRUE)
ggsave(paste0(opts$output, "Difference_valcano.pdf"),
       p_volcano2, width = opts$width, height = opts$height, unit = "mm")
cat("Done. Output:", paste0(opts$output, "Difference_valcano.pdf"), "\n")
