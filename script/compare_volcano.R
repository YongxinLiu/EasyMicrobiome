#!/usr/bin/env Rscript

# Copyright 2024-2026 Defeng Bai <baidefeng@caas.cn>
# If used this script, please cited:
# Bai, et al. 2025. EasyMetagenome: A User‐Friendly and Flexible Pipeline for
# Shotgun Metagenomic Analysis in Microbiome Research. iMeta 4: e70001.
# https://doi.org/10.1002/imt2.70001
#
# 更新
# 2024/11/12：增加利用 MaAsLin2 差异分析结果绘制火山图代码
# 2025/11/27：规范代码
# 2026/5/22：更新依赖包安装
# 2026/11/17：整合 MaAsLin2/MaAsLin3 两种输出格式，自动检测列名

# 程序功能：火山图展示组间差异（兼容 MaAsLin2 / MaAsLin3）
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
              help = "MaAsLin2/3 overall difference result [default %default]"),
  make_option(c("-d", "--trend"), type = "character",
              default = "result/metaphlan4/MaAsLin3_result/MaAsLin3_enriched_depleted.csv",
              help = "MaAsLin2/3 enriched/depleted trend table [default %default]"),
  make_option(c("-o", "--output"), type = "character",
              default = "result/metaphlan4/",
              help = "Output directory [default %default]"),
  make_option(c("-w", "--width"), type = "numeric", default = 89,
              help = "Figure width (mm) [default %default]"),
  make_option(c("-e", "--height"), type = "numeric", default = 59,
              help = "Figure height (mm) [default %default]"),
  make_option(c("-m", "--model"), type = "character", default = NULL,
              help = "For MaAsLin3: keep only this Model (e.g. abundance); default keep all non-NA"),
  make_option(c("-g", "--group_value"), type = "character", default = NULL,
              help = "For MaAsLin3: keep only this Group_value (e.g. Centenarians); default keep all")
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

## ========================= 4. 读取输入（自动兼容 MaAsLin2/3） =========================

# 工具函数：去掉行名，仅保留数据列
read_no_rownames <- function(file) {
  df <- read.csv(file, header = TRUE, check.names = FALSE)
  # 若第一列是无名列/行名（如空字符串），去掉
  if (ncol(df) > 0 && (colnames(df)[1] == "" || is.na(colnames(df)[1])))
    df <- df[, -1, drop = FALSE]
  df
}

# 工具函数：找列名（忽略大小写、空格、点、下划线）
find_col <- function(df, patterns) {
  nm <- colnames(df)
  norm <- gsub("[[:space:]_.-]", "", tolower(nm))
  for (p in patterns) {
    p_norm <- gsub("[[:space:]_.-]", "", tolower(p))
    idx <- which(norm == p_norm)
    if (length(idx) > 0) return(nm[idx[1]])
  }
  # 再做一次模糊匹配
  for (p in patterns) {
    p_norm <- gsub("[[:space:]_.-]", "", tolower(p))
    idx <- grep(p_norm, norm, fixed = TRUE)
    if (length(idx) > 0) return(nm[idx[1]])
  }
  return(NA_character_)
}

# ---- 读取趋势表 ----
data_trend <- read_no_rownames(opts$trend)

# 找 Feature 列
trend_feature <- find_col(data_trend, c("Feature", "feature"))
if (is.na(trend_feature)) stop("趋势表找不到 Feature 列")
# 找 Association 列（MaAsLin2: "NPC association" / "NPC.association"；MaAsLin3: "Association"）
trend_assoc <- find_col(data_trend, c("Association", "NPC association",
                                      "NPC.association", "NPCassociation"))
if (is.na(trend_assoc)) stop("趋势表找不到 Association / NPC association 列")
cat("趋势表：Feature =", trend_feature, "；Association =", trend_assoc, "\n")

# 若 MaAsLin3 且用户指定 Group_value，先过滤
if (!is.null(opts$group_value) && "Group_value" %in% colnames(data_trend)) {
  data_trend <- data_trend[data_trend$Group_value == opts$group_value, ]
  cat("趋势表按 Group_value =", opts$group_value, "过滤后行数：", nrow(data_trend), "\n")
}

data_trend_use <- data.frame(
  species     = data_trend[[trend_feature]],
  Association = data_trend[[trend_assoc]],
  stringsAsFactors = FALSE
)
# 统一 Association 取值（MaAsLin3 可能是 not_significant）
data_trend_use$Association <- as.character(data_trend_use$Association)
data_trend_use$Association[data_trend_use$Association %in% c("not_significant", "not significant", "ns")] <- "not_significant"

# 若趋势表同一物种多行（多组比较），按物种去重，保留第一条
data_trend_use <- data_trend_use[!duplicated(data_trend_use$species), ]

# ---- 读取 overall 表 ----
data_MWAS <- read_no_rownames(opts$input)

mw_feature <- find_col(data_MWAS, c("Feature", "feature"))
if (is.na(mw_feature)) stop("overall 表找不到 Feature 列")
mw_beta <- find_col(data_MWAS, c("Beta", "coef", "Coefficient"))
mw_fdr  <- find_col(data_MWAS, c("FDR", "qvalue", "q_value", "padj", "adj.P.Val"))
mw_p    <- find_col(data_MWAS, c("P", "pvalue", "p_value", "P.Value"))
if (is.na(mw_beta) || is.na(mw_fdr)) stop("overall 表找不到 Beta 或 FDR 列")
cat("overall 表：Feature =", mw_feature, "；Beta =", mw_beta, "；FDR =", mw_fdr, "\n")

# MaAsLin3：若指定 Model，先过滤
if (!is.null(opts$model) && "Model" %in% colnames(data_MWAS)) {
  data_MWAS <- data_MWAS[data_MWAS$Model == opts$model, ]
  cat("overall 表按 Model =", opts$model, "过滤后行数：", nrow(data_MWAS), "\n")
}
# MaAsLin3：若指定 Group_value，先过滤
if (!is.null(opts$group_value) && "Group_value" %in% colnames(data_MWAS)) {
  data_MWAS <- data_MWAS[data_MWAS$Group_value == opts$group_value, ]
  cat("overall 表按 Group_value =", opts$group_value, "过滤后行数：", nrow(data_MWAS), "\n")
}

data_MWAS_use <- data.frame(
  species = data_MWAS[[mw_feature]],
  Beta    = suppressWarnings(as.numeric(data_MWAS[[mw_beta]])),
  FDR     = suppressWarnings(as.numeric(data_MWAS[[mw_fdr]])),
  stringsAsFactors = FALSE
)
# 若同一物种多行（多模型/多组），保留 FDR 最小的一行
data_MWAS_use <- data_MWAS_use[!is.na(data_MWAS_use$Beta) & !is.na(data_MWAS_use$FDR), ]
if (any(duplicated(data_MWAS_use$species))) {
  data_MWAS_use <- data_MWAS_use[order(data_MWAS_use$species, data_MWAS_use$FDR), ]
  data_MWAS_use <- data_MWAS_use[!duplicated(data_MWAS_use$species), ]
}

# ---- merge ----
data_vol <- merge(data_trend_use, data_MWAS_use, by = "species")
data_vol <- as.data.frame(data_vol)
data_vol$padj2 <- -log10(data_vol$FDR)
data_vol <- data_vol[!is.na(data_vol$padj2) & !is.na(data_vol$Beta), ]

cat("合并后参与绘图的物种数：", nrow(data_vol), "\n")
print(table(data_vol$Association, useNA = "ifany"))

if (nrow(data_vol) == 0) {
  stop("合并后没有可用数据。请检查：\n",
       "  - 趋势表与 overall 表的 Feature 名称是否一致；\n",
       "  - MaAsLin3 是否需用 -m abundance / -g <Group_value> 指定模型与分组；\n",
       "  - FDR/Beta 是否全为 NA。")
}

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
  labs(x = "Coef. (by MaAsLin)", y = bquote(atop(-Log[10]~italic(FDR)))) +
  theme_bw()

## ========================= 6. 标注 Top 5 =========================
top_n <- 5
up_data1 <- dplyr::filter(data_vol, Association == "enriched")
up_data2 <- dplyr::arrange(up_data1, dplyr::desc(padj2))
up_data_5 <- head(up_data2, top_n)

down_data1 <- dplyr::filter(data_vol, Association == "depleted")
down_data2 <- dplyr::arrange(down_data1, dplyr::desc(padj2))
down_data_5 <- head(down_data2, top_n)

p_volcano2 <- p_volcano1
if (nrow(up_data_5) > 0) {
  p_volcano2 <- p_volcano2 +
    geom_text_repel(data = up_data_5,
                    aes(x = -Beta, y = padj2, label = species),
                    size = 2, max.overlaps = 20)
}
if (nrow(down_data_5) > 0) {
  p_volcano2 <- p_volcano2 +
    geom_text_repel(data = down_data_5,
                    aes(x = -Beta, y = padj2, label = species),
                    size = 2, max.overlaps = 20)
}
p_volcano2 <- p_volcano2 +
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
out_file <- paste0(opts$output, "Difference_valcano.pdf")
ggsave(out_file, p_volcano2, width = opts$width, height = opts$height, unit = "mm")
cat("Done. Output:", out_file, "\n")