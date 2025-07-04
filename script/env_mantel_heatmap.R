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
  make_option(c("-i", "--input"), type="character", default="result2/tax/otutab2.txt",
              help="OTU composition [default %default]"),
  make_option(c("-n", "--env"), type = "character", default = "result2/tax/env_amplicon.txt", 
              help = "Environment variables"),
  make_option(c("-o", "--output"), type="character", default="result2/tax/",
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
# 基于CRAN安装R包，检测没有则安装
p_list = c("ggplot2", "dplyr", "vegan")
for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
  library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

# 基于github安装
library(devtools)
if(!requireNamespace("ggradar", quietly = TRUE))
  install_github("Hy4m/linkET", force = TRUE)

# install.packages('ggtern')
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(vegan)))
suppressWarnings(suppressMessages(library(linkET)))

# Load and transform data
## OTU data
df <- read.table(opts$input, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
df <- data.frame(t(df))

## Environmental data
env <- read.table(opts$env, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
env <- as.data.frame(env)

# Perform Mantel test and create custom labels
df_mantel <- mantel_test(df, env,
                         spec_select = list(
                           "East" = 1:7,
                           "North" = 8:9,
                           "Outside" = 10:11,
                           "South" = 12
                         )) %>%
  mutate(
    df_r = cut(r, breaks = c(-Inf, 0.25, 0.5, Inf),
               labels = c("< 0.25", "0.25 - 0.5", ">= 0.5")),
    df_p = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
               labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))
  )

# 连线在右上(The connection is in the top right)
# Save the plot
# pdf("results/mantel_plot_amplicon.pdf", width = 10, height = 8)  # Adjust width and height as needed
pdf(paste0(opts$output, "mantel_test_right.pdf"), width = 10, height = 8)  # Adjust width and height as needed
# Plotting
qcorrplot(correlate(env), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(
    aes(
      xend = .xend + 1.25,  # Define end of line for connections
      yend = .yend + 0.5,
      colour = df_p,
      size = df_r
    ),
    data = df_mantel,
    curvature = 0.1
  ) +
  geom_diag_label(
    mapping = aes(y = .y + 0.05),  # Define label positions on the diagonal
    hjust = 0.15
  ) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +  # Color scale for squares
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(
    size = guide_legend(
      title = "Mantel's r",
      override.aes = list(colour = "grey35"),
      order = 2
    ),
    colour = guide_legend(
      title = "Mantel's p",
      override.aes = list(size = 3),
      order = 1
    ),
    fill = guide_colorbar(title = "Pearson's r", order = 3)
  ) +
  theme(
    axis.text.y = element_blank()  # Remove y-axis text
  )
dev.off()  # Close the PDF device


# 连线在左下(The connection is in the bottom left)
# Save the plot
#pdf("results/mantel_plot2_amplicon2.pdf", width = 12, height = 8)  # Adjust width and height as needed
pdf(paste0(opts$output, "mantel_test_left.pdf"), width = 10, height = 8)  # Adjust width and height as needed
# Plotting
qcorrplot(correlate(env), type = "upper", diag = FALSE) +
  geom_square() +
  geom_couple(
    aes(
      xend = .xend - 0.5,  # Define end of line for connections
      yend = .yend - 0.25,
      colour = df_p,
      size = df_r
    ),
    data = df_mantel,
    curvature = 0.1
  ) +
  geom_diag_label(
    mapping = aes(y = .y + 0.05),  # Define label positions on the diagonal
    hjust = 0.25
  ) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +  # Color scale for squares
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(
    size = guide_legend(
      title = "Mantel's r",
      override.aes = list(colour = "grey35"),
      order = 2
    ),
    colour = guide_legend(
      title = "Mantel's p",
      override.aes = list(size = 3),
      order = 1
    ),
    fill = guide_colorbar(title = "Pearson's r", order = 3)
  ) +
  theme(
    axis.text.y = element_blank()  # Remove y-axis text
  )
dev.off()  # Close the PDF device

