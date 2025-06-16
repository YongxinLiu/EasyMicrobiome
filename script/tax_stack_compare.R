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
  make_option(c("-i", "--input"), type="character", default="result2/tax/data_practice53_16S.txt",
              help="Taxonomy composition [default %default]"),
  make_option(c("-c", "--compare"), type = "character", default = "Illumina-PacBio", 
              help = "Compared group"),
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
# install.packages('ggtern')
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggprism)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(dplyr)))

compare <- opts$compare

# 2.将宽数据转换为长数据绘图
df <- read.table(opts$input, sep = "\t", header = TRUE, check.names = FALSE)

# 宽数据转换为长数据
df_long <- df %>%
  pivot_longer(cols = c("feces1", "feces2", "feces3", "feces4", "feces5", 
                        "feces6"),  # 指定要转换的所有列
               names_to = "Samples",  # 新列名
               values_to = "values") #%>%  # 存储数据的列名
  #mutate(values2 = (values / sum(values)) * 100)  # 计算 values2 列，百分比

# 创建新的数据框
#df_long_final <- df_long %>%
#  select(Genus = Genus, Samples, values, values2, Group)
df_long_final <- df_long %>%
  select(Genus = Genus, Samples, values, Group)

# 可以将数据保存为新的文件
# write.table(df_long_final, file = "result2/tax/converted_data.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# 绘制堆叠柱状图并增加facet标签
p1 <- ggplot(df_long_final, aes(x = Samples, y = values, fill = Genus)) +
  geom_col(width = 0.75, position = position_stack(vjust = 1)) +  # 堆叠柱状图
  scale_fill_manual(values = c("#d2da93","#5196d5","#00ceff","#ff630d","#e50719",
                               "#e5acd7","#77aecd","#ec8181","#dfc6a5","#35978b",
                               "#d27e43","#a6bddb","#8DD3C7","#fe5094","#fdbb84","#f94e54",
                               "#CAB2D6", "#FFFF99",  "#FFED6F", "#BC80BD",
                               "#ffad00","#36999d","#00fc8d","#b64aa0","#9b82e1",
                               "#66c2a4","#8c96c6","#8d342e","#8a4984","#fa9fb5",
                               "#66c2a4")) +  # 自定义颜色
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +  # y轴扩展和限制
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.001)) +  # y轴扩展和限制
  scale_x_discrete(expand = c(0.1, 0.1)) +  # x轴扩展
  coord_flip() +  # 横向堆叠柱状图
  labs(x = "", y = "Relative abundance", title = "") +  # 设置标签
  theme_prism(base_fontface = "plain", base_family = "serif", base_size = 16, base_line_size = 0.8, axis_text_angle = 0) +  # 设置主题
  theme_bw() +
  theme(axis.title.y = element_text(face = "bold", size = 14, color = "black", family = "sans")) +  # y轴标题样式
  facet_wrap(~Group, scales = "fixed", ncol = 2, labeller = labeller(group = label_value)) +  # 保持y轴相同
  theme(strip.text = element_text(face = "bold", size = 16, color = "black"),  # 控制facet标签的样式
        legend.position = "bottom") # 将图例放置在下方

# 输出图形
p1
ggsave(paste0(opts$output, compare, "_", "stacked_composition.pdf"), plot = p1, width = 10, height = 8, dpi = 600)

