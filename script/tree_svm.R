rm(list=ls())
if (!requireNamespace("ggtreeExtra", quietly = TRUE)) {
  remotes::install_github("YuLab-SMU/ggtreeExtra")
}
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(ggnewscale)
library(phangorn)
library(RColorBrewer)
# install.packages("ape")  # 如果你还没有安装ape包
library(ape)
##R语言构建进化树
# # 读入OTU数据、tax表
# OTU <- read.table("otutab.txt", sep = "\t",  row.names = 1,stringsAsFactors =FALSE, check.names =FALSE,header=1)
# tax <- read.table("taxonomy.txt", sep = "\t",stringsAsFactors =FALSE, check.names =FALSE,header=1)
# 
# # 计算每一行的和并添加到数据框
# OTU$sum <- rowSums(OTU)
# 
# # 按每行的和降序排列
# OTU1 <- OTU[order(OTU$sum, decreasing = TRUE), ]
# 
# # 删除sum列
# OTU1 <- OTU1[, -which(names(OTU1) == "sum")]
# 
# # 取出丰度排名前100的OTU，减少展示
# OTU <- OTU1[1:100, ]
# 
# #使用bray curtis方法计算距离矩阵
# otu_dist <- vegdist(OTU,method = 'bray')
# 
# head(otu_dist)
# 
# # 使用类平均法进行聚类
# otu_hc1 <- hclust(otu_dist,method="average")
# 
# # 使用默认函数绘制进化树
# plot(as.dendrogram(otu_hc1),type="rectangle",horiz=T)
# 
# # 将聚类结果转成系统发育格式
# otu_tree <- as.phylo(otu_hc1)
# 
# #将数据导出为newick格式文件
# write.tree(phy=otu_tree, file="tree.nwk")


#读取树文件（本树文件由fasttree构建）
tree <- read.tree('D:/EasyAmplicon_paper_materials/PacBio/result/tree/otus.nwk')
#在没有确定祖先节点的情况下，用进化距离中点给树定根，使结果更对称、更有生物学解释力
tree <- midpoint(tree)


#热图数据
df <- read.delim("D:/EasyAmplicon_paper_materials/PacBio/result/tree/annotation.txt", header = TRUE, check.names = FALSE)
df_phylum <- df %>% select(OTUID, Phylum)

# 给树加分组：按 Phylum 把 tip 分组
tree <- groupOTU(tree, split(df_phylum$OTUID, df_phylum$Phylum))
# 使用 ggtree 绘图，并根据分组上色树枝
p1 <- ggtree(tree, layout = "fan", open.angle = 30, branch.length = 'none',
             aes(color = group), size = 0.3) +
  scale_color_manual(values = colorset2) +     # 设置颜色
  geom_tiplab(size = 2, align = TRUE, linetype = "dotted") +  # 添加标签
  theme(legend.title = element_blank())        # 去掉图例标题
p1


# 添加第二圈（按 Sample 上色）
sample_colors <- c(
  feces = "#edae11",
  plaque = "#0f6657",
  saliva = "#c74732"
)

p2 <- p1 + 
  geom_fruit(
    data = df_long,
    geom = geom_tile,
    mapping = aes(y = OTUID, x = Sample, fill = Sample, alpha = Abundance),
    offset = 0.15, size = 0.1,
    axis.params = list(axis = "x", text.angle = -90, text.size = 2)
  ) +
  scale_fill_manual(values = sample_colors) + 
  scale_alpha(range = c(0.1, 6)) +
  theme(legend.position = "right")

p2

## ggsave保存图
ggsave("tree.pdf", plot = p2, width = 9.5, height = 6, dpi = 600)
