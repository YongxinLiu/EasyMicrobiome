rm(list=ls())

# 安装和加载包
# install.packages('ggtern')
library(ggplot2)
library(ggtern)
library(reshape2)
library(tidyverse)

# 载入数据
# Load data

data <- read.table("result/tax/sum_g.txt",header = TRUE,row.names = 1,sep = "\t")
group <- read.delim("result/metadata.txt", sep = '\t', stringsAsFactors = FALSE)

# # 计算所有功能相对丰度并且以百分比的形式表示
# data = as.data.frame(t(t(data)/colSums(data,na=T)*100))

#筛选相对丰度大于1%的功能
data  =  data %>% filter(apply(data,1,mean) > 1)
data = t(data)# 行列转换
data = data[rownames(data) != "All", ]  # 删除行名为 "All" 的行
# 假设data和group是两个已经存在的数据框
data1 <- cbind(data, group)
data1$group = as.factor(data1$Group)

# 筛选出KC组和MC组的数据，并保留group列
data_subset <- data1 %>%
  filter(group %in% c("feces", "plaque")) %>%
  select(group, all_of(names(data1)[!names(data1) %in% "group"])) %>%
  select(group, everything())  # 选择group列和所有其他列

# 方差齐，进行t检验
dfg = data_subset %>%
  select_if(is.numeric) %>%
  map_df(~ broom::tidy(t.test(. ~ group,data = data_subset)), .id = 'var')
# p矫正“holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”
dfg$p.value = p.adjust(dfg$p.value,"none")# p未矫正
dfg = dfg %>% filter(p.value < 0.05)# 筛选p值显著的功能0.05

# 方差不齐wilcox检验
# dfg = data_subset %>%
#   select_if(is.numeric) %>%
#   map_df(~ broom::tidy(wilcox.test(. ~ group,data = data_subset)), .id = 'var')
# dfg$p.value <- p.adjust(dfg$p.value,"BH")
# dfg = dfg %>% filter(p.value < 0.05)

# 筛选差异显著的功能，将宽数据转化为长数据，计算平均值
abun.bar = data_subset[,c(dfg$var,"group")] %>% 
  gather(variable,value,-group) %>% #将宽数据转化为长数据
  group_by(variable,group) %>% 
  summarise(Mean = mean(value))
# 筛选dfg结果中的信息
diff.mean <- dfg[,c("var","estimate","conf.low","conf.high","p.value")]
diff.mean$group <- c(ifelse(diff.mean$estimate >0,levels(data1$group)[1],
                            levels(data1$group)[2]))
diff.mean <- diff.mean[order(diff.mean$estimate,decreasing = TRUE),]


# 先绘制左侧条形图
col <- c("#edae11","#0f6657")
# 按照功能的因子排序
abun.bar$variable  =  factor(abun.bar$variable,levels = rev(diff.mean$var))
p1 = ggplot(abun.bar,aes(variable,Mean,fill = group)) +
  scale_fill_manual(values=col)+
  scale_x_discrete(limits = levels(diff.mean$var)) +
  coord_flip() +
  xlab("") +
  ylab("Mean proportion (%)") +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size=12,face = "bold"),
        axis.text=element_text(colour='black',size=10,face = "bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=12,face = "bold",colour = "black",
                                 margin = margin(r = 20)),
        legend.position = "top",# 图例位置设置上方
        legend.direction = "horizontal",
        legend.key.width = unit(0.8,"cm"),
        legend.key.height = unit(0.5,"cm"))

p1
# 添加柱状图的背景
for (i in 1:(nrow(diff.mean) - 1)) 
  p1 <- p1 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf, 
                      fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
# 将柱状图添加
p1 = p1+geom_bar(stat = "identity",position = "dodge",
                 width = 0.7,colour = "black")
p1

# 右侧散点图绘制
diff.mean$var = factor(diff.mean$var,levels = levels(abun.bar$variable))
diff.mean$p.value = as.numeric(diff.mean$p.value)
diff.mean$p.value = round(diff.mean$p.value,3)# 保留3位小数
diff.mean$p.value = as.character(diff.mean$p.value)

p2 = ggplot(diff.mean,aes(var,estimate,fill = group)) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size=12,face = "bold"),
        axis.text=element_text(colour='black',size=10,face = "bold"),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 15,face = "bold",colour = "black",hjust = 0.5)) +
  scale_x_discrete(limits = levels(diff.mean$var)) +
  coord_flip() +
  xlab("") +
  ylab("Difference in mean proportions (%)") +
  labs(title="95% confidence intervals") 
p2
for (i in 1:(nrow(diff.mean) - 1)) 
  p2 <- p2 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf, 
                      fill = ifelse(i %% 2 == 0, 'white', 'gray95'))

p2 <- p2 +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                position = position_dodge(0.8), width = 0.5, linewidth = 0.5) +# 误差线
  geom_point(shape = 21,size = 3) +
  scale_fill_manual(values=col) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black')

p2 

# 最右侧p值文本添加   
p3 <- ggplot(diff.mean,aes(var,estimate,fill = group)) +
  geom_text(aes(y = 0,x = var),label = diff.mean$p.value,
            hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
  geom_text(aes(x = nrow(diff.mean)/2 +0.5,y = 0.85),label = "P-value (corrected)",
            srt = 90,fontface = "bold",size = 5) +
  coord_flip() +
  ylim(c(0,1)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
p3

## 画到这里有了三个图p1、p2、p3
## 我们只需要把它们拼接在一起就好啦
library(patchwork)
p <- p1 + p2 + p3 + plot_layout(widths = c(4,6,2))

p
## ggsave保存图
ggsave("stamp.pdf", plot = p, width = 9.5, height = 6, dpi = 600)
