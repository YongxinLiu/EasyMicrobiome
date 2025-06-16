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
  make_option(c("-i", "--input"), type="character", default="result2/tax/PacBio_data.txt",
              help="OTU/ASV table [default %default]"),
  make_option(c("-g", "--group"), type = "character", default = "result2/tax/PacBio_metadata.txt", 
              help = "Group information [default %default]"),
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
p_list = c("reshape2","ggplot2","ggprism","dplyr","plyr","caret",
           "randomForest","PRROC","pROC","yardstick","patchwork","cols4all",
           "openxlsx","tidyverse")
for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
  library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

# 加载R包 Load the package
suppressWarnings(suppressMessages(library("reshape2")))
suppressWarnings(suppressMessages(library("ggplot2")))
suppressWarnings(suppressMessages(library("ggprism")))
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("plyr")))
suppressWarnings(suppressMessages(library("caret")))
suppressWarnings(suppressMessages(library("randomForest")))
suppressWarnings(suppressMessages(library("PRROC")))
suppressWarnings(suppressMessages(library("ROCR")))
suppressWarnings(suppressMessages(library("pROC")))
suppressWarnings(suppressMessages(library("yardstick")))
suppressWarnings(suppressMessages(library("patchwork")))
suppressWarnings(suppressMessages(library("cols4all")))
suppressWarnings(suppressMessages(library("openxlsx")))
suppressWarnings(suppressMessages(library("tidyverse")))

# 载入设置和函数, 这里主要用到了里面的main_theme绘图
#source("scripts/stat_plot_functions.R")
#source("scripts/randomforest.crossvalidation.R")
# 交叉验证的函数
#For windows system#
##ramdomforest.crossvalidation.r##
##Begin##
rfcv1 <-
  function (trainx, trainy, cv.fold = 5, scale = "log", step = 0.5,
            mtry = function(p) max(1, floor(sqrt(p))), recursive = FALSE,
            ...)
  {
    classRF <- is.factor(trainy)
    n <- nrow(trainx)
    p <- ncol(trainx)
    if (scale == "log") {
      k <- floor(log(p, base = 1/step))
      n.var <- round(p * step^(0:(k - 1)))
      same <- diff(n.var) == 0
      if (any(same))
        n.var <- n.var[-which(same)]
      if (!1 %in% n.var)
        n.var <- c(n.var, 1)
    }
    else {
      n.var <- seq(from = p, to = 1, by = step)
    }
    k <- length(n.var)
    cv.pred <- vector(k, mode = "list")
    for (i in 1:k) cv.pred[[i]] <- rep(0,length(trainy))
    if (classRF) {
      f <- trainy
    }
    else {
      f <- factor(rep(1:5, length = length(trainy))[order(order(trainy))])
    }
    nlvl <- table(f)
    idx <- numeric(n)
    for (i in 1:length(nlvl)) {
      idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold,
                                                  length = nlvl[i]))
    }
    res=list()
    for (i in 1:cv.fold) {
      all.rf <- randomForest(trainx[idx != i, , drop = FALSE],
                             trainy[idx != i],importance = TRUE)
      aa = predict(all.rf,trainx[idx == i, , drop = FALSE],type="prob")
      cv.pred[[1]][idx == i] <- as.numeric(aa[,2])
      impvar <- (1:p)[order(all.rf$importance[, 3], decreasing = TRUE)]
      res[[i]]=impvar
      for (j in 2:k) {
        12
        imp.idx <- impvar[1:n.var[j]]
        sub.rf <- randomForest(trainx[idx != i, imp.idx,
                                      drop = FALSE], trainy[idx != i]
        )
        bb <- predict(sub.rf,trainx[idx ==i,imp.idx, drop = FALSE],type="prob")
        cv.pred[[j]][idx == i] <- as.numeric(bb[,2])
        if (recursive) {
          impvar <- (1:length(imp.idx))[order(sub.rf$importance[,
                                                                3], decreasing = TRUE)]
        }
        NULL
      }
      NULL
    }
    if (classRF) {
      error.cv <- sapply(cv.pred, function(x) mean(factor(ifelse(x>0.5,1,0))!=trainy))
    }
    else {
      error.cv <- sapply(cv.pred, function(x) mean((trainy - x)^2))
    }
    names(error.cv) <- names(cv.pred) <- n.var
    list(n.var = n.var, error.cv = error.cv, predicted = cv.pred,res=res)
  }
##End##


mytheme = theme_classic() + 
  theme(text = element_text(family = "sans", size = 10))+
  theme(#legend.position="none",
    legend.text = element_text(size=8),
    legend.title = element_blank(), 
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(size=10, colour="black", family = "sans", angle = 0), 
    axis.text.x = element_text(size=10, colour="black", family = "sans", angle = 0, hjust = 0),
    axis.title= element_text(size=12, family = "sans"),
    strip.text.x = element_text(size=10, angle = 0),
    strip.text.y = element_text(size=10, angle = 0),
    panel.border = element_rect(colour = "black"),
    plot.title = element_text(size=10, angle = 0),
    strip.background.x = element_rect(fill = "#E5E4E2", colour = "black", size = 0.5),
    legend.position = c(0.85, 0.65),
  )+
  theme(axis.text.x=element_text(angle=0,vjust=1, hjust=0.6))+
  theme(axis.line = element_line(size = 0.2, colour = "black"))


# Diagnostic model using all filtered species features
# metadata 
#design <- read.table(file = "data/PacBio_metadata.txt", sep = "\t", header = T, row.names=1)
design <- read.table(opts$group, sep = "\t", header = T, row.names=1)

# 60 samples
#df_species <- read.table(file = "data/PacBio_data.txt", sep = "\t", header = T, check.names = FALSE)
df_species <- read.table(opts$input, sep = "\t", header = T, check.names = FALSE)

# sum of Species
# 计算每个Species微生物相对丰度之和，避免有重复Species统计
# data_species <- aggregate(.~ Species, data = df_species, sum)
# rownames(data_species) = data_species$Species
# data_species = data_species[, -1]
# data_species_ra = apply(data_species, 2, function(x) x/100)

df_species <- df_species[, c(5, 7:136)]

data_species <- aggregate(.~ Genus, data = df_species, sum)
rownames(data_species) = data_species$Genus
data_species = data_species[, -1]
#data_species_ra = apply(data_species, 2, function(x) x/100)

#1.微生物物种prevalence > 5%
# 创建一个空向量用于存放每一行的0个数数据
zero_counts <- vector("integer", nrow(data_species))
# 循环遍历每一行数据
for (i in 1:nrow(data_species)) {
  # 初始化计数器
  count <- 0
  # 循环遍历当前行数据的每个元素
  for (j in 1:ncol(data_species)) {
    # 判断当前元素是否为 0
    if (data_species[i, j] == 0) {
      count <- count + 1  # 计数器加一
    }
  }
  # 将当前行的0个数存放到结果向量中
  zero_counts[i] <- count
}
# 输出结果向量
zero_count = as.data.frame(zero_counts)
data_species2 = data_species
data_species2$zero_counts = zero_count$zero_counts
data_species2$all_counts = 130
data_species2$sample_percent = round(1-data_species2$zero_counts/data_species2$all_counts, 6)
data_species3 = data_species2 %>% filter(data_species2$sample_percent >= 0.05)
data_species3 = data_species3[, -c(131, 132, 133)]

#2.在样品占比大于5%的菌中，看是否在每个样品中对应的细菌丰度是否都超过0.01%，选取相对丰度超过0.01%的菌
data_species3 = apply(data_species3, 2, function(x) x/sum(x))
data_species3 = as.data.frame(data_species3)
count_t_values = apply(data_species3, 1, function(x)sum(x>=0.0001))
count_t_values = as.data.frame(count_t_values)
data_species3$count_t_values = count_t_values$count_t_values
data_species3$all_counts = 130
data_species3$t_percent = round(data_species3$count_t_values/data_species3$all_counts, 6)
data_species4 = data_species3 %>% filter(data_species3$t_percent >= 0.05)
data_species4 = data_species4[, -c(131, 132, 133)]

# 数据先进行log10转换，然后z-score标准化用于后续分析
data_species5 = log10(data_species4 + 1e-05)

# Data split
data_species6 = apply(data_species5, 1, function(x){
  return((x-mean(x))/sd(x))
})
data_species6 = t(data_species6)
#write.csv(data6, "results/rf_model_species_used.csv")

# 选取前面经过去重后的数据进行分析
otutab = data_species6
design2 = design

# Select by manual set group
if (TRUE){
  sub_design = subset(design2, Group %in% c("Male","Female")) 
  sub_design$group  = factor(sub_design$Group, levels=c("Male","Female"))
}
idx = rownames(sub_design) %in% colnames(otutab)
sub_design = sub_design[idx,]
sub_otutab = otutab[,rownames(sub_design)]

# Create data partition
# 将数据划分为训练集和测试集，# 这里大概按照7：3的比例划分训练集和测试集，60例样本的70%大约为42个，剩余18个样本约占30%
otutab_t_species = as.data.frame(t(sub_otutab))
# Set classification info.
otutab_t_species$group = factor(sub_design$Group, levels = c("Male","Female"))
otutab_t_species = na.omit(otutab_t_species)
row.name = rownames(otutab_t_species)
# 60 samples
set.seed = 515
sam.row.name = sample(row.name, 91, replace = FALSE)
train_data_species = otutab_t_species[sam.row.name, ]
unique_rows_df1 <- setdiff(rownames(otutab_t_species), rownames(train_data_species))
test_data_species <- otutab_t_species[unique_rows_df1, ]
#test_data_species = setdiff(otutab_t_species, train_data_species)

# Model training
# load data
dat1_species <- train_data_species
conf_species <- as.data.frame(dat1_species$group)
rownames(conf_species) <- rownames(dat1_species)
colnames(conf_species) <- "Group"
conf_species$sample <- rownames(conf_species)
conf_species <- as.data.frame(conf_species)
dat2_species <- dat1_species
conf2_species <- conf_species
conf2_species$Group = as.factor(as.character(conf2_species$Group))
outcome_species = conf2_species$Group
outcome_species <- sub("Male","0",outcome_species)
outcome_species <- sub("Female","1",outcome_species)
outcome_species <-as.factor(outcome_species)
dat_species <- dat2_species
X_species <- as.data.frame(dat_species)
X_species$outcome_species = outcome_species
X_species <- X_species[, -500]

######5*10_crossvalidation####
set.seed(999)
result_species <- replicate(5, rfcv1(X_species[,-ncol(X_species)], X_species$outcome_species, cv.fold=10,step=0.9), simplify=FALSE)
error.cv <- sapply(result_species, "[[", "error.cv")
matplot(result_species[[1]]$n.var, cbind(rowMeans(error.cv), error.cv), type="l",
        lwd=c(2, rep(1, ncol(error.cv))), col=1, lty=1, log="x",
        xlab="Number of variables", ylab="CV Error")
error.cv.cbm <- cbind(rowMeans(error.cv), error.cv)
cutoff <- min (error.cv.cbm[,1])+sd(error.cv.cbm[,1])
error.cv.cbm[error.cv.cbm[,1] < cutoff,]
#abline(v=24,col="pink",lwd=2)

error.cv.table <- error.cv.cbm[error.cv.cbm[,1] < cutoff,]
optimal = min(as.numeric(rownames(error.cv.table)))

#optimal = 24
error.cv.cbm2 <- as.data.frame(error.cv.cbm)
error.cv.cbm2$num <- rownames(error.cv.cbm2)
n.var = error.cv.cbm2$num
n.var = as.numeric(n.var)
error.cv = error.cv.cbm2[,1:5]
colnames(error.cv) = paste('err',1:5,sep='.')
err.mean = apply(error.cv,1,mean)
allerr = data.frame(num=n.var,err.mean=err.mean,error.cv)
allerr = as.data.frame(allerr)
# 对横坐标进行排序
#write.table(allerr, file = "results/model_RF/Species_rfcv_5_10_new.txt", sep = "\t", quote = F, row.names = T, col.names = T)
#write.table(allerr, file = "results/model_amplicon/Species_rfcv_5_10_new.txt", sep = "\t", quote = F, row.names = T, col.names = T)

#allerr <- read.table(file = "results/model_RF/Species_rfcv_5_10_new.txt", sep = "\t", header = T, row.names=1)
#allerr <- read.table(file = "results/model_amplicon/Species_rfcv_5_10_new.txt", sep = "\t", header = T, row.names=1)
mytheme3 = theme_bw() + theme(text = element_text(family = "sans", size = 7))+
  theme(legend.position="none",
        legend.text = element_text(size=14),
        legend.title = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=14, colour="black", family = "sans", angle = 0), 
        axis.text.x = element_text(size=14, colour="black", family = "sans", angle = 0, hjust = 0),
        axis.title= element_text(size=14),
        strip.text.x = element_text(size=14, angle = 0),
        strip.text.y = element_text(size=14, angle = 0),
        plot.title = element_text(size=14, angle = 0),
        strip.background.x = element_rect(fill = "#E5E4E2", colour = "black", size = 0.2))+
  theme(axis.text.x=element_text(angle=0,vjust=1, hjust=0.6))+
  theme(axis.line = element_line(size = 0.1, colour = "black"))

p01_species = ggplot(allerr, aes(x=allerr$num)) + 
  geom_line(data = allerr, aes(x = allerr$num, y = allerr$err.1), colour = 'grey') +
  geom_line(data = allerr, aes(x = allerr$num, y = allerr$err.2), colour = 'grey') +
  geom_line(data = allerr, aes(x = allerr$num, y = allerr$err.3), colour = 'grey') +
  geom_line(data = allerr, aes(x = allerr$num, y = allerr$err.4), colour = 'grey') +
  geom_line(data = allerr, aes(x = allerr$num, y = allerr$err.5), colour = 'grey') +
  geom_line(data = allerr, aes(x = allerr$num, y = allerr$err.mean), colour = 'black') + 
  geom_vline(xintercept = optimal, colour='black', lwd=0.36, linetype="dashed") + 
  #geom_hline(yintercept = 0.05976941, colour='black', lwd=0.36, linetype="dashed") +
  geom_hline(yintercept = cutoff, colour='black', lwd=0.36, linetype="dashed") +
  mytheme3+
  coord_trans(x = "log2") +
  scale_x_continuous(breaks = c(10, 30, 50, 100, 200, 400)) + # , max(allerr$num)
  labs(#title=paste('Training set (n = ', dim(train_data_species)[1],')', sep = ''),
    x='Number of species ', y='Cross-validation error rate') +
  annotate("text", x = optimal, y = max(allerr$err.mean), label=paste("optimal = ", optimal, sep="")) +
  #main_theme+ 
  theme_bw() + theme(panel.background = element_blank(),
                     panel.grid.major =element_blank(),
                     panel.grid.minor = element_blank(),
                     legend.position = "none",
                     axis.title= element_text(size=10, family = "sans"))
#ggsave("results/model_RF/Species_rfcv_5_10_top6.pdf",p01_species,width = 5,height = 3.2)
#ggsave("results/model_amplicon/Species_rfcv_5_10_top24.pdf",p01_species,width = 5,height = 3.2)
#ggsave("results/model_amplicon/Species_rfcv_5_10_top24.pdf",p01_species,width = 5,height = 3.2)
ggsave(paste0(opts$output, "Species_rfcv01.pdf"), plot = p01_species, width = 5, height = 3.2, dpi = 600)
#p01_species


#####pick 32 marker by corossvalidation#######
k=1
b <- matrix(0,ncol=499,nrow=50)
for(i in 1:5){
  for(j in 1:10){
    b[k,]<-result_species[[i]]$res[[j]]
    k=k+1
  }}
mlg.list<-b[,1:optimal]
list<-c()
k=1
for(i in 1:optimal){
  for(j in 1:50){
    list[k]<-mlg.list[j,i]
    k=k+1
  }}
mlg.sort<-as.matrix(table(list))
mlg.sort<-mlg.sort[rev(order(mlg.sort[,1])),]
pick_species<- as.numeric(names(head(mlg.sort,optimal)))
tmp= X_species[,-ncol( X_species)]
mlg.pick.species<-colnames(tmp)[pick_species]
#write.table(mlg.pick.species,"results/model_RF/cross_validation_pick_6_in_species.txt",
#            sep="\t",quote=F)
#write.table(mlg.pick.species,"results/model_amplicon/cross_validation_pick_6_in_species.txt",
#            sep="\t",quote=F)
write.table(mlg.pick.species,paste0(opts$output, "cross_validation_pick_6_in_species.txt"),
            sep="\t",quote=F)

## train.set
## 对比疾病组和健康对照组预测为疾病的概率
#train1_species <- X_species[,c(pick_species,612)]
train1_species <- X_species[,c(pick_species,500)]
train1_species <-data.frame(train1_species)
set.seed(32)
train1.rf_species <- randomForest(outcome_species~., data =train1_species,
                                  importance = TRUE)
train1.pre_species <- predict(train1.rf_species,type="prob")
p.train_species <- train1.pre_species[,2]
#boxplot(p.train~outcome,col=c(3,4),main="Probability of Patients")
#write.table(p.train_species,"results/model_RF/species.cross_validation.6makr.predict.in.train.txt",
#            sep="\t",quote=F)
#write.table(p.train_species,"results/model_amplicon/species.cross_validation.24makr.predict.in.train.txt",
#            sep="\t",quote=F)
write.table(p.train_species,paste0(opts$output, "species.cross_validation.24makr.predict.in.train.txt"),
            sep="\t",quote=F)

train1_pre2_species <- data.frame(outcome_species, p.train_species)
train1_pre2_species$outcome_species <- as.factor(train1_pre2_species$outcome_species)
train1_pre2_species$outcome_species <- sub("0","Male",train1_pre2_species$outcome_species)
train1_pre2_species$outcome_species <- sub("1","Female",train1_pre2_species$outcome_species)
compaired = list(c("Male", "Female"))

# Mean Decrease Accuracy是指在随机森林中，通过计算特征的重要度，来评估每个特征的重要程度。其重要度的计算是基于在每个随机森林的决策树中，每个特征点在随机化之前和之后所降低的准确性。Mean decrease accuracy是选择特征重要性的一种有效方法，可以帮助我们在各种机器学习问题中筛选出最重要的特征。
varImpPlot(train1.rf_species, main = "Top feature importance", n.var = optimal)
#write.table(train1.rf_species$confusion, file = "results/model_amplicon/Species_confusion_rf2.txt", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(train1.rf_species$confusion, paste0(opts$output, "Species_confusion_rf2.txt"), sep = "\t", quote = F, row.names = T, col.names = T)
imp_species = as.data.frame(round(importance(train1.rf_species), 2))
imp_species = imp_species[order(imp_species$MeanDecreaseAccuracy, decreasing = F),]
#write.table(imp_species, file = "results/model_RF/Species_imp_rf2.txt", sep = "\t", quote = F, row.names = T, col.names = T)
#write.table(imp_species, file = "results/model_amplicon/Species_imp_rf2.txt", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(imp_species, paste0(opts$output, "Species_imp_rf2.txt"), sep = "\t", quote = F, row.names = T, col.names = T)
