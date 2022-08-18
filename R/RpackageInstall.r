site= "https://mirrors.tuna.tsinghua.edu.cn/CRAN"
old <- options(BioC_mirror=c("http://bioconductor.statistik.tu-dortmund.de/"))
#, "https://mirrors.nju.edu.cn/bioconductor/",  "http://mirrors.ustc.edu.cn/bioc/","https://mirrors.tuna.tsinghua.edu.cn/bioconductor"))


local({r = getOption("repos")  
r["CRAN"] = "http://mirrors.tuna.tsinghua.edu.cn/CRAN/"   
# r["BioC_mirror"] = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor"
options(repos=r)}) 


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

cran <- c("RColorBrewer", "gplots", "agricolae","optparse", "plotrix","igraph", 
"psych","sqldf","amap", "randomForest", "gridExtra", "reshape2", "ggplot2", 
"ggrepel", "pheatmap","ggbeeswarm","cowplot","plyr","stringr","grid","VennDiagram", 
"UpSetR","dplyr","showtext", "vegan", "knitr",  "scatterplot3d", "ggfortify", 
"gridExtra", "survival", "survminer", "RColorBrewer", "readr", "data.table", 
"WGCNA","devtools", "statmod", "mvoutlier", "mclust",  "penalized", "cluster", 
"KernSmooth", "mgcv", "ROCR",  "tidyverse", "ggthemes", "corrplot","BiocManager", 
"factoextra", "ggpubr",
"circlize", "vegan", "data.tree", "biomformat", "robCompositions", "multcompView", 
"scales", "devtools","Rcpp", "RcppArmadillo", "vegan", "reshape2", 
"gridExtra", "phyloseq",  "markovchain", "picante", "ggalluvial", 'huge', 'pulsar','VGAM', 'glmnet',
'rgexf')


a = rownames(installed.packages())

for(i in cran) {if(! i %in% a) BiocManager::install(i, update=F)}

if (!"FEAST" %in% a){
  devtools::install_github("cozygene/FEAST")
}

if (!"amplicon" %in% a){
    devtools::install_github("microbiota/amplicon")
}

if(!requireNamespace("microeco", quietly = T))
  install_github("ChiLiubio/microeco")

if(!requireNamespace("SpiecEasi", quietly = T))
  install_github("zdk123/SpiecEasi")


# 手动安装
# https://cran.r-project.org/web/packages/GUniFrac/index.html
# GUniFrac

options(old)
