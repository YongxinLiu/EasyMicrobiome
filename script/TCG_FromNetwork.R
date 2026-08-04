#!/usr/bin/env Rscript

# Copyright 2026 De-feng Bai <baidefeng@caas.cn>

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse",
                   repos = "https://cloud.r-project.org")
}

library(optparse)

option_list <- list(
  
  make_option(
    c("-g", "--gtdb"),
    type = "character",
    help = "GTDBTK annotation file"
  ),
  
  make_option(
    "--group1_r",
    type = "character",
    help = "Group1 FastSpar correlation matrix"
  ),
  
  make_option(
    "--group1_p",
    type = "character",
    help = "Group1 FastSpar pvalue matrix"
  ),
  
  make_option(
    "--group2_r",
    type = "character",
    help = "Group2 FastSpar correlation matrix"
  ),
  
  make_option(
    "--group2_p",
    type = "character",
    help = "Group2 FastSpar pvalue matrix"
  ),
  
  make_option(c("-o", "--output"), type="character", default="coverm/",
              help="Output directory [default %default]")
)

opts <- parse_args(
  OptionParser(
    option_list = option_list
  )
)


# 检查并安装软件包

options(
  repos = c(
    CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
  )
)

options(
  BioC_mirror =
    "https://mirrors.tuna.tsinghua.edu.cn/bioconductor"
)

# =====================================================
# Install and load packages
# =====================================================

cran_packages <- c(
  "reshape2",
  "ggplot2",
  "dplyr",
  "purrr",
  "igraph"
)

bioc_packages <- c(
  "impute",
  "preprocessCore",
  "WGCNA"
)

# ---------- CRAN ----------
for(pkg in cran_packages){
  
  if(!requireNamespace(pkg, quietly = TRUE)){
    
    cat("Installing CRAN package:", pkg, "\n")
    
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
  
  suppressPackageStartupMessages(
    library(pkg,
            character.only = TRUE)
  )
}

# ---------- Bioconductor ----------
if(length(bioc_packages) > 0){
  
  if(!requireNamespace("BiocManager",
                       quietly = TRUE)){
    
    install.packages(
      "BiocManager",
      repos = "https://cloud.r-project.org"
    )
  }
  
  for(pkg in bioc_packages){
    
    if(!requireNamespace(pkg,
                         quietly = TRUE)){
      
      cat("Installing Bioconductor package:", pkg, "\n")
      
      BiocManager::install(
        pkg,
        ask = FALSE,
        update = FALSE
      )
    }
    
    suppressPackageStartupMessages(
      library(pkg,
              character.only = TRUE)
    )
  }
}

# library(reshape2)
# library(ggplot2)
# library(dplyr)
# library(purrr)
# library(ggplot2)
# library(igraph)
# library(WGCNA)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

#BiocManager::install("impute")

#BiocManager::install("preprocessCore")

#----
#0. local function
FastsparToABC=function(Fastspar_R,Fastspar_P)
{
  tmpdata=Fastspar_R
  tmpdata[lower.tri(tmpdata,diag = T)]=NA #set lower.tri and diag = NA
  tmpdata=melt(tmpdata) # melt the matrix
  tmpdata$node2=rep(colnames(Fastspar_R,length(colnames(Fastspar_R)))) # add node2 id
  tmpdata=tmpdata[which(is.na(tmpdata$value)==FALSE),]
  tmpR_abc=tmpdata
  tmpdata=Fastspar_P
  tmpdata[lower.tri(tmpdata,diag = T)]=NA #set lower.tri and diag = NA
  tmpdata=melt(tmpdata) # melt the matrix
  tmpdata$node2=rep(colnames(Fastspar_P,length(colnames(Fastspar_P)))) # add node2 id
  tmpdata=tmpdata[which(is.na(tmpdata$value)==FALSE),]
  tmpP_abc=tmpdata
  myresult=data.frame(Source=tmpP_abc$variable,Target=tmpP_abc$node2,FastsparR=tmpR_abc$value,FastsparP=tmpP_abc$value,stringsAsFactors = F)
  myresult$Source=as.character(myresult$Source)
  return(myresult)
}

FastsparABC_furprocess=function(myFastsparABC)
{
  myFastsparABC$Adjp=p.adjust(myFastsparABC$FastsparP,method = "BH")
  myFastsparABC$NorP="P"
  myFastsparABC$NorP[myFastsparABC$FastsparR<0]="N"
  myFastsparABC$R2=myFastsparABC$FastsparR^2
  for(i in 1:nrow(myFastsparABC))
  {
    if(myFastsparABC$Source[i]<myFastsparABC$Target[i])
    {
      myFastsparABC$link[i]=paste(myFastsparABC$Source[i],myFastsparABC$Target[i],sep = "/")
    }else{
      myFastsparABC$link[i]=paste(myFastsparABC$Target[i],myFastsparABC$Source[i],sep = "/")
    }
  }
  myFastsparABC$link_R=paste(myFastsparABC$link,myFastsparABC$NorP,sep = "=")
  return(myFastsparABC)
}

#----
#1. load GTDBTK result
# myallGTDBTK=read.table("GTDBTK2.txt",header = T,sep = "\t",row.names = 1,check.names = F)

#----
#2. load fastspar result
# myallBTFastsparR=read.table('R_Centenarians2.txt',header = T,row.names = 1,sep = "\t",check.names = F)
# myallBTFastsparP=read.table('P_Centenarians2.txt',header = T,row.names = 1,sep = "\t",check.names = F)
# myallMTFastsparR=read.table('R_Young2.txt',header = T,row.names = 1,sep = "\t",check.names = F)
# myallMTFastsparP=read.table('P_Young2.txt',header = T,row.names = 1,sep = "\t",check.names = F)


# =========================================================
# Read GTDBTK result
# =========================================================

myallGTDBTK <- read.table(
  opts$gtdb,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# =========================================================
# Read FastSpar results
# =========================================================

myallBTFastsparR <- read.table(
  opts$group1_r,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

myallBTFastsparP <- read.table(
  opts$group1_p,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

myallMTFastsparR <- read.table(
  opts$group2_r,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

myallMTFastsparP <- read.table(
  opts$group2_p,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)


myalledge=list()
myalledge$BT=FastsparToABC(Fastspar_R = myallBTFastsparR,Fastspar_P = myallBTFastsparP)
myalledge$BT=FastsparABC_furprocess(myalledge$BT)
myalledge$Tumor=FastsparToABC(Fastspar_R = myallMTFastsparR,Fastspar_P = myallMTFastsparP)
myalledge$Tumor=FastsparABC_furprocess(myalledge$Tumor)
class(myalledge$BT$FastsparR) # check! should be numeric
class(myalledge$BT$FastsparP) # check! should be numeric

#----
#3. get signficiant edges
mysigedge=c()
mytestedge=myalledge# 
mysigedge.list=list()
alledge.list=list() # for bar plot
for (i in 1:2)
{
  mytestedge[[i]]$NorP="P"
  mytestedge[[i]]$NorP[mytestedge[[i]]$FastsparR<0]="N"
  mytestedge[[i]]$NorP[which(mytestedge[[i]]$Adjp>=0.05)]="U"
  ind=which(mytestedge[[i]]$Adjp<0.05)
  print(length(ind))
  tmpedge=mytestedge[[i]][ind,]
  print(length(unique(c(tmpedge$Source,tmpedge$Target))))
  mysigedge=c(mysigedge,tmpedge$link_R)
  mysigedge.list[[i]]=tmpedge
  alledge.list[[i]]=mytestedge[[i]]
}
mysigedge_count=table(mysigedge)
table(mysigedge_count)

c=alledge.list%>%reduce(full_join,by="link")
d=c[,c("link","NorP.x","NorP.y")]
d[is.na(d)]="U"
e=table(paste(d$NorP.x,d$NorP.y,d$NorP,sep = ""))
e=data.frame(e)
e$Freq_ratio=e$Freq/sum(e$Freq)*100
e=e[order(e$Freq),]
e$Var1=factor(e$Var1,levels = e$Var1)
ggplot(data = e,aes(x=Var1,y=Freq))+
  geom_bar(stat = "identity",aes(fill=Var1))+
  coord_flip()+
  labs(x="",y="Number")+
  theme_bw()+
  theme(text=element_text(size=40))+theme(legend.position = "none")  #--> figure S2

#----
#4 get stable correlations
a=names(mysigedge_count)[which(mysigedge_count==2)] # =2 -> in all 2 groups
alink=gsub("=.*","",a)
mysigedge_stable=data.frame(
  Source=gsub("/.*","",alink),
  Target=gsub(".*/","",alink),
  NorP=gsub(".*=","",a),
  link=alink
)
mysigedge_stable$NorPv2=ifelse(mysigedge_stable$NorP=="P",1,-1)
mysignode_stable=data.frame(node=unique(c(mysigedge_stable$Source,mysigedge_stable$Target)))
mysignode_stable$GTDBTK=myallGTDBTK$classification[match(mysignode_stable$node,myallGTDBTK$ID)]

#write.table(mysigedge_stable,"372stable.edge2.txt",sep = "\t",quote = F,row.names = F)
#write.table(mysignode_stable,"372stable.node2.txt",sep = "\t",quote = F,row.names = F)

# write.table(mysigedge_stable,
#            file = "stable.edge2.txt",
#            output_dir = opts$output)
# 
# write.table(mysignode_stable,
#            file = "stable.node2.txt",
#            output_dir = opts$output)

write.table(
  mysigedge_stable,
  file = file.path(opts$output, "stable.edge2.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  mysignode_stable,
  file = file.path(opts$output, "stable.node2.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


