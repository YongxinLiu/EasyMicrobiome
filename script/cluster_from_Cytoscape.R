# #!/usr/bin/env Rscript
# 
# if (!requireNamespace("optparse", quietly = TRUE)) {
#   install.packages("optparse",
#                    repos = "https://cloud.r-project.org")
# }
# 
# library(optparse)
# 
# 
# # 检查并安装软件包
# 
# options(
#   repos = c(
#     CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
#   )
# )
# 
# options(
#   BioC_mirror =
#     "https://mirrors.tuna.tsinghua.edu.cn/bioconductor"
# )
# 
# # =====================================================
# # Install and load packages
# # =====================================================
# 
# cran_packages <- c(
#   "reshape2",
#   "ggplot2",
#   "dplyr",
#   "purrr",
#   "igraph"
# )
# 
# bioc_packages <- c(
#   "impute",
#   "preprocessCore",
#   "WGCNA"
# )
# 
# # ---------- CRAN ----------
# for(pkg in cran_packages){
#   
#   if(!requireNamespace(pkg, quietly = TRUE)){
#     
#     cat("Installing CRAN package:", pkg, "\n")
#     
#     install.packages(
#       pkg,
#       repos = "https://cloud.r-project.org"
#     )
#   }
#   
#   suppressPackageStartupMessages(
#     library(pkg,
#             character.only = TRUE)
#   )
# }
# 
# # ---------- Bioconductor ----------
# if(length(bioc_packages) > 0){
#   
#   if(!requireNamespace("BiocManager",
#                        quietly = TRUE)){
#     
#     install.packages(
#       "BiocManager",
#       repos = "https://cloud.r-project.org"
#     )
#   }
#   
#   for(pkg in bioc_packages){
#     
#     if(!requireNamespace(pkg,
#                          quietly = TRUE)){
#       
#       cat("Installing Bioconductor package:", pkg, "\n")
#       
#       BiocManager::install(
#         pkg,
#         ask = FALSE,
#         update = FALSE
#       )
#     }
#     
#     suppressPackageStartupMessages(
#       library(pkg,
#               character.only = TRUE)
#     )
#   }
# }
# 
# 
# #----
# #5.load ccluster from Cytoscape
# mysignode_stable_ccluster=read.table("./coverm/ccCluster.txt",header = T,sep = "\t",check.names = F,row.names = 1)
# mysignode_stable$CCcluster=mysignode_stable_ccluster$ccCluster[match(mysignode_stable$node,rownames(mysignode_stable_ccluster))]
# mysignode_stable$CCcluster=paste("C",mysignode_stable$CCcluster,sep = "")
# #5.1 find clusters in C1
# C1genome=mysignode_stable$node[which(mysignode_stable$CCcluster=="C1")]
# ind=which(mysigedge_stable$Source%in%C1genome&mysigedge_stable$Target%in%C1genome)
# tmpdata=graph.data.frame(mysigedge_stable[ind,],directed = F)
# tmpdata=get.adjacency(tmpdata,attr="NorPv2",sparse = F)
# tmpdist=1-tmpdata
# tmpcluster=hclust(as.dist(tmpdist),method = "average")
# hc=hclust(as.dist(tmpdist),method = "average")
# b=plot(hc,hang=-1,cex=0.8)
# mynamicmods=cutreeDynamic(dendro = hc,distM=tmpdist)
# table(mynamicmods)
# 
# mynamicmods=data.frame(ID=colnames(tmpdata),g=mynamicmods)
# mynamicmods$g=plyr::mapvalues(mynamicmods$g,from=c(0:25),to=LETTERS)
# mynamicmods$g=paste("C1",mynamicmods$g,sep = "")
# mynamicmods$ID_raw=rownames(myallGTDBTK)[match(mynamicmods$ID,myallGTDBTK$ID)]



#!/usr/bin/env Rscript

# =====================================================
# Install and load packages
# =====================================================

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages(
    "optparse",
    repos = "https://cloud.r-project.org"
  )
}

library(optparse)

options(
  repos = c(
    CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
  )
)

options(
  BioC_mirror =
    "https://mirrors.tuna.tsinghua.edu.cn/bioconductor"
)

cran_packages <- c(
  "reshape2",
  "ggplot2",
  "dplyr",
  "purrr",
  "igraph",
  "plyr"
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
    library(
      pkg,
      character.only = TRUE
    )
  )
}

# ---------- Bioconductor ----------
if(!requireNamespace(
  "BiocManager",
  quietly = TRUE
)){
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

for(pkg in bioc_packages){
  
  if(!requireNamespace(
    pkg,
    quietly = TRUE
  )){
    
    cat(
      "Installing Bioconductor package:",
      pkg,
      "\n"
    )
    
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
  
  suppressPackageStartupMessages(
    library(
      pkg,
      character.only = TRUE
    )
  )
}

# =====================================================
# Parameters
# =====================================================

option_list <- list(
  
  make_option(
    c("--ccCluster"),
    type = "character",
    default = "./coverm/ccCluster.txt",
    help = "ccCluster.txt"
  ),
  
  make_option(
    c("--node"),
    type = "character",
    default = "./coverm/372stable.node.txt",
    help = "stable.node2.txt"
  ),
  
  make_option(
    c("--edge"),
    type = "character",
    default = "./coverm/372stable.edge.txt",
    help = "stable.edge2.txt"
  ),
  
  make_option(
    c("--gtdb"),
    type = "character",
    default = "./coverm/GTDBTK.txt",
    help = "myallGTDBTK.txt"
  ),
  
  make_option(
    c("--cluster"),
    type = "character",
    default = "C1",
    help = "Cluster ID (default=C1)"
  ),
  
  make_option(
    c("-o","--output"),
    type = "character",
    default = "./coverm",
    help = "Output directory"
  )
)

opts <- parse_args(
  OptionParser(
    option_list = option_list
  )
)

dir.create(
  opts$output,
  recursive = TRUE,
  showWarnings = FALSE
)

# =====================================================
# Read files
# =====================================================

cat("Reading files ...\n")

mysignode_stable_ccluster <- read.table(
  opts$ccCluster,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  row.names = 1
)

mysignode_stable <- read.table(
  opts$node,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

mysigedge_stable <- read.table(
  opts$edge,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

myallGTDBTK <- read.table(
  opts$gtdb,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  row.names = 1
)

# =====================================================
# Merge cluster information
# =====================================================

mysignode_stable$CCcluster <-
  mysignode_stable_ccluster$ccCluster[
    match(
      mysignode_stable$node,
      rownames(mysignode_stable_ccluster)
    )
  ]

mysignode_stable$CCcluster <-
  paste0(
    "C",
    mysignode_stable$CCcluster
  )

# =====================================================
# Extract target cluster
# =====================================================

cat(
  "Processing cluster:",
  opts$cluster,
  "\n"
)

target_genome <-
  mysignode_stable$node[
    mysignode_stable$CCcluster ==
      opts$cluster
  ]

if(length(target_genome) < 2){
  
  stop(
    "Less than 2 genomes found in ",
    opts$cluster
  )
}

ind <- which(
  mysigedge_stable$Source %in%
    target_genome &
    mysigedge_stable$Target %in%
    target_genome
)

if(length(ind) == 0){
  
  stop(
    "No edges found in ",
    opts$cluster
  )
}

# =====================================================
# Build adjacency matrix
# =====================================================

cat("Building network ...\n")

tmpgraph <- graph.data.frame(
  mysigedge_stable[ind,],
  directed = FALSE
)

tmpdata <- get.adjacency(
  tmpgraph,
  attr = "NorPv2",
  sparse = FALSE
)

tmpdist <- 1 - tmpdata

# =====================================================
# Hierarchical clustering
# =====================================================

cat("Running hclust ...\n")

hc <- hclust(
  as.dist(tmpdist),
  method = "average"
)

# =====================================================
# Save clustering tree
# =====================================================

pdf(
  file.path(
    opts$output,
    paste0(
      opts$cluster,
      "_hclust_tree.pdf"
    )
  ),
  width = 10,
  height = 8
)

plot(
  hc,
  hang = -1,
  cex = 0.8,
  main = paste(
    opts$cluster,
    "Hierarchical Clustering"
  )
)

dev.off()

# =====================================================
# Dynamic Tree Cut
# =====================================================

cat(
  "Running dynamic tree cut ...\n"
)

mynamicmods <- cutreeDynamic(
  dendro = hc,
  distM = tmpdist
)

print(
  table(mynamicmods)
)

# =====================================================
# Generate module table
# =====================================================

module_letters <- c(
  LETTERS,
  paste0(
    rep(LETTERS, each = 26),
    LETTERS
  )
)

module_ids <-
  sort(unique(mynamicmods))

module_map <- setNames(
  module_letters[
    seq_along(module_ids)
  ],
  module_ids
)

mynamicmods <- data.frame(
  ID = colnames(tmpdata),
  Module = module_map[
    as.character(mynamicmods)
  ],
  stringsAsFactors = FALSE
)

mynamicmods$Module <-
  paste0(
    opts$cluster,
    "_",
    mynamicmods$Module
  )

if("ID" %in% colnames(myallGTDBTK)){
  
  mynamicmods$ID_raw <-
    rownames(myallGTDBTK)[
      match(
        mynamicmods$ID,
        myallGTDBTK$ID
      )
    ]
}

# =====================================================
# Save results
# =====================================================

write.table(
  mynamicmods,
  file = file.path(
    opts$output,
    paste0(
      opts$cluster,
      "_dynamic_modules.txt"
    )
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nFinished!\n"
)

cat(
  "Tree : ",
  file.path(
    opts$output,
    paste0(
      opts$cluster,
      "_hclust_tree.pdf"
    )
  ),
  "\n"
)

cat(
  "Module table : ",
  file.path(
    opts$output,
    paste0(
      opts$cluster,
      "_dynamic_modules.txt"
    )
  ),
  "\n"
)

