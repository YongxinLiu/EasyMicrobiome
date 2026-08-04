#!/usr/bin/env Rscript

#############################
## Install & Load packages
#############################

# 检查并安装软件包
packages <- c(
  "optparse",
  "Hmisc",
  "igraph"
)

for(pkg in packages){
  
  if(!require(pkg,
              character.only = TRUE,
              quietly = TRUE)){
    
    install.packages(
      pkg,
      repos="https://cloud.r-project.org"
    )
    
    library(pkg,
            character.only = TRUE)
    
  }
  
}

# 定义输入参数
option_list <- list(
  
  make_option("--abundance",
              type="character",
              help="MAG abundance table"),
  
  make_option("--taxonomy",
              type="character",
              help="taxonomy table"),
  
  make_option("--metadata",
              type="character",
              help="metadata table"),
  
  make_option("--group1",
              default="Centenarians",
              type="character"),
  
  make_option("--group2",
              default="Young",
              type="character"),
  
  make_option("--abundance_cutoff",
              default=0.0007,
              type="double"),
  
  make_option("--prevalence_cutoff",
              default=0,
              type="double"),
  
  make_option("--outdir",
              default="result",
              type="character")
  
)

opt <- parse_args(
  OptionParser(option_list=option_list)
)

dir.create(
  opt$outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

# 读取数据

MAG_table_raw <-
  read.table(
    opt$abundance,
    header=TRUE,
    row.names=1,
    sep="\t",
    check.names=FALSE
  )

tax_table_raw <-
  read.table(
    opt$taxonomy,
    header=TRUE,
    row.names=1,
    sep="\t",
    check.names=FALSE
  )

meta_table_raw <-
  read.table(
    opt$metadata,
    header=TRUE,
    row.names=1,
    sep="\t",
    check.names=FALSE
  )


# 按组提取样本

MAG_table_group1 <-
  MAG_table_raw[,
                rownames(
                  meta_table_raw[
                    meta_table_raw$Group==opt$group1,])
  ]

MAG_table_group2 <-
  MAG_table_raw[,
                rownames(
                  meta_table_raw[
                    meta_table_raw$Group==opt$group2,])
  ]

MAG_tab_group1 <- t(MAG_table_group1)

MAG_tab_group2 <- t(MAG_table_group2)


# 计算Prevalence

prev_group1 <- colMeans(MAG_tab_group1>0)

prev_group2 <- colMeans(MAG_tab_group2>0)

write.table(
  data.frame(prevalence=prev_group1),
  file=file.path(opt$outdir,
                 "Prevalence_group1.tsv"),
  sep="\t",
  quote=FALSE
)

write.table(
  data.frame(prevalence=prev_group2),
  file=file.path(opt$outdir,
                 "Prevalence_group2.tsv"),
  sep="\t",
  quote=FALSE
)


# 绘制 Prevalence 曲线数据

calculate_prevalence <- function(data,group){
  
  threshold <- seq(0.01,1,0.01)
  
  prevalence <-
    sapply(
      threshold,
      function(x)
        sum(data>x)
    )
  
  data.frame(
    threshold,
    prevalence,
    group
  )
  
}

curve1 <-
  calculate_prevalence(
    prev_group1,
    opt$group1)

curve2 <-
  calculate_prevalence(
    prev_group2,
    opt$group2)

curve <-
  rbind(
    curve1,
    curve2)

write.table(
  curve,
  file=file.path(
    opt$outdir,
    "Prevalence_curve.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

# 计算 Relative abundance
total_relative_abundance <-
  colSums(t(MAG_table_raw))/
  sum(colSums(t(MAG_table_raw)))

write.table(
  data.frame(
    Relative_abundance=
      total_relative_abundance),
  file=file.path(
    opt$outdir,
    "Relative_abundance.tsv"),
  sep="\t",
  quote=FALSE
)

# Relative abundance 曲线
calculate_rela <- function(data){
  
  threshold <- seq(
    0,
    max(data),
    by=0.00025
  )
  
  species <-
    sapply(
      threshold,
      function(x)
        sum(data>x)
    )
  
  data.frame(
    threshold,
    species
  )
  
}

rela <-
  calculate_rela(
    total_relative_abundance)

write.table(
  rela,
  file=file.path(
    opt$outdir,
    "Relative_abundance_curve.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

# MAG筛选
is.above.prev <-
  prev_group1>=opt$prevalence_cutoff &
  prev_group2>=opt$prevalence_cutoff

is.above.abundance <-
  total_relative_abundance>=
  opt$abundance_cutoff

target_MAG <-
  rownames(MAG_table_raw)[
    is.above.prev &
      is.above.abundance
  ]

cat(
  "Selected MAG:",
  length(target_MAG),
  "\n"
)

# 输出筛选后的MAG
target_MAG_table <-
  MAG_table_raw[target_MAG,]

write.table(
  target_MAG_table,
  file=file.path(
    opt$outdir,
    "Selected_MAG.tsv"),
  sep="\t",
  quote=FALSE
)

# 输出FastSpar输入文件

write.table(
  MAG_table_group1[target_MAG,],
  file=file.path(
    opt$outdir,
    paste0(
      "MAG_for_FastSpar_",
      opt$group1,
      ".txt")),
  sep="\t",
  quote=FALSE
)

write.table(
  MAG_table_group2[target_MAG,],
  file=file.path(
    opt$outdir,
    paste0(
      "MAG_for_FastSpar_",
      opt$group2,
      ".txt")),
  sep="\t",
  quote=FALSE
)



