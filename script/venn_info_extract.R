#!/usr/bin/env Rscript

# -------------------------------
# 自动检查并安装缺失包
# -------------------------------
packages <- c("optparse", "VennDiagram", "gtools", "grid")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}


############################
## 参数定义
############################

option_list <- list(
  make_option(c("-t", "--taxonomy"),
              type = "character",
              help = "taxonomy.spf file"),
  
  make_option(c("-m", "--metadata"),
              type = "character",
              help = "metadata.txt file"),
  
  make_option(c("-g", "--groups"),
              type = "character",
              help = "GroupID list, comma separated"),
  
  make_option(c("-k", "--k_common"),
              type = "character",
              default = "2",
              help = "k groups for common species, e.g. 2 or 2,3"),
  
  make_option(c("-o", "--outdir"),
              type = "character",
              default = "results",
              help = "Output directory"),
  
  make_option(c("--max_venn"),
              type = "integer",
              default = 4,
              help = "Max group number for Venn plot")
)

opt <- parse_args(OptionParser(option_list = option_list))

############################
## 参数解析
############################

groups <- unlist(strsplit(opt$groups, ","))
k_common <- as.integer(unlist(strsplit(opt$k_common, ",")))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

############################
## 读入数据
############################

data_s <- read.table(opt$taxonomy,
                     header = TRUE,
                     sep = "\t",
                     stringsAsFactors = FALSE)

metadata <- read.table(opt$metadata,
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE)

############################
## 筛选 metadata
############################

sub_metadata <- metadata[metadata$GroupID %in% groups, ]
rownames(sub_metadata) <- sub_metadata$SampleID

############################
## 物种矩阵处理
############################

# data_s2 <- data_s[, c(7, 9:ncol(data_s))]

if ("Species" %in% colnames(data_s)) {
  # Species 数据
  data_s2 <- data_s[, c(7, 9:ncol(data_s))]
} else {
  # KO 数据
  data_s2 <- data_s
}

# 统一第一列列名为 Feature
colnames(data_s2)[1] <- "Feature"

# data_s3 <- aggregate(. ~ Species,
#                      data = data_s2,
#                      FUN = sum)

data_s3 <- aggregate(. ~ Feature,
                     data = data_s2,
                     FUN = sum)

#rownames(data_s3) <- data_s3$Species
rownames(data_s3) <- data_s3$Feature
data_s3 <- data_s3[, -1]

data_s4 <- data_s3[, colnames(data_s3) %in% rownames(sub_metadata)]

############################
## 分组求平均
############################

data_s4_t <- as.data.frame(t(data_s4))
data_s4_t$GroupID <- sub_metadata[rownames(data_s4_t), "GroupID"]

group_mean <- aggregate(. ~ GroupID,
                        data = data_s4_t,
                        FUN = mean)

rownames(group_mean) <- group_mean$GroupID
group_mean <- group_mean[, -1]

############################
## 二值化
############################

bin_mat <- as.data.frame(t(group_mean))
bin_mat[] <- as.integer(bin_mat > 0)

############################
## Venn 图（<= max_venn）
############################

if (ncol(bin_mat) <= opt$max_venn) {
  
  venn_list <- lapply(colnames(bin_mat), function(g) {
    rownames(bin_mat)[bin_mat[, g] == 1]
  })
  names(venn_list) <- colnames(bin_mat)
  
  venn.plot <- venn.diagram(
    x = venn_list,
    filename = NULL,
    fill = c("#4DBBD5", "#E64B35", "#00A087", "#3C5488")[1:ncol(bin_mat)],
    alpha = 0.5,
    cex = 1.3,
    cat.cex = 1.2,
    cat.pos = seq(-30, 30, length.out = ncol(bin_mat)),
    cat.dist = rep(0.05, ncol(bin_mat)),
    main = "Feature Presence Venn"
  )
  
  pdf(file.path(opt$outdir, "Venn.pdf"), width = 6, height = 6)
  grid.draw(venn.plot)
  dev.off()
}

############################
## 特有 & 共有物种
############################

group_species <- lapply(colnames(bin_mat), function(g) {
  rownames(bin_mat)[bin_mat[, g] == 1]
})
names(group_species) <- colnames(bin_mat)

## 特有物种
for (g in names(group_species)) {
  
  others <- setdiff(names(group_species), g)
  
  unique_sp <- setdiff(
    group_species[[g]],
    unique(unlist(group_species[others]))
  )
  
  write.csv(
    data.frame(Feature = unique_sp),
    file = file.path(opt$outdir, paste0("Unique_", g, ".csv")),
    row.names = FALSE
  )
}

## 指定 k 组共有
for (k in k_common) {
  
  if (k > length(group_species)) next
  
  combs <- combinations(
    n = length(group_species),
    r = k,
    v = names(group_species)
  )
  
  for (i in 1:nrow(combs)) {
    
    gs <- combs[i, ]
    
    common_sp <- Reduce(intersect, group_species[gs])
    
    write.csv(
      data.frame(Feature = common_sp),
      file = file.path(
        opt$outdir,
        paste0("Common_", paste(gs, collapse = "_"), ".csv")
      ),
      row.names = FALSE
    )
  }
}
