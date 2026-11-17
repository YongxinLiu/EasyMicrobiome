
#!/usr/bin/env Rscript

# ============================================================
# Copyright 2024-2026 Defeng Bai <baidefeng@caas.cn>
#
# EasyMetagenome
# MaAsLin3 differential abundance analysis
#
# Bai, et al. 2025. EasyMetagenome: A User-Friendly and
# Flexible Pipeline for Shotgun Metagenomic Analysis in
# Microbiome Research. iMeta 4: e70001.
# https://doi.org/10.1002/imt2.70001
#
# 功能：
# 1. 使用 MaAsLin3 进行微生物组多变量关联分析
# 2. 同时分析 abundance 和 prevalence
# 3. 支持固定效应和随机效应
# 4. 支持丰度、流行率、FDR 和 FC 统计
# 5. 输出总体差异分析结果
# 6. 输出 Group 分组的 enriched/depleted 结果
# 7. 保留 MaAsLin3 模型错误信息
# 8. 修复样本过滤后丰度矩阵与 Metadata 错位的问题
# 9. 区分模型不可估计与模型统计量缺失
#
# 更新：
# 2024/11/12：增加 MaAsLin2 方法
# 2025/11/27：规范代码
# 2026/05/22：更新软件包
# 2026/11/17：升级为 MaAsLin3
# ============================================================


# ============================================================
# 0. 清理环境
# ============================================================

rm(list = ls())

options(warn = 1)

set.seed(123)


# ============================================================
# 1. 加载命令行参数解析包
# ============================================================

if (!requireNamespace("optparse", quietly = TRUE)) {
  
  install.packages(
    "optparse",
    repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"
  )
  
}

suppressPackageStartupMessages(
  library(optparse)
)


# ============================================================
# 2. 参数设置
# ============================================================

option_list <- list(
  
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "result/metaphlan4/Species.txt",
    help = "Species relative abundance table [default %default]"
  ),
  
  make_option(
    c("-g", "--metadata"),
    type = "character",
    default = "result/metadata.txt",
    help = "Metadata file [default %default]"
  ),
  
  make_option(
    c("-a", "--min_abundance"),
    type = "numeric",
    default = 0,
    help = "Minimum abundance [default %default]"
  ),
  
  make_option(
    c("-p", "--min_prevalence"),
    type = "numeric",
    default = 0.1,
    help = "Minimum prevalence [default %default]"
  ),
  
  make_option(
    c("-s", "--max_significance"),
    type = "numeric",
    default = 0.05,
    help = "Maximum significance threshold [default %default]"
  ),
  
  make_option(
    c("-n", "--normalization"),
    type = "character",
    default = "TSS",
    help = "Normalization: TSS, CLR, NONE [default %default]"
  ),
  
  make_option(
    c("-t", "--transform"),
    type = "character",
    default = "LOG",
    help = "Transformation: LOG, PLOG, NONE [default %default]"
  ),
  
  make_option(
    c("-m", "--analysis_method"),
    type = "character",
    default = "LM",
    help = "Legacy MaAsLin2 parameter; retained for compatibility, ignored by MaAsLin3"
  ),
  
  make_option(
    c("-r", "--random_effects"),
    type = "character",
    default = "none",
    help = "Random effects, comma-separated [default %default]"
  ),
  
  make_option(
    c("-f", "--fixed_effects"),
    type = "character",
    default = "Group",
    help = "Fixed effects, comma-separated [default %default]"
  ),
  
  make_option(
    c("-c", "--correction"),
    type = "character",
    default = "BH",
    help = "Multiple testing correction [default %default]"
  ),
  
  make_option(
    c("-b", "--reference"),
    type = "character",
    default = "",
    help = "Reference levels, e.g. Group,Control;Sex,F [default %default]"
  ),
  
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "result/metaphlan4/",
    help = "Output directory [default %default]"
  )
  
)

opts <- parse_args(
  OptionParser(option_list = option_list)
)


# ============================================================
# 3. 参数检查
# ============================================================

if (!file.exists(opts$input)) {
  stop("输入丰度文件不存在：", opts$input)
}

if (!file.exists(opts$metadata)) {
  stop("Metadata 文件不存在：", opts$metadata)
}

if (opts$min_prevalence < 0 ||
    opts$min_prevalence > 1) {
  stop("min_prevalence 必须在 0-1 之间")
}

if (opts$min_abundance < 0) {
  stop("min_abundance 不能小于 0")
}

if (opts$max_significance <= 0 ||
    opts$max_significance > 1) {
  stop("max_significance 必须在 0-1 之间")
}

normalization <- toupper(opts$normalization)

transform_method <- toupper(opts$transform)

if (!normalization %in% c("TSS", "CLR", "NONE")) {
  
  stop(
    "MaAsLin3 normalization 仅支持 TSS、CLR、NONE；",
    "不支持 MaAsLin2 的 CSS 或 TMM"
  )
  
}

if (!transform_method %in% c("LOG", "PLOG", "NONE")) {
  
  stop(
    "MaAsLin3 transform 仅支持 LOG、PLOG、NONE；",
    "不支持 MaAsLin2 的 LOGIT 或 AST"
  )
  
}

cat(
  "MaAsLin3 使用内置线性模型和逻辑回归模型。\n",
  "原 analysis_method 参数：",
  opts$analysis_method,
  "，不再参与模型设置。\n"
)


# ============================================================
# 4. 安装和加载 MaAsLin3
# ============================================================

cran_repo <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN"

# ------------------------------------------------------------
# 4.1 安装 remotes
# ------------------------------------------------------------

if (!requireNamespace("remotes", quietly = TRUE)) {
  
  install.packages(
    "remotes",
    repos = cran_repo
  )
  
}

# ------------------------------------------------------------
# 4.2 从 GitHub 安装 MaAsLin3
# ------------------------------------------------------------

if (!requireNamespace("maaslin3", quietly = TRUE)) {
  
  cat("正在从 GitHub 安装 MaAsLin3，请稍候...\n")
  
  remotes::install_github(
    "biobakery/maaslin3",
    dependencies = TRUE,
    upgrade = "never"
  )
  
}

# ------------------------------------------------------------
# 4.3 检查安装结果
# ------------------------------------------------------------

if (!requireNamespace("maaslin3", quietly = TRUE)) {
  
  stop(
    "MaAsLin3 安装失败，请检查 GitHub 网络连接、",
    "R 版本及依赖包环境"
  )
  
}

# ------------------------------------------------------------
# 4.4 加载 MaAsLin3
# ------------------------------------------------------------

suppressPackageStartupMessages(
  library(maaslin3)
)

cat(
  "MaAsLin3 version: ",
  as.character(packageVersion("maaslin3")),
  "\n"
)

cat("MaAsLin3 加载成功！\n")


# ============================================================
# 5. 创建输出目录
# ============================================================

output_dir <- opts$output

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# MaAsLin3 独立输出目录
maaslin3_dir <- file.path(
  output_dir,
  "MaAsLin3_output"
)

dir.create(
  maaslin3_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 6. 读取 Metadata
# ============================================================

cat("\n========== 读取 Metadata ==========\n")

metadata <- read.delim(
  opts$metadata,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"SampleID" %in% colnames(metadata)) {
  stop("Metadata 中缺少 SampleID 列")
}

if (!"Group" %in% colnames(metadata)) {
  stop("Metadata 中缺少 Group 列")
}

metadata$SampleID <- as.character(
  metadata$SampleID
)

if (anyDuplicated(metadata$SampleID)) {
  stop("Metadata 中存在重复 SampleID")
}

rownames(metadata) <- metadata$SampleID

# SampleID 仅用于匹配，不作为模型协变量
metadata$SampleID <- NULL

# 将字符型 metadata 转换为分类变量
# 数值型变量保留为连续变量
for (col in colnames(metadata)) {
  
  if (is.character(metadata[[col]])) {
    
    metadata[[col]] <- factor(
      metadata[[col]],
      levels = unique(
        metadata[[col]][!is.na(metadata[[col]])]
      )
    )
    
  }
  
}

# Group 必须作为分类变量
metadata$Group <- factor(
  as.character(metadata$Group),
  levels = unique(
    as.character(
      metadata$Group[!is.na(metadata$Group)]
    )
  )
)

if (nlevels(metadata$Group) < 2) {
  stop("Group 至少需要包含两个有效分组")
}

cat("Metadata 样本数：", nrow(metadata), "\n")

cat(
  "Group 分组：",
  paste(levels(metadata$Group), collapse = " / "),
  "\n"
)


# ============================================================
# 7. 读取物种相对丰度表
# ============================================================

cat("\n========== 读取物种丰度表 ==========\n")

ra <- read.delim(
  opts$input,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"Taxonomy" %in% colnames(ra)) {
  stop("丰度表中缺少 Taxonomy 列")
}

if (anyDuplicated(ra$Taxonomy)) {
  stop("Taxonomy 中存在重复物种名称")
}

# 检查共同样本
common_samples <- intersect(
  rownames(metadata),
  colnames(ra)
)

if (length(common_samples) < 3) {
  stop("丰度表与 Metadata 的共同样本数不足")
}

cat("共同样本数：", length(common_samples), "\n")

# 按照 metadata 的顺序排列样本
common_samples <- rownames(metadata)[
  rownames(metadata) %in% common_samples
]

metadata <- metadata[
  common_samples,
  ,
  drop = FALSE
]

# 提取丰度矩阵
abundance_matrix <- as.matrix(
  ra[, common_samples, drop = FALSE]
)

# 转换为数值型
suppressWarnings(
  storage.mode(abundance_matrix) <- "numeric"
)

if (anyNA(abundance_matrix)) {
  stop("丰度表中存在 NA 或非数值数据，请检查输入文件")
}

if (any(!is.finite(abundance_matrix))) {
  stop("丰度表中存在 Inf 或 NaN")
}

if (any(abundance_matrix < 0)) {
  stop("丰度表中存在负数")
}

rownames(abundance_matrix) <- as.character(
  ra$Taxonomy
)

if (nrow(abundance_matrix) == 0 ||
    ncol(abundance_matrix) == 0) {
  stop("丰度矩阵为空")
}

# 转换为样本 x 物种格式
input_data_original <- as.data.frame(
  t(abundance_matrix),
  check.names = FALSE
)

# ============================================================
# 8. 处理物种名称
# ============================================================

# MaAsLin3 建模时使用合法且唯一的列名
# 保存原始物种名称用于结果还原

original_features <- colnames(input_data_original)

model_features <- make.names(
  original_features,
  unique = TRUE
)

feature_map <- data.frame(
  model_feature = model_features,
  original_feature = original_features,
  stringsAsFactors = FALSE
)

colnames(input_data_original) <- model_features

input_data <- input_data_original


# ============================================================
# 9. 设置固定效应和随机效应
# ============================================================

fixed_effects <- trimws(
  unlist(
    strsplit(
      opts$fixed_effects,
      ",",
      fixed = TRUE
    )
  )
)

fixed_effects <- fixed_effects[
  nzchar(fixed_effects)
]

if (length(fixed_effects) == 0) {
  stop("至少需要指定一个固定效应")
}

if (!all(fixed_effects %in% colnames(metadata))) {
  
  stop(
    "以下固定效应在 Metadata 中不存在：",
    paste(
      setdiff(fixed_effects, colnames(metadata)),
      collapse = ", "
    )
  )
  
}

random_effects <- character(0)

if (!tolower(opts$random_effects) %in% c("none", "")) {
  
  random_effects <- trimws(
    unlist(
      strsplit(
        opts$random_effects,
        ",",
        fixed = TRUE
      )
    )
  )
  
  random_effects <- random_effects[
    nzchar(random_effects)
  ]
  
}

if (!all(random_effects %in% colnames(metadata))) {
  
  stop(
    "以下随机效应在 Metadata 中不存在：",
    paste(
      setdiff(random_effects, colnames(metadata)),
      collapse = ", "
    )
  )
  
}

model_variables <- unique(
  c(fixed_effects, random_effects)
)

# ------------------------------------------------------------
# 9.1 删除模型变量中存在缺失值的样本
# ------------------------------------------------------------

complete_index <- complete.cases(
  metadata[, model_variables, drop = FALSE]
)

metadata <- metadata[
  complete_index,
  ,
  drop = FALSE
]

# 关键修复：
# 丰度矩阵必须同步过滤，并严格按照 metadata 的样本顺序排列。
# 后续所有组内丰度统计都使用这个过滤后的矩阵。
input_data <- input_data[
  rownames(metadata),
  ,
  drop = FALSE
]

if (nrow(metadata) < 3) {
  stop("删除缺失值后有效样本数不足")
}

# 删除未使用的因子水平
metadata <- droplevels(metadata)

if (nlevels(metadata$Group) < 2) {
  stop("删除缺失值后 Group 不足两个分组")
}

# 再次检查样本顺序
stopifnot(
  identical(
    rownames(input_data),
    rownames(metadata)
  )
)

cat("\n固定效应：\n")
print(fixed_effects)

cat("\n随机效应：\n")
print(random_effects)

cat("\n最终分析样本数：", nrow(metadata), "\n")

cat("\n最终各组样本数量：\n")
print(table(metadata$Group))


# ============================================================
# 10. 设置参考分组
# ============================================================

# 示例：
# -b "Group,Control"
# -b "Group,Control;Sex,F"
#
# MaAsLin3 对分类变量的系数以参考水平为比较基准

reference_string <- trimws(opts$reference)

if (nzchar(reference_string)) {
  
  reference_pairs <- strsplit(
    reference_string,
    ";",
    fixed = TRUE
  )[[1]]
  
  for (pair in reference_pairs) {
    
    pair_parts <- trimws(
      strsplit(
        pair,
        ",",
        fixed = TRUE
      )[[1]]
    )
    
    if (length(pair_parts) != 2) {
      stop("reference 参数格式错误，应为 Variable,Level")
    }
    
    variable_name <- pair_parts[1]
    
    reference_level <- pair_parts[2]
    
    if (!variable_name %in% colnames(metadata)) {
      stop("参考变量不存在：", variable_name)
    }
    
    if (!is.factor(metadata[[variable_name]])) {
      stop("参考变量不是分类变量：", variable_name)
    }
    
    if (!reference_level %in% levels(metadata[[variable_name]])) {
      
      stop(
        "参考水平不存在：",
        variable_name,
        " = ",
        reference_level
      )
      
    }
    
    # 在 metadata 中同步设置参考水平
    metadata[[variable_name]] <- stats::relevel(
      metadata[[variable_name]],
      ref = reference_level
    )
    
  }
  
}

cat("\nGroup 参考水平：", levels(metadata$Group)[1], "\n")

cat(
  "Group 比较水平：",
  paste(levels(metadata$Group)[-1], collapse = " / "),
  "\n"
)


# ============================================================
# 11. 运行 MaAsLin3
# ============================================================

cat("\n========== 开始 MaAsLin3 分析 ==========\n")

# ------------------------------------------------------------
# 11.1 组装 MaAsLin3 参数
# ------------------------------------------------------------

maaslin_args <- list(
  
  # 样本 x 物种丰度矩阵
  input_data = input_data,
  
  # 样本 metadata
  input_metadata = metadata,
  
  # 输出目录
  output = maaslin3_dir,
  
  # 固定效应
  fixed_effects = fixed_effects,
  
  # 随机效应
  random_effects = if (length(random_effects) > 0) {
    random_effects
  } else {
    NULL
  },
  
  # 参考水平
  reference = if (nzchar(reference_string)) {
    reference_string
  } else {
    NULL
  },
  
  # 丰度过滤
  min_abundance = opts$min_abundance,
  
  # 流行率过滤
  min_prevalence = opts$min_prevalence,
  
  # 标准化方法
  normalization = normalization,
  
  # 转换方法
  transform = transform_method,
  
  # 多重检验校正
  correction = opts$correction,
  
  # 显著性阈值
  max_significance = opts$max_significance,
  
  # 连续变量标准化
  standardize = TRUE,
  
  # 相对丰度组成效应校正
  median_comparison_abundance = TRUE,
  
  # Prevalence 不进行中位数校正
  median_comparison_prevalence = FALSE,
  
  # 生成汇总图
  plot_summary_plot = TRUE,
  
  # 生成关联图
  plot_associations = TRUE,
  
  # 最多绘制的关联图数量
  max_pngs = 30,
  
  # 日志输出等级
  verbosity = "WARN"
)

# ------------------------------------------------------------
# 11.2 根据当前 MaAsLin3 版本添加可选参数
# ------------------------------------------------------------

available_args <- names(
  formals(maaslin3::maaslin3)
)

# augment 对 prevalence logistic 模型可能有帮助
# 不能解决 abundance 模型中某组完全没有非零样本的问题
if ("augment" %in% available_args) {
  maaslin_args$augment <- TRUE
}

# 保存模型对象，便于后续诊断
if ("save_models" %in% available_args) {
  maaslin_args$save_models <- TRUE
}

# ------------------------------------------------------------
# 11.3 检查参数兼容性
# ------------------------------------------------------------

unknown_args <- setdiff(
  names(maaslin_args),
  available_args
)

if (length(unknown_args) > 0) {
  
  stop(
    "当前 MaAsLin3 版本不支持以下参数：",
    paste(unknown_args, collapse = ", ")
  )
  
}

# ------------------------------------------------------------
# 11.4 执行分析
# ------------------------------------------------------------

fit <- do.call(
  maaslin3::maaslin3,
  maaslin_args
)

cat("\nMaAsLin3 分析完成！\n")


# ============================================================
# 12. 读取 MaAsLin3 结果
# ============================================================

result_file <- file.path(
  maaslin3_dir,
  "all_results.tsv"
)

if (!file.exists(result_file)) {
  stop("未找到 MaAsLin3 输出文件：", result_file)
}

maaslin_results <- read.delim(
  result_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 检查必要字段
required_result_cols <- c(
  "feature",
  "metadata",
  "value",
  "coef",
  "stderr",
  "pval_individual",
  "qval_individual",
  "model"
)

if (!all(required_result_cols %in% colnames(maaslin_results))) {
  
  stop(
    "MaAsLin3 结果缺少必要字段：",
    paste(
      setdiff(
        required_result_cols,
        colnames(maaslin_results)
      ),
      collapse = ", "
    )
  )
  
}

# 兼容不同版本的 MaAsLin3
# 部分版本可能没有输出 error、N 或 N_not_zero 字段
optional_result_cols <- c(
  "error",
  "N",
  "N_not_zero"
)

for (col in optional_result_cols) {
  
  if (!col %in% colnames(maaslin_results)) {
    
    maaslin_results[[col]] <- NA
    
  }
  
}

# ------------------------------------------------------------
# 12.1 还原原始物种名称
# ------------------------------------------------------------

feature_index <- match(
  maaslin_results$feature,
  feature_map$model_feature
)

maaslin_results$feature_original <-
  feature_map$original_feature[feature_index]

# 检查是否存在无法匹配的物种名称
if (any(is.na(maaslin_results$feature_original))) {
  
  warning(
    "部分物种名称未能还原，将保留 MaAsLin3 输出的名称"
  )
  
  maaslin_results$feature_original[
    is.na(maaslin_results$feature_original)
  ] <- maaslin_results$feature[
    is.na(maaslin_results$feature_original)
  ]
  
}


# ============================================================
# 13. 计算 FC、置信区间和分组丰度
# ============================================================

cat("\n========== 汇总统计结果 ==========\n")

# ------------------------------------------------------------
# 13.1 关键修复：使用过滤后的丰度矩阵
# ------------------------------------------------------------

# input_data_original 是过滤缺失样本之前的矩阵。
# 如果直接使用它，可能导致丰度向量与 metadata 的分组索引长度不同。
#
# 此处必须使用 input_data，并按照 metadata 的样本顺序排列。
#
# input_data 的结构：
# 行 = 样本
# 列 = 物种

raw_abundance <- input_data[
  rownames(metadata),
  ,
  drop = FALSE
]

# 再次检查样本顺序
stopifnot(
  identical(
    rownames(raw_abundance),
    rownames(metadata)
  )
)

# ------------------------------------------------------------
# 13.2 定义安全数值转换函数
# ------------------------------------------------------------

get_numeric <- function(x) {
  
  suppressWarnings(
    as.numeric(as.character(x))
  )
  
}

# ------------------------------------------------------------
# 13.3 定义安全均值函数
# ------------------------------------------------------------

safe_mean <- function(x) {
  
  x <- as.numeric(x)
  
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  mean(x)
  
}

# ------------------------------------------------------------
# 13.4 定义安全非零样本数函数
# ------------------------------------------------------------

safe_nonzero <- function(x) {
  
  x <- as.numeric(x)
  
  sum(
    is.finite(x) & x > 0
  )
  
}

# ------------------------------------------------------------
# 13.5 初始化结果列表
# ------------------------------------------------------------

result_list <- vector(
  "list",
  nrow(maaslin_results)
)

# ------------------------------------------------------------
# 13.6 按模型结果逐行统计
# ------------------------------------------------------------

for (i in seq_len(nrow(maaslin_results))) {
  
  current <- maaslin_results[i, ]
  
  feature <- as.character(
    current$feature_original
  )
  
  variable <- as.character(
    current$metadata
  )
  
  model_type <- as.character(
    current$model
  )
  
  # ----------------------------------------------------------
  # 获取当前 Feature 对应的原始丰度
  # ----------------------------------------------------------
  
  feature_model_name <- feature_map$model_feature[
    match(
      feature,
      feature_map$original_feature
    )
  ]
  
  if (length(feature_model_name) == 0 ||
      is.na(feature_model_name)) {
    
    feature_abundance <- rep(
      NA_real_,
      nrow(metadata)
    )
    
  } else {
    
    feature_abundance <- as.numeric(
      raw_abundance[[feature_model_name]]
    )
    
  }
  
  # 检查丰度向量和 Metadata 样本数是否一致
  if (length(feature_abundance) != nrow(metadata)) {
    
    stop(
      "Feature 丰度向量与 Metadata 样本数不一致：",
      feature
    )
    
  }
  
  # 获取当前变量
  if (!variable %in% colnames(metadata)) {
    
    stop(
      "MaAsLin3 结果中的变量不在 Metadata 中：",
      variable
    )
    
  }
  
  variable_data <- metadata[[variable]]
  
  # ----------------------------------------------------------
  # 初始化统计量
  # ----------------------------------------------------------
  
  N1 <- NA_real_
  N2 <- NA_real_
  
  Mean1 <- NA_real_
  Mean2 <- NA_real_
  
  FC <- NA_real_
  FC_lower <- NA_real_
  FC_upper <- NA_real_
  
  Group_value <- NA_character_
  Reference <- NA_character_
  
  # ----------------------------------------------------------
  # 13.7 分类变量
  # ----------------------------------------------------------
  
  if (is.factor(variable_data)) {
    
    Group_value <- as.character(
      current$value
    )
    
    # 参考水平以实际 Metadata 因子水平为准
    Reference <- levels(variable_data)[1]
    
    # 只有非参考水平才进行组间统计
    if (Group_value %in% levels(variable_data) &&
        Group_value != Reference) {
      
      # 比较组
      group1_index <- variable_data == Group_value
      
      # 参考组
      group2_index <- variable_data == Reference
      
      # 避免 NA 逻辑值参与索引
      group1_index[is.na(group1_index)] <- FALSE
      
      group2_index[is.na(group2_index)] <- FALSE
      
      # ------------------------------------------------------
      # 组内非零样本数
      # ------------------------------------------------------
      
      N1 <- safe_nonzero(
        feature_abundance[group1_index]
      )
      
      N2 <- safe_nonzero(
        feature_abundance[group2_index]
      )
      
      # ------------------------------------------------------
      # 组内总样本的平均丰度
      # ------------------------------------------------------
      
      Mean1 <- safe_mean(
        feature_abundance[group1_index]
      )
      
      Mean2 <- safe_mean(
        feature_abundance[group2_index]
      )
      
    }
    
  } else {
    
    # --------------------------------------------------------
    # 13.8 连续变量
    # --------------------------------------------------------
    
    # 对于连续变量，保留原脚本的统计方式
    N1 <- safe_nonzero(feature_abundance)
    
    Mean1 <- safe_mean(feature_abundance)
    
  }
  
  # ----------------------------------------------------------
  # 13.9 提取模型统计量
  # ----------------------------------------------------------
  
  beta <- get_numeric(current$coef)
  
  se <- get_numeric(current$stderr)
  
  p_value <- get_numeric(current$pval_individual)
  
  fdr_value <- get_numeric(current$qval_individual)
  
  # ----------------------------------------------------------
  # 13.10 仅对 abundance 模型计算 FC
  # ----------------------------------------------------------
  
  # LOG 转换为 log2，因此 FC = 2^Beta
  # 如果使用 PLOG 或 NONE，不能直接将 2^Beta 解释为原始丰度倍数变化。
  
  if (tolower(model_type) == "abundance" &&
      transform_method == "LOG") {
    
    if (is.finite(beta)) {
      
      FC <- 2^beta
      
      if (is.finite(se) && se >= 0) {
        
        FC_lower <- 2^(beta - 1.96 * se)
        
        FC_upper <- 2^(beta + 1.96 * se)
        
      }
      
    }
    
  }
  
  # ----------------------------------------------------------
  # 13.11 汇总结果
  # ----------------------------------------------------------
  
  result_list[[i]] <- data.frame(
    
    Variable = variable,
    
    Model = model_type,
    
    Group_value = Group_value,
    
    Reference = Reference,
    
    Feature = feature,
    
    N1 = N1,
    
    N2 = N2,
    
    Mean1 = Mean1,
    
    Mean2 = Mean2,
    
    Beta = beta,
    
    SE = se,
    
    P = p_value,
    
    FDR = fdr_value,
    
    FC = FC,
    
    FC_lower = FC_lower,
    
    FC_upper = FC_upper,
    
    stringsAsFactors = FALSE,
    
    check.names = FALSE
    
  )
  
}

# 合并所有结果
res <- do.call(
  rbind,
  result_list
)

# 按照变量、模型和 P 值排序
res <- res[
  order(
    res$Variable,
    res$Model,
    res$P,
    na.last = TRUE
  ),
]

rownames(res) <- NULL


# ============================================================
# 14. 输出总体差异分析结果
# ============================================================

output_overall <- file.path(
  output_dir,
  "MaAsLin3_overall_difference.csv"
)

# 保留原有 16 列结构
# na = "NA" 可以避免缺失值以空白形式输出

write.csv(
  res,
  output_overall,
  row.names = FALSE,
  na = "NA"
)

cat(
  "总体差异分析结果：",
  output_overall,
  "\n"
)


# ============================================================
# 15. 输出 Group 的 enriched/depleted 结果
# ============================================================

# 仅提取 Group 的 abundance 模型
plot.data <- res[
  res$Variable == "Group" &
    tolower(res$Model) == "abundance",
  c(
    "Feature",
    "Group_value",
    "Reference",
    "P",
    "FDR",
    "FC"
  ),
  drop = FALSE
]

# 去除无法进行显著性判断的结果
plot.data <- plot.data[
  complete.cases(
    plot.data[, c("P", "FDR", "FC")]
  ),
  ,
  drop = FALSE
]

# 根据 FC 判断富集或耗竭
# FC > 1：Group_value 相对于 Reference 丰度升高
# FC < 1：Group_value 相对于 Reference 丰度降低

plot.data$`NPC association` <- ifelse(
  
  plot.data$FDR < opts$max_significance,
  
  ifelse(
    plot.data$FC > 1,
    "enriched",
    ifelse(
      plot.data$FC < 1,
      "depleted",
      "opposite directions"
    )
  ),
  
  "not associated"
  
)

# 输出
output_enriched <- file.path(
  output_dir,
  "MaAsLin3_enriched_depleted.csv"
)

write.csv(
  plot.data,
  output_enriched,
  row.names = FALSE,
  na = "NA"
)

cat(
  "Enriched/depleted 结果：",
  output_enriched,
  "\n"
)


# ============================================================
# 16. 输出模型诊断信息
# ============================================================

cat("\n========== 模型诊断 ==========\n")

# ------------------------------------------------------------
# 16.1 提取模型错误信息
# ------------------------------------------------------------

model_error <- as.character(
  maaslin_results$error
)

# 将空字符串和字符 NA 转换为真正的 NA
model_error[
  is.na(model_error) |
    model_error == "" |
    model_error == "NA"
] <- NA_character_

# ------------------------------------------------------------
# 16.2 识别模型状态
# ------------------------------------------------------------

diagnostic_beta <- get_numeric(
  maaslin_results$coef
)

diagnostic_se <- get_numeric(
  maaslin_results$stderr
)

diagnostic_p <- get_numeric(
  maaslin_results$pval_individual
)

diagnostic_fdr <- get_numeric(
  maaslin_results$qval_individual
)

fit_status <- rep(
  "OK",
  nrow(maaslin_results)
)

# 存在模型错误
has_error <- !is.na(model_error)

fit_status[has_error] <- "Model_error"

# Beta 缺失
idx <- fit_status == "OK" &
  !is.finite(diagnostic_beta)

fit_status[idx] <- "Not_estimable"

# SE 缺失
idx <- fit_status == "OK" &
  !is.finite(diagnostic_se)

fit_status[idx] <- "Missing_SE"

# P 值缺失
idx <- fit_status == "OK" &
  !is.finite(diagnostic_p)

fit_status[idx] <- "Missing_P"

# FDR 缺失
idx <- fit_status == "OK" &
  !is.finite(diagnostic_fdr)

fit_status[idx] <- "Missing_FDR"

# ------------------------------------------------------------
# 16.3 整理诊断表
# ------------------------------------------------------------

diagnostic_result <- data.frame(
  
  Variable = maaslin_results$metadata,
  
  Model = maaslin_results$model,
  
  Group_value = maaslin_results$value,
  
  Feature = maaslin_results$feature_original,
  
  Beta = diagnostic_beta,
  
  SE = diagnostic_se,
  
  P = diagnostic_p,
  
  FDR = diagnostic_fdr,
  
  Model_N = get_numeric(maaslin_results$N),
  
  Model_N_not_zero = get_numeric(
    maaslin_results$N_not_zero
  ),
  
  Fit_status = fit_status,
  
  Fit_error = model_error,
  
  stringsAsFactors = FALSE
  
)

# ------------------------------------------------------------
# 16.4 添加每组总样本数和非零样本数
# ------------------------------------------------------------

diagnostic_result$Group_N1 <- NA_real_

diagnostic_result$Group_N2 <- NA_real_

diagnostic_result$N1 <- NA_real_

diagnostic_result$N2 <- NA_real_

diagnostic_result$Mean1 <- NA_real_

diagnostic_result$Mean2 <- NA_real_

for (i in seq_len(nrow(diagnostic_result))) {
  
  variable <- as.character(
    diagnostic_result$Variable[i]
  )
  
  group_value <- as.character(
    diagnostic_result$Group_value[i]
  )
  
  feature <- as.character(
    diagnostic_result$Feature[i]
  )
  
  if (!variable %in% colnames(metadata)) {
    next
  }
  
  if (!feature %in% feature_map$original_feature) {
    next
  }
  
  model_feature <- feature_map$model_feature[
    match(
      feature,
      feature_map$original_feature
    )
  ]
  
  values <- as.numeric(
    raw_abundance[[model_feature]]
  )
  
  variable_data <- metadata[[variable]]
  
  if (is.factor(variable_data) &&
      group_value %in% levels(variable_data)) {
    
    reference_value <- levels(variable_data)[1]
    
    index1 <- variable_data == group_value
    
    index2 <- variable_data == reference_value
    
    index1[is.na(index1)] <- FALSE
    
    index2[is.na(index2)] <- FALSE
    
    diagnostic_result$Group_N1[i] <- sum(index1)
    
    diagnostic_result$Group_N2[i] <- sum(index2)
    
    diagnostic_result$N1[i] <- safe_nonzero(
      values[index1]
    )
    
    diagnostic_result$N2[i] <- safe_nonzero(
      values[index2]
    )
    
    diagnostic_result$Mean1[i] <- safe_mean(
      values[index1]
    )
    
    diagnostic_result$Mean2[i] <- safe_mean(
      values[index2]
    )
    
  }
  
}

# ------------------------------------------------------------
# 16.5 输出模型诊断文件
# ------------------------------------------------------------

output_diagnostic <- file.path(
  output_dir,
  "MaAsLin3_model_diagnostics.csv"
)

write.csv(
  diagnostic_result,
  output_diagnostic,
  row.names = FALSE,
  na = "NA"
)

# 单独输出模型存在问题的结果
problem_result <- diagnostic_result[
  diagnostic_result$Fit_status != "OK",
  ,
  drop = FALSE
]

write.csv(
  problem_result,
  file.path(
    output_dir,
    "MaAsLin3_model_problem.csv"
  ),
  row.names = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# 16.6 汇总模型错误信息
# ------------------------------------------------------------

error_summary <- diagnostic_result[
  !is.na(diagnostic_result$Fit_error),
  ,
  drop = FALSE
]

if (nrow(error_summary) > 0) {
  
  error_summary <- as.data.frame(
    table(
      Fit_status = error_summary$Fit_status,
      Fit_error = error_summary$Fit_error
    )
  )
  
  error_summary <- error_summary[
    error_summary$Freq > 0,
    ,
    drop = FALSE
  ]
  
  colnames(error_summary)[3] <- "Feature_number"
  
} else {
  
  error_summary <- data.frame(
    Fit_status = character(0),
    Fit_error = character(0),
    Feature_number = integer(0)
  )
  
}

write.csv(
  error_summary,
  file.path(
    output_dir,
    "MaAsLin3_error_summary.csv"
  ),
  row.names = FALSE,
  na = "NA"
)

cat("\n模型状态分布：\n")

print(
  table(
    diagnostic_result$Fit_status,
    useNA = "ifany"
  )
)


# ============================================================
# 17. 统计分析结果
# ============================================================

cat("\n========== 分析结果统计 ==========\n")

cat(
  "总关联结果数：",
  nrow(res),
  "\n"
)

cat(
  "Abundance 结果数：",
  sum(tolower(res$Model) == "abundance"),
  "\n"
)

cat(
  "Prevalence 结果数：",
  sum(tolower(res$Model) == "prevalence"),
  "\n"
)

cat(
  "显著 Abundance 结果数：",
  sum(
    tolower(res$Model) == "abundance" &
      !is.na(res$FDR) &
      res$FDR < opts$max_significance
  ),
  "\n"
)

cat(
  "显著 Prevalence 结果数：",
  sum(
    tolower(res$Model) == "prevalence" &
      !is.na(res$FDR) &
      res$FDR < opts$max_significance
  ),
  "\n"
)

cat(
  "Enriched 数量：",
  sum(
    plot.data$`NPC association` == "enriched",
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Depleted 数量：",
  sum(
    plot.data$`NPC association` == "depleted",
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "模型异常结果数：",
  nrow(problem_result),
  "\n"
)


# ============================================================
# 18. 保存 R 对象
# ============================================================

saveRDS(
  list(
    result_summary = res,
    enriched_depleted = plot.data,
    model_diagnostics = diagnostic_result,
    model_problem = problem_result,
    maaslin3_fit = fit,
    metadata = metadata,
    input_data = input_data,
    feature_map = feature_map
  ),
  file = file.path(
    output_dir,
    "MaAsLin3_analysis.rds"
  )
)


# ============================================================
# 19. 分析结束
# ============================================================

cat("\n========================================\n")
cat("MaAsLin3 差异分析全部完成！\n")
cat("========================================\n")

cat("输入丰度表：", opts$input, "\n")

cat("输入 Metadata：", opts$metadata, "\n")

cat("输出目录：", output_dir, "\n")

cat("MaAsLin3 原始结果：", result_file, "\n")

cat("总体结果：", output_overall, "\n")

cat("Enriched/depleted：", output_enriched, "\n")

cat(
  "模型诊断：",
  output_diagnostic,
  "\n"
)

cat("========================================\n")
