#!/usr/bin/env Rscript

# ============================================================
# Copyright 2024-2026 Defeng Bai <baidefeng@caas.cn>
#
# EasyMetagenome 风格的 MaAsLin3 差异分析脚本
#
# 输入文件：
#   1. Species.txt：第一列为 Taxonomy，后续列为样本丰度
#   2. metadata.txt：包含 SampleID 和分组/协变量信息
#
# 本脚本会自动将 Species.txt 转换为 MaAsLin3 所需的
# “样本为行、物种为列”的数据格式，并按 SampleID 匹配 metadata。
#
# 运行方式：
#   Rscript MaAsLin3_EasyMetagenome.R --help
#
# 示例：
#   Rscript MaAsLin3_EasyMetagenome.R \
#     -i result/metaphlan4/Species.txt \
#     -g result/metaphlan4/metadata.txt \
#     -o result/metaphlan4/MaAsLin3_result \
#     -f '~ Group' \
#     --reference 'Group,Young'
# ============================================================

# 清理当前 R 环境，避免已有对象影响本次分析
rm(list = ls())

# 设置随机数种子，尽可能保证分析过程可重复
set.seed(1)

# ============================================================
# 1. 检查并加载依赖包
# ============================================================

# 本脚本只需要 optparse 和 maaslin3，不再安装 MaAsLin2 及其他无关依赖
required_packages <- c("optparse", "maaslin3")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "缺少以下 R 包：", paste(missing_packages, collapse = ", "), "\n",
      "请先安装依赖包。\n",
      "安装 optparse：install.packages('optparse')\n",
      "安装 MaAsLin3：if (!requireNamespace('BiocManager', quietly=TRUE)) ",
      "install.packages('BiocManager'); ",
      "BiocManager::install('biobakery/maaslin3')"
    )
  )
}

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(maaslin3))

# ============================================================
# 2. 设置命令行参数
# ============================================================

# 所有参数均可在脚本运行时通过命令行修改。
# 路径参数支持相对路径和绝对路径；Windows 下建议使用正斜杠 /。

option_list <- list(
  make_option(
    c("-i", "--input"), type = "character",
    default = "result/metaphlan4/Species.txt",
    help = "物种丰度表，第一列为 Taxonomy，后续列为样本 [default: %default]"
  ),
  make_option(
    c("-g", "--metadata"), type = "character",
    default = "result/metadata.txt",
    help = "Metadata 文件，必须包含 SampleID 和模型所需变量 [default: %default]"
  ),
  make_option(
    c("-o", "--output"), type = "character",
    default = "result/metaphlan4/MaAsLin3_result",
    help = "MaAsLin3 结果输出目录 [default: %default]"
  ),
  make_option(
    c("--taxonomy_column"), type = "character",
    default = "Taxonomy",
    help = "丰度表中的物种名称列 [default: %default]"
  ),
  make_option(
    c("--sample_id_column"), type = "character",
    default = "SampleID",
    help = "Metadata 中的样本 ID 列 [default: %default]"
  ),
  make_option(
    c("-f", "--formula"), type = "character",
    default = "~ Group",
    help = "MaAsLin3 模型公式，例如 '~ Group' 或 '~ Group + age + sex' [default: %default]"
  ),
  make_option(
    c("--reference"), type = "character",
    default = "Group,Young",
    help = "分类变量参考水平，格式为 '变量名,参考水平'；多个变量用分号分隔 [default: %default]"
  ),
  make_option(
    c("-a", "--min_abundance"), type = "numeric",
    default = 0,
    help = "最小丰度过滤阈值 [default: %default]"
  ),
  make_option(
    c("-p", "--min_prevalence"), type = "numeric",
    default = 0.1,
    help = "最小流行率过滤阈值，范围通常为 0-1 [default: %default]"
  ),
  make_option(
    c("--max_prevalence"), type = "numeric",
    default = 1.01,
    help = "最大流行率过滤阈值 [default: %default]"
  ),
  make_option(
    c("--zero_threshold"), type = "numeric",
    default = 0,
    help = "小于或等于该值的丰度按 0 处理 [default: %default]"
  ),
  make_option(
    c("--min_variance"), type = "numeric",
    default = 0,
    help = "最小方差过滤阈值 [default: %default]"
  ),
  make_option(
    c("-s", "--max_significance"), type = "numeric",
    default = 0.05,
    help = "显著性阈值（FDR/q-value）[default: %default]"
  ),
  make_option(
    c("-n", "--normalization"), type = "character",
    default = "TSS",
    help = "标准化方法：TSS、CLR、CSS、NONE、TMM [default: %default]"
  ),
  make_option(
    c("-t", "--transform"), type = "character",
    default = "LOG",
    help = "数据转换方法：LOG、LOGIT、AST、NONE [default: %default]"
  ),
  make_option(
    c("-c", "--correction"), type = "character",
    default = "BH",
    help = "多重检验校正方法，例如 BH、holm、bonferroni [default: %default]"
  ),
  make_option(
    c("--median_comparison_abundance"), type = "logical",
    default = TRUE,
    help = "abundance 模型是否进行中位数比较 TRUE/FALSE [default: %default]"
  ),
  make_option(
    c("--median_comparison_prevalence"), type = "logical",
    default = FALSE,
    help = "prevalence 模型是否进行中位数比较 TRUE/FALSE [default: %default]"
  ),
  make_option(
    c("--standardize"), type = "logical",
    default = TRUE,
    help = "是否对连续协变量标准化 TRUE/FALSE [default: %default]"
  ),
  make_option(
    c("--augment"), type = "logical",
    default = TRUE,
    help = "是否对 prevalence 模型使用 augment TRUE/FALSE [default: %default]"
  ),
  make_option(
    c("--save_models"), type = "logical",
    default = TRUE,
    help = "是否保存模型对象 TRUE/FALSE [default: %default]"
  ),
  make_option(
    c("--warn_prevalence"), type = "logical",
    default = TRUE,
    help = "是否检查 prevalence 关联可能受丰度关联影响 TRUE/FALSE [default: %default]"
  ),
  make_option(
    c("--max_pngs"), type = "integer",
    default = 100,
    help = "最多输出多少张关联图 [default: %default]"
  )
)

# 解析用户输入的命令行参数
opt <- parse_args(OptionParser(option_list = option_list))

# ============================================================
# 3. 打印本次分析参数并检查输入文件
# ============================================================

cat("\n============================================================\n")
cat("MaAsLin3 差异分析参数\n")
cat("============================================================\n")
cat("丰度表：", opt$input, "\n")
cat("Metadata：", opt$metadata, "\n")
cat("输出目录：", opt$output, "\n")
cat("模型公式：", opt$formula, "\n")
cat("参考水平：", opt$reference, "\n")
cat("标准化方法：", opt$normalization, "\n")
cat("转换方法：", opt$transform, "\n")
cat("最小丰度：", opt$min_abundance, "\n")
cat("最小流行率：", opt$min_prevalence, "\n")
cat("显著性阈值：", opt$max_significance, "\n")
cat("============================================================\n\n")

if (!file.exists(opt$input)) {
  stop("找不到丰度表文件：", opt$input)
}

if (!file.exists(opt$metadata)) {
  stop("找不到 Metadata 文件：", opt$metadata)
}

# 创建输出目录；若目录已经存在则继续使用
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================
# 4. 读取物种丰度表
# ============================================================

# Species.txt 的预期格式：
# Taxonomy  C1  C2  C3  Y1  Y2  Y3
# species1  1.2 0   3.4 0   0   0.5
#
# 这里不将样本列名交给 R 自动修改，避免样本 ID 发生变化。

taxa_table <- read.delim(
  opt$input,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  #quote = "\""
  #comment.char = ""
  quote = '"',
  comment.char = ""
)

if (!(opt$taxonomy_column %in% colnames(taxa_table))) {
  stop(
    "丰度表中找不到物种名称列：", opt$taxonomy_column,
    "。当前列名为：", paste(colnames(taxa_table), collapse = ", ")
  )
}

# 提取物种名称，并检查空值和重复值
feature_names <- as.character(taxa_table[[opt$taxonomy_column]])

if (anyNA(feature_names) || any(feature_names == "")) {
  stop("丰度表的物种名称列存在 NA 或空字符串，请先检查输入文件。")
}

if (anyDuplicated(feature_names)) {
  duplicated_features <- unique(feature_names[duplicated(feature_names)])
  stop(
    "丰度表存在重复物种名称，请先合并或去重：",
    paste(duplicated_features, collapse = ", ")
  )
}

# 删除 Taxonomy 列，剩余列全部作为样本丰度
abundance_df <- taxa_table[
  , colnames(taxa_table) != opt$taxonomy_column,
  drop = FALSE
]

if (ncol(abundance_df) == 0) {
  stop("丰度表中没有检测到样本丰度列。")
}

# 将每个样本丰度列转换为 numeric，防止字符型数据进入模型
abundance_df[] <- lapply(
  abundance_df,
  function(x) suppressWarnings(as.numeric(as.character(x)))
)

abundance_matrix <- as.matrix(abundance_df)

if (anyNA(abundance_matrix)) {
  stop("丰度表中存在无法转换为数值的内容或 NA，请检查 Species.txt。")
}

if (any(!is.finite(abundance_matrix))) {
  stop("丰度表中存在 Inf 或 -Inf，请检查 Species.txt。")
}

if (any(abundance_matrix < 0)) {
  stop("丰度表中存在负数，MaAsLin3 不接受负丰度。")
}

# 设置物种名称为矩阵行名
rownames(abundance_matrix) <- feature_names

cat("丰度表读取成功：\n")
cat("  物种数量：", nrow(abundance_matrix), "\n")
cat("  样本数量：", ncol(abundance_matrix), "\n\n")

# ============================================================
# 5. 读取 Metadata
# ============================================================

metadata <- read.delim(
  opt$metadata,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  #quote = "\""
  #comment.char = ""
  quote = '"',
  comment.char = ""
)

if (!(opt$sample_id_column %in% colnames(metadata))) {
  stop(
    "Metadata 中找不到样本 ID 列：", opt$sample_id_column,
    "。当前列名为：", paste(colnames(metadata), collapse = ", ")
  )
}

sample_ids <- as.character(metadata[[opt$sample_id_column]])
original_metadata_sample_ids <- sample_ids

if (anyNA(sample_ids) || any(sample_ids == "")) {
  stop("Metadata 的样本 ID 列存在 NA 或空字符串。")
}

if (anyDuplicated(sample_ids)) {
  duplicated_samples <- unique(sample_ids[duplicated(sample_ids)])
  stop(
    "Metadata 中存在重复 SampleID：",
    paste(duplicated_samples, collapse = ", ")
  )
}

rownames(metadata) <- sample_ids
metadata[[opt$sample_id_column]] <- NULL

# 将字符型 metadata 转换为 factor，便于 MaAsLin3 识别分类变量。
# 数值型列（如 age、reads）保持 numeric，不进行转换。
for (column_name in colnames(metadata)) {
  if (is.character(metadata[[column_name]])) {
    metadata[[column_name]] <- factor(metadata[[column_name]])
  }
}

# ============================================================
# 6. 根据 --reference 设置分类变量参考水平
# ============================================================

# 参考水平参数格式示例：
#   --reference 'Group,Young'
# 多个变量示例：
#   --reference 'Group,Young;Sex,Female'
#
# 当使用 formula 时，MaAsLin3 会将分类变量的第一个 factor level
# 作为参考水平。因此这里直接调整 metadata 中 factor 的 levels 顺序。

if (!is.null(opt$reference) && nzchar(trimws(opt$reference))) {
  reference_items <- strsplit(opt$reference, ";", fixed = TRUE)[[1]]

  for (item in reference_items) {
    item <- trimws(item)
    pair <- strsplit(item, ",", fixed = TRUE)[[1]]

    if (length(pair) != 2) {
      stop(
        "--reference 参数格式错误：", item,
        "。正确格式示例：Group,Young；多个变量用分号分隔。"
      )
    }

    variable_name <- trimws(pair[1])
    reference_level <- trimws(pair[2])

    if (!(variable_name %in% colnames(metadata))) {
      stop("参考变量不存在于 Metadata 中：", variable_name)
    }

    # 参考变量转为 factor，并将用户指定水平放在第一个位置
    current_values <- as.character(metadata[[variable_name]])
    current_levels <- unique(current_values[!is.na(current_values)])

    if (!(reference_level %in% current_levels)) {
      stop(
        "指定的参考水平不存在：", variable_name, " = ", reference_level,
        "。该变量现有水平为：", paste(current_levels, collapse = ", ")
      )
    }

    new_levels <- c(reference_level, setdiff(current_levels, reference_level))
    metadata[[variable_name]] <- factor(current_values, levels = new_levels)
  }
}

# ============================================================
# 7. 匹配丰度表和 Metadata 样本
# ============================================================

# 仅保留两个文件中都存在的样本
original_abundance_sample_ids <- colnames(abundance_matrix)

common_samples <- intersect(
  original_abundance_sample_ids,
  rownames(metadata)
)

if (length(common_samples) == 0) {
  stop("丰度表和 Metadata 没有匹配的样本，请检查样本 ID 是否一致。")
}

# 按照丰度表原始样本顺序排列，保证矩阵与 metadata 完全对应
common_samples <- colnames(abundance_matrix)[
  colnames(abundance_matrix) %in% common_samples
]

abundance_matrix <- abundance_matrix[
  , common_samples,
  drop = FALSE
]

metadata <- metadata[
  common_samples,
  ,
  drop = FALSE
]

if (!identical(colnames(abundance_matrix), rownames(metadata))) {
  stop("样本匹配后丰度表与 Metadata 的样本顺序仍不一致。")
}

# MaAsLin3 要求样本为行、特征为列，因此在这里转置
input_data <- as.data.frame(
  t(abundance_matrix),
  check.names = FALSE
)

input_metadata <- metadata

if (!identical(rownames(input_data), rownames(input_metadata))) {
  stop("转置后丰度表与 Metadata 的样本顺序不一致。")
}

cat("样本匹配完成：\n")
cat("  匹配样本数：", nrow(input_data), "\n")
cat("  未匹配丰度表样本数：",
    length(setdiff(original_abundance_sample_ids, common_samples)), "\n")
cat("  未匹配 Metadata 样本数：",
    length(setdiff(original_metadata_sample_ids, common_samples)), "\n\n")

# ============================================================
# 8. 检查模型公式和变量
# ============================================================

model_formula <- tryCatch(
  as.formula(opt$formula),
  error = function(e) stop("模型公式无法解析：", conditionMessage(e))
)

# 提取公式中涉及的变量，并检查是否都存在于 Metadata
formula_variables <- all.vars(model_formula)
missing_formula_variables <- setdiff(formula_variables, colnames(input_metadata))

if (length(missing_formula_variables) > 0) {
  stop(
    "模型公式中以下变量不存在于 Metadata：",
    paste(missing_formula_variables, collapse = ", ")
  )
}

# 检查公式涉及的变量是否有缺失值
# 为避免模型内部因缺失值删样本后造成统计汇总错位，
# 这里先统一删除公式变量存在缺失值的样本。
if (length(formula_variables) > 0) {
  complete_rows <- complete.cases(
    input_metadata[, formula_variables, drop = FALSE]
  )

  if (any(!complete_rows)) {
    cat("由于模型变量存在缺失值，删除样本数：", sum(!complete_rows), "\n")
    input_metadata <- input_metadata[complete_rows, , drop = FALSE]
    input_data <- input_data[rownames(input_metadata), , drop = FALSE]
  }
}

if (nrow(input_data) < 2) {
  stop("删除缺失样本后有效样本数少于 2，无法进行分析。")
}

# 去除全零物种，避免无信息特征进入模型
nonzero_features <- colSums(input_data, na.rm = TRUE) > 0
input_data <- input_data[, nonzero_features, drop = FALSE]

if (ncol(input_data) == 0) {
  stop("过滤后没有非零物种，无法进行分析。")
}

# 保存原始样本-物种丰度矩阵，用于后续计算各组非零数和均值。
# 必须使用最终分析样本集合及其顺序，不能再使用未过滤的原始矩阵。
raw_abundance <- input_data

cat("最终进入模型的数据：\n")
cat("  样本数：", nrow(input_data), "\n")
cat("  物种数：", ncol(input_data), "\n")
cat("  模型公式：", paste(deparse(model_formula), collapse = " "), "\n\n")

# ============================================================
# 9. 运行 MaAsLin3
# ============================================================

cat("开始运行 MaAsLin3，请耐心等待...\n")

# 将参数集中管理，方便以后扩展和调整。
# min_abundance/min_prevalence 是 MaAsLin3 自身的特征过滤参数，
# 不等同于在原始丰度表上简单按比例过滤。
maaslin_args <- list(
  input_data = input_data,
  input_metadata = input_metadata,
  output = opt$output,
  formula = model_formula,
  min_abundance = opt$min_abundance,
  min_prevalence = opt$min_prevalence,
  max_prevalence = opt$max_prevalence,
  zero_threshold = opt$zero_threshold,
  min_variance = opt$min_variance,
  max_significance = opt$max_significance,
  normalization = toupper(opt$normalization),
  transform = toupper(opt$transform),
  correction = opt$correction,
  standardize = opt$standardize,
  median_comparison_abundance = opt$median_comparison_abundance,
  median_comparison_prevalence = opt$median_comparison_prevalence,
  warn_prevalence = opt$warn_prevalence,
  augment = opt$augment,
  max_pngs = opt$max_pngs,
  save_models = opt$save_models
)

# 根据当前安装版本的 maaslin3() 函数参数进行筛选，
# 避免不同版本间少数可选参数差异导致 unused argument 错误。
maaslin3_formals <- names(formals(maaslin3::maaslin3))
maaslin_args <- maaslin_args[names(maaslin_args) %in% maaslin3_formals]

fit_out <- do.call(maaslin3::maaslin3, maaslin_args)

cat("\nMaAsLin3 模型运行完成。\n")

# ============================================================
# 10. 检查并读取 MaAsLin3 原始结果
# ============================================================

all_results_file <- file.path(opt$output, "all_results.tsv")

if (!file.exists(all_results_file)) {
  stop(
    "MaAsLin3 运行结束，但没有找到 all_results.tsv：",
    all_results_file,
    "。请检查 MaAsLin3 日志文件和输出目录。"
  )
}

# 先读取文件开头，检查文件是否为空以及实际使用的分隔符。
# 正常的 MaAsLin3 all_results.tsv 应当是制表符分隔文件。
file_preview <- tryCatch(
  readLines(all_results_file, n = 5, warn = FALSE),
  error = function(e) character(0)
)

if (length(file_preview) == 0 || all(!nzchar(trimws(file_preview)))) {
  stop(
    "all_results.tsv 文件为空或无法读取：", all_results_file,
    "。请检查 MaAsLin3 是否成功完成模型拟合。"
  )
}

# 定义一个稳健的结果文件读取函数：
# 1. 首先按 TSV 读取；
# 2. 如果只读到一列，则尝试逗号或分号分隔；
# 3. 清除列名中的 UTF-8 BOM 和首尾空格；
# 4. 不擅自更改原始列名中的大小写。
read_maaslin_results <- function(file) {

  # 优先按 MaAsLin3 标准 TSV 格式读取
  result <- tryCatch(
    read.delim(
      file,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      #quote = "\""
      #comment.char = "",
      quote = '"',
      comment.char = "",
      fill = TRUE
    ),
    error = function(e) NULL
  )

  # 如果 TSV 读取后只有一列，说明文件可能不是标准制表符分隔
  if (is.null(result) || ncol(result) <= 1) {

    # 尝试逗号分隔
    result_comma <- tryCatch(
      read.csv(
        file,
        header = TRUE,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        quote = "\"",
        comment.char = "",
        fill = TRUE
      ),
      error = function(e) NULL
    )

    if (!is.null(result_comma) && ncol(result_comma) > 1) {
      result <- result_comma
    }
  }

  # 如果仍然只有一列，再尝试分号分隔
  if (is.null(result) || ncol(result) <= 1) {
    result_semicolon <- tryCatch(
      read.delim(
        file,
        header = TRUE,
        sep = ";",
        check.names = FALSE,
        stringsAsFactors = FALSE,
        quote = "\"",
        comment.char = "",
        fill = TRUE
      ),
      error = function(e) NULL
    )

    if (!is.null(result_semicolon) && ncol(result_semicolon) > 1) {
      result <- result_semicolon
    }
  }

  if (is.null(result)) {
    stop("无法解析 MaAsLin3 结果文件：", file)
  }

  # 清理列名中的 BOM 和空格，避免隐藏字符导致字段匹配失败
  clean_names <- trimws(names(result))
  clean_names <- sub("^\ufeff", "", clean_names)
  names(result) <- clean_names

  result
}

all_results <- read_maaslin_results(all_results_file)

# 兼容部分版本或不同导出方式下的列名大小写差异
names(all_results) <- tolower(trimws(names(all_results)))

# 兼容少数旧版或经过二次整理的结果文件中的常见字段别名
column_aliases <- c(
  "pval" = "pval_individual",
  "p_value" = "pval_individual",
  "pvalue" = "pval_individual",
  "qval" = "qval_individual",
  "q_value" = "qval_individual",
  "qvalue" = "qval_individual",
  "stderr" = "stderr",
  "standard_error" = "stderr",
  "std_error" = "stderr"
)

for (old_name in names(column_aliases)) {
  new_name <- unname(column_aliases[old_name])

  if (old_name %in% names(all_results) &&
      !(new_name %in% names(all_results))) {
    names(all_results)[names(all_results) == old_name] <- new_name
  }
}

# 检查必需字段
required_result_columns <- c(
  "feature", "metadata", "model", "coef", "stderr",
  "pval_individual", "qval_individual"
)

missing_result_columns <- setdiff(
  required_result_columns,
  colnames(all_results)
)

# 如果文件字段不完整，尝试使用 maaslin3() 返回对象中的结果表。
# 不同版本返回对象的结构可能略有差异，因此先检查对象结构。
if (length(missing_result_columns) > 0) {

  returned_results <- NULL

  if (is.list(fit_out) && "results" %in% names(fit_out)) {
    returned_results <- fit_out$results
  }

  if (is.data.frame(returned_results)) {
    names(returned_results) <- tolower(trimws(names(returned_results)))

    # 兼容返回对象中可能出现的字段别名
    for (old_name in names(column_aliases)) {
      new_name <- unname(column_aliases[old_name])

      if (old_name %in% names(returned_results) &&
          !(new_name %in% names(returned_results))) {
        names(returned_results)[
          names(returned_results) == old_name
        ] <- new_name
      }
    }

    returned_missing <- setdiff(
      required_result_columns,
      colnames(returned_results)
    )

    if (length(returned_missing) == 0) {
      all_results <- returned_results
      missing_result_columns <- character(0)
      message(
        "提示：all_results.tsv 字段不完整，已改用 fit_out$results 中的结果。"
      )
    }
  }
}

# 如果仍然无法识别，输出详细诊断信息，而不是只报告缺少字段
if (length(missing_result_columns) > 0) {

  cat("\n========== all_results.tsv 文件诊断 ==========\n")
  cat("文件路径：", normalizePath(
    all_results_file,
    winslash = "/",
    mustWork = FALSE
  ), "\n")

  cat("\n实际读取到的列名：\n")
  print(colnames(all_results))

  cat("\n文件前 5 行原始内容：\n")
  cat(paste(file_preview, collapse = "\n"), "\n")

  cat("\n缺少字段：\n")
  print(missing_result_columns)

  cat("==============================================\n")

  stop(
    "无法识别 MaAsLin3 结果文件结构。请检查上方打印的实际列名和文件前 5 行。"
  )
}

# 不同 MaAsLin3 版本可能没有某些可选诊断字段；
# 没有时补为 NA，保证后续诊断表可以正常生成。
optional_result_columns <- c(
  "value", "name", "error", "n", "n_not_zero"
)

for (column_name in optional_result_columns) {
  if (!(column_name %in% colnames(all_results))) {
    all_results[[column_name]] <- NA
  }
}

# 将结果中的特征名称映射回输入文件中的原始 Taxonomy 名称。
# MaAsLin3 或 R 在部分情况下会将特征名转换为合法变量名，
# 因此同时尝试原始名称和 make.names() 后的名称。
feature_map <- data.frame(
  feature_internal = make.names(colnames(raw_abundance), unique = TRUE),
  feature_original = colnames(raw_abundance),
  stringsAsFactors = FALSE
)

all_results$feature_original <- all_results$feature

match_internal <- match(all_results$feature, feature_map$feature_internal)
match_original <- match(all_results$feature, feature_map$feature_original)

use_internal <- !is.na(match_internal)
all_results$feature_original[use_internal] <-
  feature_map$feature_original[match_internal[use_internal]]

use_original <- !use_internal & !is.na(match_original)
all_results$feature_original[use_original] <-
  feature_map$feature_original[match_original[use_original]]

# 将数值字段转换为 numeric，兼容某些版本将数值读成字符的情况
numeric_result_columns <- c(
  "coef", "stderr", "pval_individual", "qval_individual", "n", "n_not_zero"
)

for (column_name in numeric_result_columns) {
  all_results[[column_name]] <- suppressWarnings(
    as.numeric(as.character(all_results[[column_name]]))
  )
}

# ============================================================
# 11. 生成总体差异结果表
# ============================================================

# 主结果表字段固定，方便后续 EasyMetagenome 其他模块调用。
# N1/N2 表示对应组中该物种的非零样本数，而不是组内总样本数。
# Mean1/Mean2 表示对应组的原始输入丰度均值。
# abundance + LOG 模型下 FC = 2^Beta。
# prevalence 模型的 Beta 是 log-odds，不将其误写为丰度 FC。

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

result_rows <- vector("list", nrow(all_results))

for (i in seq_len(nrow(all_results))) {
  one <- all_results[i, , drop = FALSE]

  variable_name <- as.character(one$metadata[1])
  feature_name <- as.character(one$feature_original[1])
  model_name <- as.character(one$model[1])
  group_value <- as.character(one$value[1])

  beta <- safe_numeric(one$coef[1])
  se <- safe_numeric(one$stderr[1])
  p_value <- safe_numeric(one$pval_individual[1])
  fdr_value <- safe_numeric(one$qval_individual[1])

  # 初始化组间统计量
  n1 <- n2 <- mean1 <- mean2 <- NA_real_
  reference_value <- NA_character_
  group_n1 <- group_n2 <- NA_integer_

  # 对分类变量，按当前 factor 的第一个水平作为参考组
  if (variable_name %in% colnames(input_metadata)) {
    variable_data <- input_metadata[[variable_name]]

    if (is.factor(variable_data) || is.character(variable_data)) {
      variable_data <- as.character(variable_data)
      observed_levels <- if (is.factor(input_metadata[[variable_name]])) {
        levels(input_metadata[[variable_name]])
      } else {
        unique(variable_data[!is.na(variable_data)])
      }

      observed_levels <- observed_levels[observed_levels %in% variable_data]

      if (length(observed_levels) >= 2 &&
          !is.na(group_value) && group_value %in% observed_levels) {
        reference_value <- observed_levels[1]

        group1_samples <- rownames(input_metadata)[
          !is.na(variable_data) & variable_data == group_value
        ]
        group2_samples <- rownames(input_metadata)[
          !is.na(variable_data) & variable_data == reference_value
        ]

        if (feature_name %in% colnames(raw_abundance)) {
          feature_abundance <- raw_abundance[, feature_name]

          values1 <- feature_abundance[group1_samples]
          values2 <- feature_abundance[group2_samples]

          group_n1 <- length(group1_samples)
          group_n2 <- length(group2_samples)
          n1 <- sum(values1 > 0, na.rm = TRUE)
          n2 <- sum(values2 > 0, na.rm = TRUE)
          mean1 <- if (length(values1) > 0) mean(values1, na.rm = TRUE) else NA_real_
          mean2 <- if (length(values2) > 0) mean(values2, na.rm = TRUE) else NA_real_
        }
      }
    }
  }

  # 只有 LOG 转换的 abundance 模型，Beta 才能直接换算为 2^Beta 的倍数变化
  fc <- fc_lower <- fc_upper <- NA_real_

  if (identical(tolower(model_name), "abundance") &&
      identical(toupper(opt$transform), "LOG") &&
      is.finite(beta)) {
    fc <- 2^beta

    if (is.finite(se)) {
      fc_lower <- 2^(beta - 1.96 * se)
      fc_upper <- 2^(beta + 1.96 * se)
    }
  }

  result_rows[[i]] <- data.frame(
    Variable = variable_name,
    Model = model_name,
    Group_value = group_value,
    Reference = reference_value,
    Feature = feature_name,
    N1 = n1,
    N2 = n2,
    Mean1 = mean1,
    Mean2 = mean2,
    Beta = beta,
    SE = se,
    P = p_value,
    FDR = fdr_value,
    FC = fc,
    FC_lower = fc_lower,
    FC_upper = fc_upper,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

if (length(result_rows) > 0) {
  overall_result <- do.call(rbind, result_rows)
} else {
  overall_result <- data.frame(
    Variable = character(), Model = character(), Group_value = character(),
    Reference = character(), Feature = character(), N1 = numeric(),
    N2 = numeric(), Mean1 = numeric(), Mean2 = numeric(), Beta = numeric(),
    SE = numeric(), P = numeric(), FDR = numeric(), FC = numeric(),
    FC_lower = numeric(), FC_upper = numeric(),
    stringsAsFactors = FALSE
  )
}

# 保证主结果表的列名和顺序固定
expected_columns <- c(
  "Variable", "Model", "Group_value", "Reference", "Feature",
  "N1", "N2", "Mean1", "Mean2", "Beta", "SE", "P", "FDR",
  "FC", "FC_lower", "FC_upper"
)

overall_result <- overall_result[, expected_columns, drop = FALSE]

# 按 FDR、P 值排序；NA 排在最后
if (nrow(overall_result) > 0) {
  overall_result <- overall_result[
    order(overall_result$FDR, overall_result$P, na.last = TRUE),
    , drop = FALSE
  ]
}

overall_file <- file.path(opt$output, "MaAsLin3_overall_difference.csv")
write.csv(
  overall_result,
  file = overall_file,
  row.names = FALSE,
  na = "NA",
  fileEncoding = "UTF-8"
)

cat("总体差异结果已输出：", overall_file, "\n")

# ============================================================
# 12. 输出模型诊断信息
# ============================================================

# 将 MaAsLin3 原始结果中的 error、N、N_not_zero 保留下来，
# 便于判断 NA 是模型不可估计、样本稀疏还是模型拟合错误导致。

diagnostics <- data.frame(
  Variable = all_results$metadata,
  Model = all_results$model,
  Group_value = all_results$value,
  Feature = all_results$feature_original,
  Beta = all_results$coef,
  SE = all_results$stderr,
  P = all_results$pval_individual,
  FDR = all_results$qval_individual,
  MaAsLin3_N = all_results$n,
  MaAsLin3_N_not_zero = all_results$n_not_zero,
  Fit_error = all_results$error,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 根据模型输出字段标记拟合状态，仅用于诊断，不修改统计结果
# 不将无法估计的 NA 强行填成 0 或其他数值。
diagnostics$Fit_status <- ifelse(
  !is.na(diagnostics$Fit_error) & diagnostics$Fit_error != "",
  "Model_error",
  ifelse(
    is.na(diagnostics$Beta), "Not_estimable",
    ifelse(
      is.na(diagnostics$SE), "Missing_SE",
      ifelse(
        is.na(diagnostics$P), "Missing_P",
        ifelse(is.na(diagnostics$FDR), "Missing_FDR", "OK")
      )
    )
  )
)

diagnostics_file <- file.path(
  opt$output,
  "MaAsLin3_model_diagnostics.csv"
)

write.csv(
  diagnostics,
  file = diagnostics_file,
  row.names = FALSE,
  na = "NA",
  fileEncoding = "UTF-8"
)

cat("模型诊断信息已输出：", diagnostics_file, "\n")

# ============================================================
# 13. 输出 enriched / depleted 结果
# ============================================================

# 仅使用 abundance 模型，并且只保留 P、FDR、FC 均可用的结果。
# FC > 1：比较组相对于参考组丰度升高（enriched）。
# FC < 1：比较组相对于参考组丰度降低（depleted）。

plot_data <- overall_result[
  !is.na(overall_result$Model) & tolower(overall_result$Model) == "abundance",
  c("Variable", "Group_value", "Reference", "Feature", "P", "FDR", "FC"),
  drop = FALSE
]

plot_data <- plot_data[
  is.finite(plot_data$P) &
    is.finite(plot_data$FDR) &
    is.finite(plot_data$FC),
  , drop = FALSE
]

plot_data$Association <- ifelse(
  plot_data$FDR <= opt$max_significance & plot_data$FC > 1,
  "enriched",
  ifelse(
    plot_data$FDR <= opt$max_significance & plot_data$FC < 1,
    "depleted",
    "not_significant"
  )
)

plot_data <- plot_data[
  order(plot_data$FDR, plot_data$P, na.last = TRUE),
  , drop = FALSE
]

plot_file <- file.path(
  opt$output,
  "MaAsLin3_enriched_depleted.csv"
)

write.csv(
  plot_data,
  file = plot_file,
  row.names = FALSE,
  na = "NA",
  fileEncoding = "UTF-8"
)

cat("Enriched/depleted 结果已输出：", plot_file, "\n")

# ============================================================
# 14. 输出分析统计汇总
# ============================================================

cat("\n============================================================\n")
cat("MaAsLin3 分析完成\n")
cat("============================================================\n")
cat("分析样本数：", nrow(input_data), "\n")
cat("分析物种数：", ncol(input_data), "\n")
cat("MaAsLin3 结果行数：", nrow(all_results), "\n")
cat("总体差异结果：", overall_file, "\n")
cat("模型诊断结果：", diagnostics_file, "\n")
cat("Enriched/depleted 结果：", plot_file, "\n")

cat("\n模型拟合状态统计：\n")
print(table(diagnostics$Fit_status, useNA = "ifany"))

cat("\n输出目录：", normalizePath(opt$output, winslash = "/", mustWork = FALSE), "\n")
cat("============================================================\n")

# 保存本次运行的关键对象，便于后续在 R 中继续分析
saveRDS(
  list(
    options = opt,
    input_data = input_data,
    input_metadata = input_metadata,
    fit_out = fit_out,
    all_results = all_results,
    overall_result = overall_result,
    diagnostics = diagnostics
  ),
  file = file.path(opt$output, "MaAsLin3_analysis.rds")
)

cat("分析对象已保存：", file.path(opt$output, "MaAsLin3_analysis.rds"), "\n")
