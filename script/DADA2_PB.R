#!/usr/bin/env Rscript

# Title: DADA2 Pipeline for PacBio Full-Length 16S Data
# 标题: 用于PacBio全长16S数据的DADA2流程
# Author: luohao
# 作者: 罗豪
# Date: 2025-07-15
# 日期: 2025-07-15
# Description: This script processes multiple PacBio CCS samples,
#              including primer removal, filtering, denoising,
#              chimera removal, and ASV table generation.
# 描述: 该脚本用于处理多个PacBio CCS样本，包括去除引物、过滤、去噪、
#       去除嵌合体以及生成ASV表。

# -------------------------------------------------------------------
# 1. Load necessary packages
# 1. 加载必要的R包
# -------------------------------------------------------------------
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))

# -------------------------------------------------------------------
# 2. Set up command-line argument parsing
# 2. 设置命令行参数解析
# -------------------------------------------------------------------
parser <- ArgumentParser(description = "DADA2 pipeline for PacBio data. (用于PacBio数据的DADA2流程)")
parser$add_argument("--input_dir", type = "character", required = TRUE,
                    help = "Path to the directory containing input FASTQ files. (输入FASTQ文件所在目录的路径)")
parser$add_argument("--output_dir", type = "character", required = TRUE,
                    help = "Path to the directory where results will be saved. (用于保存结果的目录路径)")
parser$add_argument("--metadata_file", type = "character", required = TRUE,
                    help = "Path to the metadata file. Must contain a 'SampleID' column. (元数据文件路径，必须包含'SampleID'列)")
parser$add_argument("--fwd_primer", type = "character", default = "AGRGTTYGATYMTGGCTCAG",
                    help = "Forward primer sequence (e.g., 27F). (正向引物序列，例如27F)")
parser$add_argument("--rev_primer", type = "character", default = "RGYTACCTTGTTACGACTT",
                    help = "Reverse primer sequence (e.g., 1492R). (反向引物序列，例如1492R)")
parser$add_argument("--min_len", type = "integer", default = 1000,
                    help = "Minimum length of reads to keep after filtering. (过滤后保留序列的最小长度)")
parser$add_argument("--max_len", type = "integer", default = 1600,
                    help = "Maximum length of reads to keep after filtering. (过滤后保留序列的最大长度)")
parser$add_argument("--max_ee", type = "double", default = 2.0,
                    help = "Maximum expected errors allowed. (允许的最大期望误差)")
parser$add_argument("--threads", type = "integer", default = 4,
                    help = "Number of threads for multithreading. (用于多线程的线程数)")
parser$add_argument("--taxonomy_db", type = "character", required = TRUE,
                    help = "Path to the taxonomy database file (e.g., SILVA, GTDB). (物种注释数据库文件路径，例如SILVA、GTDB)")

args <- parser$parse_args()

# -------------------------------------------------------------------
# 3. Set up paths and environment
# 3. 设置路径和环境
# -------------------------------------------------------------------
cat("INFO: Setting up paths and reading metadata...\nINFO: 正在设置路径并读取元数据...\n")
path <- args$input_dir
path_out <- args$output_dir
metadata_path <- args$metadata_file

# Create output directories
# 创建输出目录
dir.create(path_out, showWarnings = FALSE, recursive = TRUE)
path_noprimers <- file.path(path_out, "noprimers")
dir.create(path_noprimers, showWarnings = FALSE)
path_filt <- file.path(path_out, "filtered")
dir.create(path_filt, showWarnings = FALSE)

# Read metadata and get sample names
# 读取元数据并获取样本名称
metadata <- read.delim(metadata_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
samples <- metadata$SampleID

# Define input file paths
# 定义输入文件路径
# Assuming files are named as SampleID.fastq.gz or SampleID.fastq
# 假设文件名为 SampleID.fastq.gz 或 SampleID.fastq
fastq_files <- file.path(path, paste0(samples, ".fastq.gz"))
if (!all(file.exists(fastq_files))) {
    fastq_files <- file.path(path, paste0(samples, ".fastq"))
    if (!all(file.exists(fastq_files))) {
        stop("Could not find all FASTQ files for the samples listed in metadata. (未能找到元数据中列出的所有样本的FASTQ文件)")
    }
}
names(fastq_files) <- samples

# -------------------------------------------------------------------
# Step 1: Remove Primers and Orient Reads
# 步骤 1: 去除引物并统一序列方向
# -------------------------------------------------------------------
cat("\nINFO: Step 1: Removing primers and orienting reads...\nINFO: 步骤 1: 正在去除引物并统一序列方向...\n")
noprimer_files <- file.path(path_noprimers, basename(fastq_files))
names(noprimer_files) <- samples

prim_out <- removePrimers(
  fastq_files,
  noprimer_files,
  primer.fwd = args$fwd_primer,
  primer.rev = dada2:::rc(args$rev_primer),
  orient = TRUE,
  verbose = TRUE
)
cat("INFO: Primer removal summary:\nINFO: 引物去除摘要:\n")
print(prim_out)

# -------------------------------------------------------------------
# Step 2: Filter and Trim
# 步骤 2: 过滤和修剪
# -------------------------------------------------------------------
cat("\nINFO: Step 2: Filtering and trimming reads...\nINFO: 步骤 2: 正在过滤和修剪序列...\n")
filt_files <- file.path(path_filt, basename(fastq_files))
names(filt_files) <- samples

filter_out <- filterAndTrim(
  noprimer_files,
  filt_files,
  minQ = 3,
  minLen = args$min_len,
  maxLen = args$max_len,
  maxN = 0,
  rm.phix = FALSE,
  maxEE = args$max_ee,
  verbose = TRUE,
  multithread = args$threads,
  compress = TRUE
)
cat("INFO: Filtering summary:\nINFO: 过滤摘要:\n")
print(filter_out)

# -------------------------------------------------------------------
# Step 3: Learn Error Model
# 步骤 3: 学习错误率模型
# -------------------------------------------------------------------
cat("\nINFO: Step 3: Learning error model...\nINFO: 步骤 3: 正在学习错误率模型...\n")
# Using PacBio-specific error function
# 使用PacBio特定的错误函数
err <- learnErrors(filt_files, multithread = args$threads, randomize = TRUE, errorEstimationFunction = dada2:::PacBioErrfun)
p_err <- plotErrors(err, nominalQ = TRUE)
ggsave(file.path(path_out, "error_model_plot.pdf"), p_err)
cat("INFO: Error model plot saved to", file.path(path_out, "error_model_plot.pdf"), "\nINFO: 错误率模型图已保存至", file.path(path_out, "error_model_plot.pdf"), "\n")

# -------------------------------------------------------------------
# Step 4: Denoise Samples
# 步骤 4: 对样本进行去噪
# -------------------------------------------------------------------
cat("\nINFO: Step 4: Denoising samples...\nINFO: 步骤 4: 正在对样本进行去噪...\n")
dada_list <- vector("list", length(samples))
names(dada_list) <- samples

for (sample in samples) {
  cat("  -> Processing sample:", sample, "\n  -> 正在处理样本:", sample, "\n")
  derep <- derepFastq(filt_files[[sample]])
  dada_list[[sample]] <- dada(derep, err = err, multithread = args$threads, BAND_SIZE = 32)
}

# -------------------------------------------------------------------
# Step 5: Construct ASV Table and Remove Chimeras
# 步骤 5: 构建ASV表并去除嵌合体
# -------------------------------------------------------------------
cat("\nINFO: Step 5: Constructing ASV table and removing chimeras...\nINFO: 步骤 5: 正在构建ASV表并去除嵌合体...\n")
seqtab <- makeSequenceTable(dada_list)
cat("INFO: Dimensions of ASV table before chimera removal:", dim(seqtab), "\nINFO: 去除嵌合体前ASV表的维度:", dim(seqtab), "\n")

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = args$threads, verbose = TRUE)
cat("INFO: Dimensions of ASV table after chimera removal:", dim(seqtab.nochim), "\nINFO: 去除嵌合体后ASV表的维度:", dim(seqtab.nochim), "\n")
cat("INFO: Proportion of non-chimeric reads:", sum(seqtab.nochim) / sum(seqtab), "\nINFO: 非嵌合体序列的比例:", sum(seqtab.nochim) / sum(seqtab), "\n")

# -------------------------------------------------------------------
# Step 6: Save Outputs
# 步骤 6: 保存输出结果
# -------------------------------------------------------------------
cat("\nINFO: Step 6: Saving final ASV table and sequences...\nINFO: 步骤 6: 正在保存最终的ASV表和序列...\n")
# Get sequences and assign standard names (ASV1, ASV2, ...)
# 获取序列并分配标准名称 (ASV1, ASV2, ...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- paste0("ASV", seq_along(asv_seqs))

# Save FASTA file
# 保存FASTA文件
asv_fasta_path <- file.path(path_out, "ASVs.fasta")
uniquesToFasta(seqtab.nochim, asv_fasta_path, ids = asv_headers)
cat("INFO: ASV sequences saved to:", asv_fasta_path, "\nINFO: ASV序列已保存至:", asv_fasta_path, "\n")

# Save ASV table
# 保存ASV表
asv_table <- as.data.frame(t(seqtab.nochim))
rownames(asv_table) <- asv_headers
colnames(asv_table) <- samples
asv_table_path <- file.path(path_out, "ASV_table.csv")
write.csv(asv_table, asv_table_path, row.names = TRUE)
cat("INFO: ASV count table saved to:", asv_table_path, "\nINFO: ASV计数表已保存至:", asv_table_path, "\n")

# -------------------------------------------------------------------
# Step 7: Assign Taxonomy
# 步骤 7: 物种注释
# -------------------------------------------------------------------
cat("\nINFO: Step 7: Assigning taxonomy...\nINFO: 步骤 7: 正在进行物种注释...\n")
taxa <- assignTaxonomy(seqtab.nochim, args$taxonomy_db, multithread = args$threads)
taxa_path <- file.path(path_out, "ASV_taxonomy.csv")
write.csv(taxa, taxa_path, row.names = TRUE)
cat("INFO: Taxonomy table saved to:", taxa_path, "\nINFO: 物种注释表已保存至:", taxa_path, "\n")

cat("\nINFO: DADA2 pipeline finished successfully!\nINFO: DADA2流程成功完成！\n")