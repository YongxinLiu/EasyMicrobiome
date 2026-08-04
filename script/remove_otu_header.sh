#!/bin/bash
###############################################################################
# Script Name : remove_otu_header.sh
# Description : Remove "#OTU ID" from the first line of a FastSpar output file.
#
# Usage:
#   bash remove_otu_header.sh input.tsv output.txt
#
# Example:
#   bash remove_otu_header.sh \
#       Centenarians_correlation.tsv \
#       R_Centenarians2.txt
###############################################################################

# 检查参数
if [ $# -ne 2 ]; then
    echo "Usage: bash $0 <input.tsv> <output.txt>"
    exit 1
fi

INPUT=$1
OUTPUT=$2

# 判断输入文件是否存在
if [ ! -f "$INPUT" ]; then
    echo "Error: Input file '$INPUT' not found!"
    exit 1
fi

# 创建输出目录
mkdir -p "$(dirname "$OUTPUT")"

# 删除第一行的 #OTU ID
sed '1s/^#OTU ID//' "$INPUT" > "$OUTPUT"

echo "Done!"
echo "Input : $INPUT"
echo "Output: $OUTPUT"