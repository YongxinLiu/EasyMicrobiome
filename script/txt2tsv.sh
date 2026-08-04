#!/bin/bash

###############################################################################
# Script Name : txt2tsv.sh
# Description : Convert whitespace-delimited txt to tab-delimited tsv
#               and rename the first header to "#OTU ID".
#
# Usage:
#   bash txt2tsv.sh input.txt output.tsv
#
# Example:
#   bash txt2tsv.sh \
#       ./coverm/MAG_for_FastSpar_Centenarians.txt \
#       ./coverm/MAG_for_FastSpar_Centenarians_fixed.tsv
###############################################################################

# 检查参数
if [ $# -ne 2 ]; then
    echo "Usage: bash $0 <input.txt> <output.tsv>"
    exit 1
fi

INPUT=$1
OUTPUT=$2

# 判断输入文件是否存在
if [ ! -f "$INPUT" ]; then
    echo "Error: Input file '$INPUT' does not exist!"
    exit 1
fi

# 创建输出目录
mkdir -p "$(dirname "$OUTPUT")"

# 转换
awk '
BEGIN{
    OFS="\t"
}
{
    # 将连续空格压缩为TAB
    $1=$1

    # 第一行为表头
    if(NR==1){
        $1="#OTU ID"
    }

    print
}
' "$INPUT" > "$OUTPUT"

echo "Done!"
echo "Input : $INPUT"
echo "Output: $OUTPUT"