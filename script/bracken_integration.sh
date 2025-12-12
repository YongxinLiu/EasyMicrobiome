#!/bin/bash

# -------------------------------
# Usage:
# bash bracken_integration.sh -t G -i temp/bracken/*.brk -o result/kraken2/bracken.G.txt
# -------------------------------

# 默认值
tax="G"
output_file="result/kraken2/bracken.${tax}.txt"
input_files=()

# 解析参数
while [[ $# -gt 0 ]]; do
  case $1 in
    -t|--tax)
      tax="$2"
      shift 2
      ;;
    -i|--input)
      shift
      # 收集后续所有参数直到遇到 -o 或结束
      while [[ $# -gt 0 && "$1" != "-o" ]]; do
        input_files+=("$1")
        shift
      done
      ;;
    -o|--output)
      output_file="$2"
      shift 2
      ;;
    *)
      echo "Unknown parameter: $1"
      exit 1
      ;;
  esac
done

# 检查输入文件是否存在
if [[ ${#input_files[@]} -eq 0 ]]; then
  echo "No input files provided!"
  exit 1
fi

# 创建输出目录
mkdir -p "$(dirname "$output_file")"

echo "Tax level: $tax"
echo "Number of input files: ${#input_files[@]}"
echo "Output file: $output_file"

# -------------------------------
# 收集所有tax信息
# -------------------------------
echo "Collecting all tax names..."
for f in "${input_files[@]}"; do
    awk 'NR>1 {print $1}' "$f"
done | sort -u > temp_bracken_all_names.data

# 输出第一列 Taxonomy
{
    echo "Taxonomy"
    cat temp_bracken_all_names.data
} > "$output_file"

# -------------------------------
# 遍历每个 brk 文件，提取 count
# -------------------------------
for f in "${input_files[@]}"; do
    sample=$(basename "$f" .brk)
    echo "Processing sample: $sample"

    # 建立 name -> count 对应表
    awk 'NR>1 {a[$1]=$6} END{for(k in a) print k"\t"a[k]}' "$f" > map.tmp

    # 当前样本列：表头 + 数据（按全集顺序）
    {
        echo "$sample"
        awk 'NR==FNR{m[$1]=$2; next} {print m[$1]+0}' \
            map.tmp temp_bracken_all_names.data
    } > col.tmp

    # 合并到总表
    paste "$output_file" col.tmp > tmp && mv tmp "$output_file"
done

# 清理临时文件
rm -f map.tmp col.tmp temp_bracken_all_names.data

echo "Done. Output saved to $output_file"

