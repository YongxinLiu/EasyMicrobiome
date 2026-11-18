
#!/usr/bin/env bash

# ============================================================
# Bracken 结果整合脚本
#
# 使用方式：
# bash bracken_integration2.sh \
#     -t S \
#     -i temp/bracken/*.brk \
#     -o result/kraken2/bracken.S.txt
#
# ============================================================

set -euo pipefail

TAX=""
OUTPUT=""
FILES=()

# ============================================================
# 1. 参数解析
# ============================================================

usage() {
    echo ""
    echo "Usage:"
    echo "  bash $0 -t <P|G|S> -i <file1.brk> [file2.brk ...] -o <output>"
    echo ""
    echo "Example:"
    echo "  bash $0 -t S -i temp/bracken/*.brk -o result/kraken2/bracken.S.txt"
    echo ""
}

while [[ $# -gt 0 ]]; do

    case "$1" in

        -t)
            TAX="$2"
            shift 2
            ;;

        -i)
            shift

            # 收集 -i 后面的所有文件，直到遇到 -o
            while [[ $# -gt 0 && "$1" != "-o" ]]; do
                FILES+=("$1")
                shift
            done
            ;;

        -o)
            OUTPUT="$2"
            shift 2
            ;;

        -h|--help)
            usage
            exit 0
            ;;

        *)
            echo "[错误] 未知参数：$1"
            usage
            exit 1
            ;;
    esac

done

# ============================================================
# 2. 参数检查
# ============================================================

if [[ -z "$TAX" || -z "$OUTPUT" || ${#FILES[@]} -eq 0 ]]; then
    echo "[错误] 缺少必要参数"
    usage
    exit 1
fi

if [[ ! "$TAX" =~ ^(P|G|S)$ ]]; then
    echo "[错误] -t 仅支持 P、G、S"
    exit 1
fi

for file in "${FILES[@]}"; do
    if [[ ! -f "$file" ]]; then
        echo "[错误] 输入文件不存在：$file"
        exit 1
    fi
done

mkdir -p "$(dirname "$OUTPUT")"

# ============================================================
# 3. 创建临时文件
# ============================================================

TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

SAMPLE_FILE="$TMP_DIR/samples.txt"
COMBINED_FILE="$TMP_DIR/combined.tsv"
TAXID_FILE="$TMP_DIR/taxids.txt"

: > "$SAMPLE_FILE"
: > "$COMBINED_FILE"

echo "============================================================"
echo "Bracken integration started"
echo "Taxonomic level : $TAX"
echo "Input files     : ${#FILES[@]}"
echo "Output file     : $OUTPUT"
echo "============================================================"

# ============================================================
# 4. 提取指定分类水平
# ============================================================

for file in "${FILES[@]}"; do

    sample=$(basename "$file" .brk)

    echo "[处理] $sample"

    printf '%s\n' "$sample" >> "$SAMPLE_FILE"

    awk -F '\t' \
        -v tax="$TAX" \
        -v sample="$sample" \
        -v file="$file" '

    BEGIN {
        OFS = "\t"
    }

    NR == 1 {

        sub(/^\xef\xbb\xbf/, "", $1)

        for (i = 1; i <= NF; i++) {

            sub(/\r$/, "", $i)

            if ($i == "name") name_col = i
            if ($i == "taxonomy_id") taxid_col = i
            if ($i == "taxonomy_lvl") level_col = i
            if ($i == "new_est_reads") reads_col = i
        }

        if (!name_col || !taxid_col || !level_col || !reads_col) {
            printf "[错误] 文件 %s 缺少必要列\n", file > "/dev/stderr"
            exit 1
        }

        next
    }

    {

        sub(/\r$/, "", $0)

        if (NF < name_col ||
            NF < taxid_col ||
            NF < level_col ||
            NF < reads_col) {
            next
        }

        if ($(level_col) == tax) {

            taxid = $(taxid_col)
            name = $(name_col)
            reads = $(reads_col)

            if (taxid == "" || name == "") next

            if (reads !~ /^[0-9]+([.][0-9]+)?$/) next

            print taxid, name, sample, reads
        }
    }

    ' "$file" >> "$COMBINED_FILE"

done

# ============================================================
# 5. 检查提取结果
# ============================================================

if [[ ! -s "$COMBINED_FILE" ]]; then
    echo "[错误] 没有提取到分类水平 $TAX 的数据"
    echo "请检查 .brk 文件中的 taxonomy_lvl 列"
    exit 1
fi

# ============================================================
# 6. 提取 taxid
# ============================================================

cut -f1 "$COMBINED_FILE" |
    LC_ALL=C sort -u > "$TAXID_FILE"

# ============================================================
# 7. 构建丰度矩阵
# ============================================================

awk -F '\t' \
    -v samples_file="$SAMPLE_FILE" \
    -v data_file="$COMBINED_FILE" \
    -v taxids_file="$TAXID_FILE" '

BEGIN {

    OFS = "\t"

    # 读取样本名称
    while ((getline sample < samples_file) > 0) {
        if (sample != "") {
            sample_count++
            samples[sample_count] = sample
        }
    }

    close(samples_file)

    # 读取丰度数据
    while ((getline line < data_file) > 0) {

        split(line, arr, "\t")

        taxid = arr[1]
        name = arr[2]
        sample = arr[3]
        reads = arr[4]

        if (!(taxid in taxonomy_name)) {
            taxonomy_name[taxid] = name
        }

        abundance[taxid, sample] += reads
    }

    close(data_file)

    # 读取 taxid
    while ((getline taxid < taxids_file) > 0) {
        if (taxid != "") {
            taxid_count++
            taxids[taxid_count] = taxid
        }
    }

    close(taxids_file)

    # 输出表头
    printf "Taxonomy"

    for (i = 1; i <= sample_count; i++) {
        printf OFS "%s", samples[i]
    }

    printf "\n"

    # 输出丰度矩阵
    for (i = 1; i <= taxid_count; i++) {

        taxid = taxids[i]

        printf "%s", taxonomy_name[taxid]

        for (j = 1; j <= sample_count; j++) {

            sample = samples[j]
            value = abundance[taxid, sample] + 0

            printf OFS "%.10g", value
        }

        printf "\n"
    }
}
' > "$OUTPUT"

# ============================================================
# 8. 输出统计信息
# ============================================================

SAMPLE_COUNT=$(wc -l < "$SAMPLE_FILE")
TAXA_COUNT=$(wc -l < "$TAXID_FILE")

echo ""
echo "============================================================"
echo "Bracken integration completed!"
echo "============================================================"
echo "Taxonomic level : $TAX"
echo "Sample number   : $SAMPLE_COUNT"
echo "Taxa number     : $TAXA_COUNT"
echo "Output file     : $OUTPUT"
echo "============================================================"
