#!/bin/bash

# ============================================================
# 按 metadata 中的分组对 eggNOG COG_category 丰度表做组内均值聚合
#
# 用法：
#   bash extract_eggnog_anno02.sh -i <input> -m <metadata> -o <output>
# 示例：
#   bash extract_eggnog_anno02.sh \
#       -i eggnog/eggnog.COG_category.raw.txt \
#       -m ./metadata.txt \
#       -o eggnog/data_venn2.txt
# ============================================================

set -euo pipefail

usage() {
    echo ""
    echo "Usage:"
    echo "  bash $0 -i <input_abundance> -m <metadata> -o <output>"
    echo ""
    echo "Options:"
    echo "  -i, --input     输入丰度表 (第1列功能ID, 其余列为样本)"
    echo "  -m, --meta      样本元数据 (第1列 SampleID, 第3列 Group)"
    echo "  -o, --output    输出文件"
    echo "  -h, --help      显示帮助"
    echo ""
    echo "Example:"
    echo "  bash $0 -i eggnog/eggnog.COG_category.raw.txt \\"
    echo "          -m ./metadata.txt \\"
    echo "          -o eggnog/data_venn2.txt"
    echo ""
}

INPUT=""
META=""
OUTPUT=""

# ============================================================
# 1. 参数解析
# ============================================================
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)   INPUT="$2";  shift 2 ;;
        -m|--meta)    META="$2";   shift 2 ;;
        -o|--output)  OUTPUT="$2"; shift 2 ;;
        -h|--help)    usage; exit 0 ;;
        *)
            echo "[错误] 未知参数：$1" >&2
            usage
            exit 1
            ;;
    esac
done

# ============================================================
# 2. 参数检查
# ============================================================
if [[ -z "$INPUT" || -z "$META" || -z "$OUTPUT" ]]; then
    echo "[错误] 缺少必要参数 -i / -m / -o" >&2
    usage
    exit 1
fi

for f in "$INPUT" "$META"; do
    if [[ ! -f "$f" ]]; then
        echo "[错误] 文件不存在：$f" >&2
        exit 1
    fi
done

mkdir -p "$(dirname "$OUTPUT")"

# ============================================================
# 3. 核心处理
# ============================================================
awk '
BEGIN {
    FS = OFS = "\t"
}

############################
# 先读取 metadata
# 第1列 SampleID, 第3列 Group
############################
FNR == NR {
    # 兼容 BOM / CRLF
    sub(/^\xef\xbb\xbf/, "", $1)
    sub(/\r$/, "", $0)

    if (FNR > 1 && $1 != "" && $3 != "") {
        sample_group[$1] = $3
        groups[$3] = 1
    }
    next
}

############################
# 处理丰度表表头
############################
FNR == 1 {
    sub(/^\xef\xbb\xbf/, "", $1)

    for (i = 2; i <= NF; i++) {
        sample = $i
        sub(/\r$/, "", sample)
        group = sample_group[sample]
        col_group[i] = group
        if (group != "") group_n[group]++
    }

    # 输出表头（第一列取原始表头）
    printf "%s", $1

    for (g in groups) {
        group_order[++ng] = g
    }
    # 排序（gawk 特有 asort；非 gawk 环境见下方说明）
    n = asort(group_order)

    for (i = 1; i <= n; i++) {
        printf OFS "%s", group_order[i]
    }
    printf "\n"
    next
}

############################
# 处理每一行数据
############################
{
    sub(/\r$/, "", $0)

    printf "%s", $1
    delete sum

    for (i = 2; i <= NF; i++) {
        g = col_group[i]
        if (g != "") sum[g] += $i
    }

    for (i = 1; i <= n; i++) {
        g = group_order[i]
        mean = (group_n[g] > 0) ? sum[g] / group_n[g] : 0
        printf OFS "%.10g", mean
    }

    printf "\n"
}
' "$META" "$INPUT" > "$OUTPUT"

echo "Done! Output written to $OUTPUT"