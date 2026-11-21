#!/bin/bash

# ============================================================
# 提取 eggNOG 注释信息，生成 Venn 分析用的存在/缺失矩阵
#
# 用法：
#   bash extract_eggnog_anno01.sh -i <input> -o <output>
# 示例：
#   bash extract_eggnog_anno01.sh \
#       -i ../temp/eggnog/output \
#       -o eggnog/data_venn1.txt
# ============================================================

set -euo pipefail

usage() {
    echo ""
    echo "Usage:"
    echo "  bash $0 -i <input_eggnog_file> -o <output_file>"
    echo ""
    echo "Example:"
    echo "  bash $0 -i ../temp/eggnog/output -o eggnog/data_venn1.txt"
    echo ""
}

INPUT=""
OUTPUT=""

# ============================================================
# 1. 参数解析
# ============================================================
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            INPUT="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
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
if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
    echo "[错误] 缺少必要参数 -i/-o" >&2
    usage
    exit 1
fi

if [[ ! -f "$INPUT" ]]; then
    echo "[错误] 输入文件不存在：$INPUT" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT")"

# ============================================================
# 3. 核心处理
# ============================================================
awk -F'\t' '
BEGIN {
    OFS = "\t"
}

# 跳过注释行（eggNOG 文件常见）
/^#/ { next }

# 读取表头
header == 0 {
    # 兼容可能的 BOM 和 CRLF
    sub(/^\xef\xbb\xbf/, "", $1)
    for (i = 1; i <= NF; i++) {
        sub(/\r$/, "", $i)
        if ($i == "query")        q = i
        if ($i == "COG_category") cog = i
        if ($i == "KEGG_ko")      kegg = i
        if ($i == "GOs")          go = i
        if ($i == "EC")           ec = i
    }

    if (!q || !cog || !kegg || !go || !ec) {
        print "Error: required columns not found!" > "/dev/stderr"
        exit 1
    }

    print "query", "COGs", "KEGG", "GOs", "ECs"
    header = 1
    next
}

# 数据行处理
{
    sub(/\r$/, "", $0)

    c = ($cog  != "-" && $cog  != "") ? 1 : 0
    k = ($kegg != "-" && $kegg != "") ? 1 : 0
    g = ($go   != "-" && $go   != "") ? 1 : 0
    e = ($ec   != "-" && $ec   != "") ? 1 : 0

    print $q, c, k, g, e
}
' "$INPUT" > "$OUTPUT"

echo "Done! Output written to $OUTPUT"