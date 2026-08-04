#!/bin/bash



# 检查参数
if [ $# -eq 0 ]; then
    echo "Usage:"
    echo "bash $0 \"Pathway1\" \"Pathway2\" ..."
    exit 1
fi

# 输入文件
ANNOTATION="humann4/KEGG_KO_Annotation.tsv"
KO_FILE="humann4/ko_stratified.tsv"

# 输出文件
SELECTED="humann4/KEGG_Selected_Pathways.tsv"
FILTERED="humann4/KEGG_Selected_Pathways_filtered.tsv"
SANKEY="humann4/data_sankey2.txt"

# 构建正则表达式
PATTERN=$(printf "%s|" "$@")
PATTERN=${PATTERN%|}

echo "Selected pathways:"
printf "  %s\n" "$@"

############################
# Step 1. 提取目标通路
############################
awk -F'\t' -v pattern="$PATTERN" '
NR==1 || $9 ~ ("^(" pattern ")$")
' "$ANNOTATION" > "$SELECTED"

############################
# Step 2. 保留样本中存在的KO
############################
awk -F'\t' '
NR==FNR{
    if(FNR>1){
        split($1,a,"|")
        ko[a[1]]=1
    }
    next
}
FNR==1 || ($1 in ko)
' "$KO_FILE" "$SELECTED" > "$FILTERED"

############################
# Step 3. 生成 Sankey 输入文件
############################
awk -F'\t' '
BEGIN{
    OFS="\t"
    print "Level2_Name","Level3_Name","KO"
}
NR>1{
    print $7,$9,$1
}
' "$FILTERED" > "$SANKEY"

echo "Done!"
echo "Output files:"
echo "  $SELECTED"
echo "  $FILTERED"
echo "  $SANKEY"
