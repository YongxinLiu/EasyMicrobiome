#处理RGI数据--热图
#处理cazy数据--热图（result/MAGdbcan3/cazy.txt）
# 遍历所有匹配的文件
for FILE in ${PREFIX}*.txt; do
    # 获取文件名的基本部分
    BASENAME=$(basename "$FILE")
    OUTPUT="ARG/${BASENAME}"

    # 提取第1列，第17列和第28列，并写入对应的输出文件
    awk -F'\t' -v OFS='\t' 'BEGIN {print "ORF_ID"} NR>1 {print $1}' "$FILE" > "$OUTPUT"
done

# 进入输出目录
cd ARG || exit

# 处理每个输出文件
for file in *.txt; do
    # 去掉回车符并更新文件
    awk -v fname="${file%.txt}" 'BEGIN {FS=OFS="\t"} {gsub(/\r/, ""); if (NR==1) print $0, fname; else print $1, $2, 1}' "$file" > "${file%.txt}_updated.txt"

    echo "$file 处理完成"
done
mkdir -p updated
mv *updated.txt updated
cd updated
# merge2table
conda activate humann3
humann_join_tables \
  --input ./ --file_name Mx_All \
  --output ./ARG.txt
sed -i 's/^[ \t]*$/1/g' ARG.txt


mkdir ARG2
# 获取ORF_ID对应的基因和抗生素注释
for FILE in ${PREFIX}*.txt; do
    # 获取文件名的基本部分
    BASENAME=$(basename "$FILE")
    OUTPUT="ARG2/${BASENAME}"
    # 提取第17列和第28列，并写入对应的输出文件
    awk -F'\t' -v OFS='\t' 'BEGIN {print "ORF_ID", "AMR_Gene_Family", "Antibiotic"} NR>1 {print $1, $17, $28}' "$FILE" > "$OUTPUT"
done


for FILE in ${PREFIX}*.txt; do
    # 获取文件名的基本部分
    BASENAME=$(basename "$FILE")
    OUTPUT="ARG2/${BASENAME}"
    # 提取第17列和第28列，并写入对应的输出文件
    awk -F'\t' -v OFS='\t' 'NR>1 {print $1, $17, $28}' "$FILE" > "$OUTPUT"
done
cd ARG2
cat ARG2/*.txt>ARG2/annotation.txt
echo -e "ORF_ID\tAMR_Gene_Family\tAntibiotic\n$(cat ARG2/annotation.txt)" > ARG2/annotation.txt
#排序ORF_ID列，手动合并ARG2/annotation.txt和ARG/updated/ARG.txt,并且将空白处替换为1，保存为ARG_F.txt
