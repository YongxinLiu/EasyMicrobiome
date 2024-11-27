#!/bin/bash

# 创建输出文件并添加表头
echo -e "Name\tCompleteness\tContamination\tGenome_Size\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies" > result/checkm2/parsed_taxonomy_filled.tsv

# 处理数据并生成最终结果
join -t $'\t' -1 1 -2 1 \
    <(tail -n +2 result/checkm2/quality_report.tsv | cut -f1-3,9 | sort -k1,1) \
    <(tail -n +2 temp/gtdb_classify/tax.bac120.summary.tsv | cut -f1-2 | sort -k1,1) | \
awk -F '\t' -v OFS="\t" '
{
    # 默认分类为 unclassified
    domain="unclassified"
    phylum="unclassified"
    class="unclassified"
    order="unclassified"
    family="unclassified"
    genus="unclassified"
    species="unclassified"

    # 将分类信息分割到数组中
    split($5, tax, ";")
    for (i in tax) {
        if (tax[i] ~ /^d__/) domain=substr(tax[i], 4)
        if (tax[i] ~ /^p__/) phylum=substr(tax[i], 4)
        if (tax[i] ~ /^c__/) class=substr(tax[i], 4)
        if (tax[i] ~ /^o__/) order=substr(tax[i], 4)
        if (tax[i] ~ /^f__/) family=substr(tax[i], 4)
        if (tax[i] ~ /^g__/) genus=substr(tax[i], 4)
        if (tax[i] ~ /^s__/) species=substr(tax[i], 4)
    }
    # 打印最终结果
    print $1, $2, $3, $4, domain, phylum, class, order, family, genus, species
}' >> result/checkm2/parsed_taxonomy_filled.tsv

# 将空白字段替换为 "unclassified"
awk -v OFS="\t" '{for(i=5;i<=11;i++) if($i=="") $i="unclassified"; print}' result/checkm2/parsed_taxonomy_filled.tsv > result/checkm2/taxonomy_merge.txt

rm result/checkm2/parsed_taxonomy_filled.tsv
# result/checkm2/taxonomy_merge.csv

