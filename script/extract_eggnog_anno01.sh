#!/bin/bash

# 固定输入输出文件
input="../temp/eggnog/output"
output="eggnog/data_venn1.txt"

awk -F'\t' '
BEGIN{
    OFS="\t"
}

# 跳过注释行（eggNOG文件常见）
/^#/ {next}

# 读取表头
header==0{
    for(i=1;i<=NF;i++){
        if($i=="query") q=i
        if($i=="COG_category") cog=i
        if($i=="KEGG_ko") kegg=i
        if($i=="GOs") go=i
        if($i=="EC") ec=i
    }

    if(!q || !cog || !kegg || !go || !ec){
        print "Error: required columns not found!" > "/dev/stderr"
        exit 1
    }

    print "query","COGs","KEGG","GOs","ECs"
    header=1
    next
}

# 数据行处理
{
    c = ($cog != "-" && $cog != "") ? 1 : 0
    k = ($kegg != "-" && $kegg != "") ? 1 : 0
    g = ($go != "-" && $go != "") ? 1 : 0
    e = ($ec != "-" && $ec != "") ? 1 : 0

    print $q,c,k,g,e
}
' "$input" > "$output"

echo "Done! Output written to $output"