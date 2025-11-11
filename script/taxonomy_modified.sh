#!/bin/bash
set -e

# Default parameter 默认参数
input=result/metaphlan4/taxonomy.spf
output=result/metaphlan4/taxonomy_modified.spf

# Function for script description and usage 脚本功能描述
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    taxonomy_modified.sh
Version:     1.0
Date:        2025/11/11
Author:      Defeng Bai(白德凤), Yong-Xin Liu(刘永鑫)
Email:       liuyongxin@caas.cn
Website:     https://github.com/YongxinLiu/EasyMicrobiome
Description: Format metaphlan4 result for GraPhlAn
Notes:
  -------------------------------------------------------------------------------
Copyright:   2016-2025 (c) Yong-Xin Liu
License:     GPL
If used this script, please cited:
Bai, Defeng, Tong Chen, Jiani Xun, Chuang Ma, Hao Luo, Haifei Yang, Chen Cao, et al. 2025. 
“EasyMetagenome: A user-friendly and flexible pipeline for shotgun metagenomic analysis in microbiome research.” 
iMeta 4: e70001. https://doi.org/10.1002/imt2.70001
-------------------------------------------------------------------------------
Version 1.0 2025/11/11
Add command line parameter for input and output

# Input files: result/metaphlan4/taxonomy.spf

# 1. MetaPhlan4 taxonomic table 物种组成表
Kingdom Phylum  Class   Order   Family  Genus   Species Strain  C1      C2      C3      Y1      Y2      Y3
k__Bacteria     p__Verrucomicrobia      c__Verrucomicrobiae     o__Verrucomicrobiales   f__Akkermansiaceae      g__Akkermansia  s__Akkermansia_muciniphila      t__SGB922>
k__Bacteria     p__Verrucomicrobia      c__Verrucomicrobiae     o__Verrucomicrobiales   f__Akkermansiaceae      g__Akkermansia  s__Akkermansia_muciniphila      t__SGB922>
k__Bacteria     p__Synergistetes        c__Synergistia  o__Synergistales        f__Synergistaceae       g__Pyramidobacter       s__Pyramidobacter_piscolens     t__SGB154>

# Output file: result/metaphlan4/taxonomy_modified.spf

# 1. Formated MetaPhlan4 taxonomic table 修改的物种组成表
Kingdom Phylum  Class   Order   Family  Genus   Species C1      C2      C3      Y1      Y2      Y3
k__Bacteria     p__Bacteroidetes        c__CFGB529      o__OFGB529      f__FGB529       g__GGB1096      s__GGB1096_SGB1408      0       0       1       0       0       0
k__Bacteria     p__Firmicutes   c__Clostridia   o__Eubacteriales        f__Eubacteriales_Family_XIII_Incertae_Sedis     g__Mogibacterium        s__Mogibacterium_sp_BX12 >
k__Bacteria     p__Firmicutes   c__Clostridia   o__Eubacteriales        f__Oscillospiraceae     g__Dysosmobacter        s__Dysosmobacter_welbionis      0       1       0>
k__Bacteria     p__Proteobacteria       c__Gammaproteobacteria  o__Enterobacterales     f__Enterobacteriaceae   g__Kluyvera     s__Kluyvera_genomosp_3  0       0       1>

OPTIONS:
-h help
-i input, default result/metaphlan4/taxonomy.spf
-o output , default result/metaphlan4/taxonomy_modified.spf
-? show help of script

Example:
  taxonomy_modified.sh -i result/metaphlan4/taxonomy.spf -o result/metaphlan4/taxonomy_modified.spf

EOF
}


# Analysis parameter 参数解析 
while getopts "i:o:" OPTION
do
	case $OPTION in
		i)
			input=$OPTARG
			;;
		o)
			output=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done



# 处理 taxonomy2.spf 文件，生成 taxonomy_merged.spf
cut -f 1-7,9- $input | awk -F '\t' '
NR==1 {
    print;  # 打印表头
    header = $0;  # 保存表头
    next;  # 跳过表头
}
{
    taxonomy = $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6; 
    species = $7;  # 获取物种
    for (i=8; i<=NF; i++) {
        # 统计出现次数，保持0，非0改为1
        if ($i > 0) {
            sum[species, i] += 1;  # 不为0则计为1
        } else {
            sum[species, i] += 0;  # 为0保持0
        }
    }
    species_seen[species] = taxonomy; 
}
END {
    for (sp in species_seen) {
        printf "%s\t%s", species_seen[sp], sp;
        for (i=8; i<=NF; i++) {
            printf "\t%d", sum[sp, i]; 
        }
        print "";
    }
}' > tmp

# 读取已生成的 taxonomy_merged.spf 文件并进行修改，输出到 taxonomy_modified.spf
awk -F '\t' '
NR==1 {
    print;  # 打印表头
    header = $0;  # 保存表头
    next;  # 跳过表头
}
{
    # 打印前7列
    printf "%s\t%s", $1, $2;  # 打印前两列
    printf "\t%s\t%s\t%s\t%s\t%s", $3, $4, $5, $6, $7;  # 打印后五列

    # 从第8列开始进行处理
    for (i=8; i<=NF; i++) {
        if ($i > 0) {
            printf "\t1";  # 非0的值改为1
        } else {
            printf "\t0";  # 0保持为0
        }
    }
    print "";  # 换行
}' tmp > $output


rm -rf tmp
