#!/bin/bash

set -e

# ==========================================
# 输入文件
# ==========================================
GTDB="../temp/gtdb_classify/tax.bac120.summary.tsv"
CHECKM="checkm2/quality_report.tsv"
ABUND="coverm/abundance.tsv"

# ==========================================
# 输出目录
# ==========================================
OUTDIR="../temp/gtdb_infer"
mkdir -p ${OUTDIR}

OUTFILE="${OUTDIR}/annotation.txt"

# ==========================================
# 去除Windows换行符
# ==========================================
sed -i 's/\r$//' ${GTDB}
sed -i 's/\r$//' ${CHECKM}
sed -i 's/\r$//' ${ABUND}

# ==========================================
# 写表头
# ==========================================
ABUND_HEADER=$(head -n 1 ${ABUND} | cut -f2-)

echo -e "ID\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tscore\tcompleteness\tcontamination\tstrain_heterogeneity\tsize\tN50\tcluster_members\t${ABUND_HEADER}" \
> ${OUTFILE}

# ==========================================
# 主程序
# ==========================================
awk -F '\t' -v OFS='\t' '

########################################################
# 读取 checkm
########################################################
FNR==NR {

    if(NR==1) next

    id=$1

    completeness=$2
    contamination=$3

    score=completeness-(5*contamination)

    size=$9
    N50=$7

    checkm[id]=score"\t"completeness"\t"contamination"\t0\t"size"\t"N50"\t1"

    next
}

########################################################
# 读取 abundance
########################################################
FILENAME==ARGV[2] {

    if(FNR==1) next

    id=$1

    abund=""

    for(i=2;i<=NF;i++){
        abund=abund"\t"$i
    }

    abundance[id]=substr(abund,2)

    next
}

########################################################
# 读取 GTDB
########################################################
FILENAME==ARGV[3] {

    if(FNR==1) next

    id=$1

    split($2,a,";")

    domain=a[1]
    phylum=a[2]
    class=a[3]
    order=a[4]
    family=a[5]
    genus=a[6]
    species=a[7]

    print \
    id,\
    domain,\
    phylum,\
    class,\
    order,\
    family,\
    genus,\
    species,\
    checkm[id],\
    abundance[id]

}
' \
${CHECKM} \
${ABUND} \
${GTDB} \
>> ${OUTFILE}

echo "======================================"
echo "annotation file generated:"
echo "${OUTFILE}"
echo "======================================"
