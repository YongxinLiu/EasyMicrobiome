#!/bin/bash
set -euo pipefail
##############################################
## Default parameters
##############################################

THREADS=20
BOOTSTRAP=1000
ITERATIONS=1000

##############################################
## Parse parameters
##############################################

usage(){
echo "
Usage:
bash run_FastSpar.sh \
   -c MAG_for_FastSpar_Centenarians_fixed.tsv \
   -y MAG_for_FastSpar_Young_fixed.tsv \
   -o FastSpar_result \
   -t 20 \
   -b 1000 \
   -i 1000
"
exit 1
}

while getopts "c:y:o:t:b:i:h" opt
do
    case $opt in
        c) CENT=$OPTARG;;
        y) YOUNG=$OPTARG;;
        o) OUTDIR=$OPTARG;;
        t) THREADS=$OPTARG;;
        b) BOOTSTRAP=$OPTARG;;
        i) ITERATIONS=$OPTARG;;
        h) usage;;
        *) usage;;
    esac
done

##############################################
## Check parameters
##############################################

[ -z "${CENT:-}" ] && usage
[ -z "${YOUNG:-}" ] && usage
[ -z "${OUTDIR:-}" ] && usage

##############################################
## Check software
##############################################

for exe in fastspar fastspar_bootstrap fastspar_pvalues parallel
do
    command -v ${exe} >/dev/null 2>&1 || {
        echo "ERROR: ${exe} not found!"
        exit 1
    }
done

mkdir -p ${OUTDIR}

##############################################
## Function
##############################################

run_fastspar(){

GROUP=$1
INPUT=$2

echo "========================================"
echo "Processing ${GROUP}"
echo "========================================"

mkdir -p ${OUTDIR}/${GROUP}
mkdir -p ${OUTDIR}/${GROUP}/bootstrap
mkdir -p ${OUTDIR}/${GROUP}/bootstrap_cor
mkdir -p ${OUTDIR}/${GROUP}/bootstrap_cov

##############################################
## Step1 correlation
##############################################

fastspar \
    --otu_table ${INPUT} \
    --correlation ${OUTDIR}/${GROUP}/${GROUP}_correlation.tsv \
    --covariance ${OUTDIR}/${GROUP}/${GROUP}_covariance.tsv \
    --iterations ${ITERATIONS} \
    --threads ${THREADS}

##############################################
## Step2 bootstrap
##############################################

fastspar_bootstrap \
    --otu_table ${INPUT} \
    --number ${BOOTSTRAP} \
    --prefix ${OUTDIR}/${GROUP}/bootstrap/${GROUP}_otu_

##############################################
## Step3 bootstrap correlation
##############################################

parallel -j ${THREADS} \
fastspar \
    --otu_table {} \
    --correlation ${OUTDIR}/${GROUP}/bootstrap_cor/{/.}_cor.tsv \
    --covariance ${OUTDIR}/${GROUP}/bootstrap_cov/{/.}_cov.tsv \
    --iterations 5 \
::: ${OUTDIR}/${GROUP}/bootstrap/*.tsv

##############################################
## Step4 pvalue
##############################################

fastspar_pvalues \
    --otu_table ${INPUT} \
    --correlation ${OUTDIR}/${GROUP}/${GROUP}_correlation.tsv \
    --prefix ${OUTDIR}/${GROUP}/bootstrap_cor/${GROUP}_otu_ \
    --permutations ${BOOTSTRAP} \
    --outfile ${OUTDIR}/${GROUP}/${GROUP}_pvalues.tsv

echo "${GROUP} finished."

}

##############################################
## Run
##############################################

run_fastspar Centenarians ${CENT}

run_fastspar Young ${YOUNG}

echo "====================================="
echo "All FastSpar analyses completed!"
echo "====================================="
