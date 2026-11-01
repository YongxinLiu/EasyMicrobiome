#!/bin/bash
# Copyright 2026 Defeng Bai <baidefeng@caas.cn>

set -euo pipefail

###############################################################
# Usage
#
# bash run_sparcc3.sh input.txt R.txt P.txt
#
###############################################################

if [ $# -ne 3 ]; then
    echo "Usage:"
    echo "bash $0 input_abundance.txt correlation.txt pvalue.txt"
    exit 1
fi

#############################
# Input
#############################

INPUT=$(realpath "$1")
R_OUT=$(realpath "$2")
P_OUT=$(realpath "$3")

#############################
# Parameters
#############################

ITER=20
BOOTSTRAP=100

#############################
# Working directory
#############################

WD=$(pwd)

#############################
# Download SparCC3
#############################

if [ ! -d "${WD}/SparCC3" ]; then
    echo "Downloading SparCC3..."
    git clone https://github.com/JCSzamosi/SparCC3.git
else
    echo "SparCC3 already exists."
fi

cd "${WD}/SparCC3"

mkdir -p data
mkdir -p example/basis_corr
mkdir -p example/pvals

#############################
# Copy input
#############################

cp "${INPUT}" data/input.txt

#############################
# Step1 Correlation
#############################

python SparCC.py \
    data/input.txt \
    -i ${ITER} \
    --cor_file=example/basis_corr/correlation.tsv \
    > example/basis_corr/sparcc.log

#############################
# Step2 Bootstrap
#############################

python MakeBootstraps.py \
    data/input.txt \
    -n ${BOOTSTRAP} \
    -t bootstrap_#.txt \
    -p example/pvals/ \
    >> example/basis_corr/sparcc.log

#############################
# Step3 Bootstrap correlations
#############################

for ((n=0;n<BOOTSTRAP;n++))
do
    python SparCC.py \
        example/pvals/bootstrap_${n}.txt \
        -i ${ITER} \
        --cor_file=example/pvals/bootstrap_cor_${n}.txt \
        >> example/basis_corr/sparcc.log
done

#############################
# Step4 P values
#############################

python PseudoPvals.py \
    example/basis_corr/correlation.tsv \
    example/pvals/bootstrap_cor_#.txt \
    ${BOOTSTRAP} \
    -o example/pvals/pvalue.tsv \
    -t two_sided \
    >> example/basis_corr/sparcc.log

#############################
# Output
#############################

cp example/basis_corr/correlation.tsv "${R_OUT}"
cp example/pvals/pvalue.tsv "${P_OUT}"

echo
echo "Finished!"
echo "Correlation : ${R_OUT}"
echo "P-value     : ${P_OUT}"
