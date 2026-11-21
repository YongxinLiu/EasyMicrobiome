#!/bin/bash
set -e

# ======================
# Default parameters 默认参数
# ======================
input="result/metaphlan4/taxonomy.spf"
output="result/metaphlan4/taxonomy_modified.spf"

# ======================
# Help document
# ======================
usage() {
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    taxonomy_modified.sh
Version:     1.0
Date:        2025/11/11
Author:      Defeng Bai(白德凤), Yong-Xin Liu(刘永鑫)
Email:       liuyongxin@caas.cn
Website:     https://github.com/YongxinLiu/EasyMicrobiome
Description:
  Format MetaPhlAn4 taxonomic table for GraPhlAn visualization.

Notes:
  - Convert MetaPhlAn4 species abundance to presence/absence (0/1)
  - Merge duplicated species entries
  - Keep taxonomy structure from Kingdom to Species

License:     GPL (2016-2025)

If you use this script, please cite:
Bai, Defeng, et al. 2025.
EasyMetagenome: A user-friendly and flexible pipeline for shotgun metagenomic analysis.
iMeta 4: e70001.

-------------------------------------------------------------------------------

Options:
  -i   Input MetaPhlAn4 taxonomy file
       default: result/metaphlan4/taxonomy.spf

  -o   Output formatted taxonomy file
       default: result/metaphlan4/taxonomy_modified.spf

  -h   Show this help message and exit
  -?   Same as -h

Example:
  bash taxonomy_modified.sh \\
    -i result/metaphlan4/taxonomy.spf \\
    -o result/metaphlan4/taxonomy_modified.spf

-------------------------------------------------------------------------------
EOF
exit 0
}

# ======================
# Parse arguments 参数解析
# ======================
while getopts ":i:o:h?" OPTION; do
  case $OPTION in
    i)
      input=$OPTARG
      ;;
    o)
      output=$OPTARG
      ;;
    h|\?)
      usage
      ;;
    :)
      echo "[ERROR] Option -$OPTARG requires an argument." >&2
      usage
      ;;
    *)
      usage
      ;;
  esac
done

# ======================
# Check input
# ======================
if [[ ! -f "$input" ]]; then
  echo "[ERROR] Input file not found: $input" >&2
  exit 1
fi

echo "[INFO] Input  file : $input"
echo "[INFO] Output file : $output"
echo "[INFO] Start processing..."
echo

# ======================
# Step 1: merge species & binarize abundance
# ======================
cut -f 1-7,9- "$input" | awk -F '\t' '
NR==1 {
    print;
    next;
}
{
    taxonomy = $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6;
    species = $7;
    for (i=8; i<=NF; i++) {
        if ($i > 0) sum[species, i] += 1;
    }
    species_seen[species] = taxonomy;
}
END {
    for (sp in species_seen) {
        printf "%s\t%s", species_seen[sp], sp;
        for (i=8; i<=NF; i++) {
            printf "\t%d", (sum[sp, i] > 0 ? 1 : 0);
        }
        print "";
    }
}' > tmp.taxonomy_merged.spf

# ======================
# Step 2: output formatted table
# ======================
awk -F '\t' '
NR==1 {
    print;
    next;
}
{
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s", $1,$2,$3,$4,$5,$6,$7;
    for (i=8; i<=NF; i++) {
        printf "\t%d", ($i > 0 ? 1 : 0);
    }
    print "";
}' tmp.taxonomy_merged.spf > "$output"

rm -f tmp.taxonomy_merged.spf

echo
echo "[DONE] taxonomy_modified.spf generated successfully."
