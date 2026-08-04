
#!/bin/bash
# set -euo pipefail

# ======================
# 参数解析
# ======================
usage() {
    echo "Usage: $0 -i <input_file> -o <output_dir>"
    echo ""
    echo "  -i    MetaPhlAn taxonomy.spf file"
    echo "  -o    Output directory"
    exit 1
}

while getopts "i:o:" opt; do
    case ${opt} in
        i ) input_file=$OPTARG ;;
        o ) outdir=$OPTARG ;;
        * ) usage ;;
    esac
done

# 参数检查
if [[ -z "${input_file:-}" || -z "${outdir:-}" ]]; then
    usage
fi

if [[ ! -f "$input_file" ]]; then
    echo "ERROR: input file not found: $input_file"
    exit 1
fi

mkdir -p "$outdir"

# ======================
# 主逻辑
# ======================
levels=("Kingdom" "Phylum" "Class" "Order" "Family" "Genus" "Species" "Strain")
sample_start=9

for i in {1..8}; do
    level=${levels[$i-1]}
    echo "Processing: $level"

    awk -v col="$i" -v sample_start="$sample_start" 'BEGIN {OFS="\t"}
        NR == 1 {
            printf "Taxonomy"
            for (j = sample_start; j <= NF; j++) {
                printf "\t%s", $j
            }
            print ""
            next
        }
        {
            key = $col
            for (j = sample_start; j <= NF; j++) {
                sums[key, j] += $j
            }
            keys[key] = 1
        }
        END {
            for (key in keys) {
                printf "%s", key
                for (j = sample_start; j <= NF; j++) {
                    printf "\t%.5f", sums[key, j]
                }
                print ""
            }
        }
    ' "$input_file" > "${outdir}/${level}.txt"
done

echo "All levels processed successfully."



