#!/bin/bash

levels=("Kingdom" "Phylum" "Class" "Order" "Family" "Genus" "Species" "Strain")
input_file="metaphlan4/taxonomy.spf"
header=$(head -n 1 "$input_file")
sample_start=9

for i in {1..8}; do
    echo "Processing: ${levels[$i-1]}"
    awk -v col=$i -v sample_start="$sample_start" -v header="$header" 'BEGIN {OFS="\t"}
    NR == 1 { 
        printf "Taxonomy\t"
        for (j = sample_start; j <= NF; j++) {
            printf "%s%s", (j == sample_start ? "" : OFS), $j
        }
        print "";
        next;
    }
    {
        key = $col  
        for (j = sample_start; j <= NF; j++) { 
            sums[key, j] += $j
        }
        keys[key] = 1;
    } 
    END {
        for (key in keys) {
            printf "%s", key
            for (j = sample_start; j <= NF; j++) {
                printf "\t%.5f", sums[key, j]
            }
            print "" 
        }
    }' "$input_file" > "metaphlan4/${levels[$i-1]}.txt" 
done
