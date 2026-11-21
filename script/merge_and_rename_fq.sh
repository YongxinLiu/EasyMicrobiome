#!/bin/bash
set -euo pipefail

# ======================
# Help document
# ======================
usage() {
cat <<EOF

merge_and_rename_fq.sh
=====================

Description:
  Link paired-end fastq files into a single directory and rename them
  according to a metadata table.

  This script creates symbolic links (NOT copies).
  Original raw data will NOT be modified.

Required input:
  - Metadata file with at least TWO columns:
      column 1: NewSampleID
      column 2: OldSampleID (fastq prefix OR directory name)

Fastq naming support:
  Mode "name" (default):
    <OldSampleID>_1.fq.gz
    <OldSampleID>_2.fq.gz

  Mode "dir":
    raw_dir/OldSampleID/*.fq.gz
    fastq names can be arbitrary but must contain R1 / R2 or 1 / 2

Usage:
  bash merge_and_rename_fq.sh -m metadata.txt -i raw_fq_dir -o output_dir [-M name|dir]

Options:
  -m   Metadata file (tab- or space-delimited, first line is header)
  -i   Directory containing raw fastq files (recursive search)
  -o   Output directory for renamed fastq symbolic links
  -M   Match mode:
         name (default) : match fastq file names
         dir            : match by directory name
  -h   Show this help message and exit

Output:
  output_dir/
    ├── SampleID_R1.fq.gz
    └── SampleID_R2.fq.gz

EOF
exit 0
}

# ======================
# Parse arguments
# ======================
mode="name"

while getopts "m:i:o:M:h" opt; do
  case $opt in
    m) metadata=$OPTARG ;;
    i) raw_dir=$OPTARG ;;
    o) out_dir=$OPTARG ;;
    M) mode=$OPTARG ;;
    h) usage ;;
    *) usage ;;
  esac
done

# ======================
# Check arguments
# ======================
[[ -z "${metadata:-}" || -z "${raw_dir:-}" || -z "${out_dir:-}" ]] && usage

mkdir -p "$out_dir"

echo "[INFO] Metadata : $metadata"
echo "[INFO] Raw dir  : $raw_dir"
echo "[INFO] Out dir  : $out_dir"
echo "[INFO] Mode     : $mode"
echo "[INFO] Start processing..."
echo

# ======================
# Main loop
# ======================
tail -n +2 "$metadata" | while read -r new old; do
    # 去掉 Windows 回车
    new=$(echo "$new" | tr -d '\r')
    old=$(echo "$old" | tr -d '\r')

    [[ -z "$new" || -z "$old" ]] && continue

    echo "[INFO] Processing: $old → $new"

    R1=""
    R2=""

    if [[ "$mode" == "name" ]]; then
        # ========= 按 fastq 文件名匹配 =========
        R1=$(find "$raw_dir" -type f -name "${old}_1.fq.gz" | head -n 1)
        R2=$(find "$raw_dir" -type f -name "${old}_2.fq.gz" | head -n 1)

    elif [[ "$mode" == "dir" ]]; then
        # ========= 按目录名匹配 =========
        sample_dir=$(find "$raw_dir" -type d -name "$old" | head -n 1)

        if [[ -z "$sample_dir" ]]; then
            echo "[WARNING] Sample directory not found: $old"
            continue
        fi

        R1=$(find "$sample_dir" -type f \( \
              -name "*_R1_*.fq.gz" -o \
              -name "*_R1.fq.gz"  -o \
              -name "*_1.fq.gz" \
            \) | head -n 1)

        R2=$(find "$sample_dir" -type f \( \
              -name "*_R2_*.fq.gz" -o \
              -name "*_R2.fq.gz"  -o \
              -name "*_2.fq.gz" \
            \) | head -n 1)
    else
        echo "[ERROR] Unknown mode: $mode"
        exit 1
    fi

    if [[ -z "$R1" || -z "$R2" ]]; then
        echo "[WARNING] Missing fastq files for $old, skip"
        continue
    fi

    ln -sf "$(realpath "$R1")" "${out_dir}/${new}_R1.fq.gz"
    ln -sf "$(realpath "$R2")" "${out_dir}/${new}_R2.fq.gz"

    echo "[OK] ${new}_R1.fq.gz  ${new}_R2.fq.gz"
done

echo
echo "[DONE] All samples processed successfully."
