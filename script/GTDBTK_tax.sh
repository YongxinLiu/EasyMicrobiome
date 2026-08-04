# #!/bin/bash
# 
# awk -F'\t' '
# BEGIN{
#     OFS="\t"
#     print "user_genome","classification","ID"
# }
# NR>1{
#     species=$8
#     sub(/^s__/,"",species)
# 
#     print $1,species,species" "$1
# }
# ' itol/tax01.txt > itol/GTDBTK_npc.txt


#!/bin/bash
# ============================================
# Convert GTDB-Tk taxonomy to iTOL annotation
# Usage:
#   bash gtdbtk2itol.sh -i input.txt -o output.txt
# ============================================

usage() {
    echo "Usage: $0 -i <input_file> -o <output_file>"
    exit 1
}

while getopts "i:o:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

[ -z "$INPUT" ] && usage
[ -z "$OUTPUT" ] && usage

if [ ! -f "$INPUT" ]; then
    echo "Error: Input file '$INPUT' not found."
    exit 1
fi

awk -F'\t' '
BEGIN{
    OFS="\t"
    print "user_genome","classification","ID"
}
NR>1{
    species=$8
    sub(/^s__/,"",species)

    print $1,species,species" "$1
}
' "$INPUT" > "$OUTPUT"

echo "Done!"
echo "Input : $INPUT"
echo "Output: $OUTPUT"