#!/bin/bash

OUTDIR="humann4"

mkdir -p ${OUTDIR}

echo "Downloading KEGG hierarchy..."

wget -q -O ${OUTDIR}/ko00001.keg \
"https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001.keg&format=htext"

echo "Parsing hierarchy..."

awk '
BEGIN{
    OFS="\t"

    print "KO",
          "Gene",
          "Description",
          "Level1_ID",
          "Level1_Name",
          "Level2_ID",
          "Level2_Name",
          "Level3_ID",
          "Level3_Name"
}

/^A/{
    line=$0
    sub(/^A/,"",line)

    match(line,/^[0-9]+/)
    L1ID=substr(line,RSTART,RLENGTH)

    sub(/^[0-9]+ /,"",line)
    L1NAME=line
}

/^B/{
    line=$0
    sub(/^B[ ]+/,"",line)

    match(line,/^[0-9]+/)
    L2ID=substr(line,RSTART,RLENGTH)

    sub(/^[0-9]+ /,"",line)
    L2NAME=line
}

/^C/{
    line=$0
    sub(/^C[ ]+/,"",line)

    match(line,/^[0-9]+/)
    L3ID=substr(line,RSTART,RLENGTH)

    sub(/^[0-9]+ /,"",line)
    sub(/ \[PATH:ko[0-9]+\]/,"",line)

    L3NAME=line
}

/^D/{
    line=$0
    sub(/^D[ ]+/,"",line)

    split(line,a," ")

    KO=a[1]

    sub(/^K[0-9]+[ ]+/,"",line)

    split(line,b,";")

    Gene=b[1]

    Desc=b[2]
    gsub(/^[ ]+/,"",Desc)

    print KO,
          Gene,
          Desc,
          L1ID,
          L1NAME,
          L2ID,
          L2NAME,
          L3ID,
          L3NAME
}
' ${OUTDIR}/ko00001.keg > ${OUTDIR}/KEGG_KO_Annotation.tsv

echo ""
echo "Finished."
echo ""
echo "Generated files:"
echo "${OUTDIR}/ko00001.keg"
echo "${OUTDIR}/KEGG_KO_Annotation.tsv"
