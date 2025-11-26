#!/bin/bash

usage()
{
cat <<EOF
${txtcyn}
Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to summary data.

${txtbld}OPTIONS${txtrst}:
	-f	Data file (with header line, 
		If 'CTCTCT' given to -c, the first column will be treated
		as the rowname and no duplicational names allowed.
		Otherwise every column in files will be treated as values.)
		${bldred}[NECESSARY]${txtrst}
	-k	If the names of your rows and columns startwith numeric value,
		this can be set to FALSE to avoid modifying these names to be
		illegal variable names. But duplicates can not be picked out.
		[${bldred}Default TRUE${txtrst}]
		Accept FALSE.
	-c	Specify the columns you want to do summary.
		[${bldred}Default all colums. This paramater is unused unless
		in the situation mentioned in -f.${txtrst}]
	-m	Get maximum-minimum values.
		[${bldred}Default FALSE, accept TRUE${txtrst}]
	-z	Is there a header[${bldred}Default TRUE${txtrst}]
		Accept FALSE.
	-e	Execute or not[${bldred}Default TRUE${txtrst}]
		Accept FALSE.
EOF
}

file=
checkNames='TRUE'
header='TRUE'
execute='TRUE'
#cols='CTCTCT'
cols='hahaha'
maxmin='FALSE'

while getopts "hf:k:z:e:m:c:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		k)
			checkNames=$OPTARG
			;;
		z)
			header=$OPTARG
			;;
		c)
			cols=$OPTARG
			;;
		m)
			maxmin=$OPTARG
			;;
		e)
			execute=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done
if [ -z $file ]; then
	usage
	exit 1
fi

midname='summary'


cat <<EOF >$file${midname}.r
if ("${cols}" == "CTCTCT"){ 
	data1 = read.table("$file", header=$header,
	sep="\t",row.names=1, comment.char="", check.names=${checkNames})
	x <- as.matrix(data1)
	cv.sum <- summary(data1)
	print(cv.sum)
	if (${maxmin}){
		print(paste0("MAX: ", max(max(data1))))
		print(paste0("MIN: ", min(min(data1))))
	}
	print(summary(colSums(data1)))
	#filesum <- "${file}${midname}"
	#write.table(as.matrix(cv.sum), file=filesum, sep="\t", col.names=F, quote=F)
}else {
	data1 = read.table("$file", header=$header,
	sep="\t",check.names=${checkNames})
	cv.sum <- summary(data1)
	print(cv.sum)
	print(summary(colSums(data1)))
}
EOF

if [ "${execute}" = 'TRUE' ]; then
	Rscript $file${midname}.r
	/bin/rm -f $file${midname}.r
fi
