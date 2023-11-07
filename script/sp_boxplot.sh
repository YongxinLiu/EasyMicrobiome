#!/bin/bash

#set -x

usage()
{
cat <<EOF
${txtcyn}

***CREATED BY Chen Tong (chentong_biology@163.com)***

Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to do boxplot using ggplot2.

fileformat for -f (suitable for data extracted from one sample, the
number of columns is unlimited. Column 'Set' is not necessary unless
you have multiple groups)

Matrix1

Name	2cell_1	2cell_2	2cell_3	2cell_4	2cell_5	2cell_6	4cell_1	4cell_2	4cell_3	4cell_4	4cell_5	4cell_6	zygote_1	zygote_2	zygote_3	zygote_4	zygote_5	zygote_6
A	8	13	14	9	19	12	3.2	5.2	5.6	3.6	7.6	4.8	0.8	1.3	1.4	0.9	1.9	1.2
B	8	13	14	9	19	12	3.2	5.2	5.6	3.6	7.6	4.8	0.8	1.3	1.4	0.9	1.9	1.2
C	8	13	14	9	19	12	3.2	5.2	5.6	3.6	7.6	4.8	0.8	1.3	1.4	0.9	1.9	1.2
D	8	13	14	9	19	12	3.2	5.2	5.6	3.6	7.6	4.8	0.8	1.3	1.4	0.9	1.9	1.2
E	8	13	14	9	19	12	3.2	5.2	5.6	3.6	7.6	4.8	0.8	1.3	1.4	0.9	1.9	1.2
F	8	13	14	9	19	12	3.2	5.2	5.6	3.6	7.6	4.8	0.8	1.3	1.4	0.9	1.9	1.2
G	8	13	14	9	19	12	3.2	5.2	5.6	3.6	7.6	4.8	0.8	1.3	1.4	0.9	1.9	1.2
H	8	13	14	9	19	12	3.2	5.2	5.6	3.6	7.6	4.8	0.8	1.3	1.4	0.9	1.9	1.2
I	8	13	14	9	19	12	3.2	5.2	5.6	3.6	7.6	4.8	0.8	1.3	1.4	0.9	1.9	1.2

Matrix2

Gene	hmC	expr	Set
NM_001003918_26622	0	83.1269257376101	TP16
NM_001011535_3260	0	0	TP16
NM_001012640_14264	0	0	TP16
NM_001012640_30427	0	0	TP16
NM_001003918_2662217393_30486	0	0	TP16
NM_001017393_30504	0	0	TP16
NM_001025241_30464	0	0	TP16
NM_001017393_30504001025241_30513	0	0	TP16

sampleGroupFile 
# (TAB separated, first column corresponds to first row of matrix)
# Group should be gave to <-F>
Sample	Group
zygote_1	zygote
zygote_2	zygote
zygote_3	zygote
zygote_4	zygote
zygote_5	zygote
zygote_6	zygote
2cell_1	2cell
2cell_2	2cell
2cell_3	2cell
2cell_4	2cell
2cell_5	2cell
2cell_6	2cell
4cell_1	4cell
4cell_2	4cell
4cell_3	4cell
4cell_4	4cell
4cell_5	4cell
4cell_6	4cell

For file using "Set" column, you can use 
boxplot.onefile.sh -f file -a Set 

fileformat when -m is true
#Default we use string "value" and "variable" to represent the data
#column and sub-class column. If you have other strings as column
#names, please give them to -d and -F.
#The "Set" column is optional.
#If you do have several groups, they can put at the "Set" column 
#with "Set" or other string as labels. The label should be given
#to parameter -a.
#Actually this format is the melted result of last format.

Matrix_melted

Gene	Sample	Group	Expr
A	zygote_1	zygote	0.8
A	zygote_2	zygote	1.3
A	zygote_3	zygote	1.4
A	zygote_4	zygote	0.9
A	zygote_5	zygote	1.9
A	zygote_6	zygote	1.2
A	2cell_1	2cell	8
A	2cell_2	2cell	13
A	2cell_3	2cell	14
A	2cell_4	2cell	9
A	2cell_5	2cell	19
A	2cell_6	2cell	12
A	4cell_1	4cell	3.2
A	4cell_2	4cell	5.2
A	4cell_3	4cell	5.6
A	4cell_4	4cell	3.6
A	4cell_5	4cell	7.6
A	4cell_6	4cell	4.8

${txtbld}OPTIONS${txtrst}:
	-f	Data file (with header line, the first row is the
 		colname, tab seperated. Multiple formats are allowed and described above)
		${bldred}[NECESSARY]${txtrst}
	-m	When true, it will skip preprocess. But the format must be
		the same as listed before (Matrix_melted).
		${bldred}[Default FALSE, accept TRUE]${txtrst}
	-d	The column represents the digital values, such as Expr.
		${bldred}[Default "value" represents the column named "value".
		This parameter can only be set when -m is TRUE.]${txtrst}
	-F	The column represents the variable information, meaning legend_variable.
		If no-subclass of X-variavle, this will be treated as X-axis variable.
		${bldred}[Default "variable" represents the column named "variable". 
		This parameter can only be set when -m is TRUE.]${txtrst}
	-I	Other columns you want to treat as ID variable columns except
		the one given to -a. Not used when <-m TRUE>.
		${bldred}[Default empty string, accept comma separated strings
		like "'Id1','Id2','Id3'" or single string "id1"]${txtrst}
	-a	Name for x-axis variable
		[${txtred}Default variable, which is an inner name, suitable 
		for data without 'Set' column. For the given example, 
		'Group' which represents groups of each gene should be 
		supplied to this parameter.
		This parameter can only be set when -m is TRUE.
		${txtrst}]
	-b	Rotation angle for x-axis value(anti clockwise)
		${bldred}[Default 0]${txtrst}
	-R	Rotate the plot from vertical to horizontal. 
		Usefull for plots with many values or very long labels at X-axis.
		${bldred}[Default FALSE]${txtrst}
	-l	Levels for legend variable
		[${txtred}Default data order,accept a string like
		"'TP16','TP22','TP23'" ***for <variable> column***.
	   	${txtrst}]
	-q	Giving one gene ID to do boxplot specifically for this gene.
		${bldred}[Default FALSE, accept a string]${txtrst}
	-Q	Giving a sampleGroup file with format specified above to
   		tell the group information for each sample.	
		When <-Q> is given, <-F> and <-a> should be one of the column 
		names of sampleGrp file.
		${bldred}[Default FALSE, accept a file name]${txtrst}
	-D	Self-define intervals for legend variable when legend is
		continuous numbers. Accept either a
		numeric vector of two or more cut points or a single number
		(greater than or equal to 2) giving the number of intervals
		into what 'x' is to be cut. This has higher priority than -l.
		[10 will generate 10 intervals or 
		"c(-1, 0, 1, 2, 5, 10)" will generate (-1,0],(0,1]...(5,10]]	
	-P	[Uppercase P] Legend position[${txtred}Default right. Accept
		top,bottom,left,none, or c(0.08,0.8) (relative to left-bottom).${txtrst}]
	-L	Levels for x-axis variable
		[${txtred}Default data order,accept a string like
		"'g','a','j','x','s','c','o','u'" ***for <Set> column***.
	   	${txtrst}]
	-B	Self-define intervals for x-axis variable. Accept either a
		numeric vector of two or more cut points or a single number
		(greater than or equal to 2) giving the number of intervals
		into what 'x' is to be cut. This has higher priority than -L.
		[10 will generate 10 intervals or 
		"c(-1, 0, 1, 2, 5, 10)" will generate (-1,0],(0,1]...(5,10]]	
	-n	Using notch (sand clock shape) or not.${txtred}[Default FALSE]${txtrst}
	-V	Do violin plot instead of boxplot.${txtred}[Default FALSE]${txtrst}
	-W	Do violin plot without inner boxplot.${txtred}[Default FALSE]${txtrst}
	-j	Do jitter plot instead of boxplot.${txtred}[Default FALSE]${txtrst}
	-J	Do jitter plot overlay with violinplot or boxplot or both.${txtred}[Default FALSE]${txtrst}
	-A	The value given to scale for violin plot.
		if "area", all violins have the same area (before trimming the tails). 
		If "count", areas are scaled proportionally to the number of observations. 
		If "width", all violins have the same maximum width. 
		'equal' is also accepted.
		${txtred}[Default 'width']${txtrst}
	-G	Wrap plots by given column. This is used to put multiple plot
		in one picture. Used when -m is TRUE, normally a string <set>
		should be suitable for this parameter.
	-g	The levels of wrapping to set the order of each group.
		${txtred}Normally the unique value of the column given to B in
		a format like <"'a','b','c','d'">.${txtrst}
	-M	The number of rows one wants when -G is used.Default NULL.
		${txtred}[one of -M and -N is enough]${txtrst}
	-N	The number of columns one wants when -G is used.Default NULL.
		${txtred}[one of -M and -N is enough]${txtrst}
	-k	Paramter for scales for facet.
		[${txtred}Optional, only used when -B is given. Default each 
		inner graph use same scale [x,y range]. 
		'free' (variable x, y ranges for each sub-plot),
		'free_x' (variable x ranges for each sub-plot),'free_y' 
		is accepted. ${txtrst}]
	-t	Title of picture[${txtred}Default empty title${txtrst}]
	-x	xlab of picture[${txtred}Default empty xlab${txtrst}]
	-y	ylab of picture[${txtred}Default empty ylab${txtrst}]
	-s	Scale y axis
		[${txtred}Default null. Accept TRUE.
		Also if the supplied number after -S is not 0, this
		parameter will be set to TRUE${txtrst}]
	-v	If scale is TRUE, give the following
		scale_y_log10()[default], coord_trans(y="log10"), 
		scale_y_continuous(trans=log2_trans()), coord_trans(y="log2"), 
	   	or other legal command for ggplot2)${txtrst}]
	-o	Exclude outliers.
		[${txtred}Exclude outliers or not, default FALSE means not.${txtrst}]
	-O	The scales for you want to zoom in to exclude outliers.
		[${txtred}Default 1.05. No recommend to change unless you know
		what you are doing.${txtrst}]
	-S	A number to add if scale is used.
		[${txtred}Default 0. If a non-zero number is given, -s is
		TRUE.${txtrst}]	
	-c	Manually set colors for each box.[${txtred}Default FALSE,
		meaning using ggplot2 default.${txtrst}]
	-C	Color for each box.[${txtred}When -c is TRUE, str in given
		format must be supplied, ususlly the number of colors should
		be equal to the number of lines.
		"'red','pink','blue','cyan','green','yellow'" or
		"rgb(255/255,0/255,0/255),rgb(255/255,0/255,255/255),rgb(0/255,0/255,255/255),
		rgb(0/255,255/255,255/255),rgb(0/255,255/255,0/255),rgb(255/255,255/255,0/255)"
		${txtrst}]
	-p	[Lowercase p] Other legal R codes for gggplot2 will be given here.
		[${txtres}Begin with '+' ${txtrst}]
	-w	The width of output picture.[${txtred}Default 20${txtrst}]
	-u	The height of output picture.[${txtred}Default 12${txtrst}] 
	-r	The resolution of output picture.[${txtred}Default 300 ppi${txtrst}]
	-E	The type of output figures.[${txtred}Default pdf, accept
		eps/ps, tex (pictex), pdf, jpeg, tiff, bmp, svg and wmf)${txtrst}]
	-z	Is there a header[${bldred}Default TRUE${txtrst}]
	-e	Execute or not[${bldred}Default TRUE${txtrst}]
	-i	Install depeneded packages[${bldred}Default FALSE${txtrst}]
EOF
}

file=
title=''
melted='FALSE'
xlab=' '
ylab=' '
xvariable=''
value='value'
variable='variable'
xtics_angle=0
xtics='TRUE'
level=""
legend_cut=""
x_level=""
x_cut=""
scaleY='FALSE'
y_add=0
scaleY_x='scale_y_log10()'
header='TRUE'
execute='TRUE'
ist='FALSE'
uwid=20
vhig=12
res=300
notch='FALSE'
par=''
outlier='FALSE'
out_scale=1.05
legend_pos='right'
color='FALSE'
ext='pdf'
violin='FALSE'
violin_nb='FALSE'
scale_violin='width'
ID_var=""
jitter='FALSE'
jitter_bp='FALSE'
colormodel='srgb'
rotate_plot='FALSE'
facet='NoMeAnInGTh_I_n_G_s'
nrow='NULL'
ncol='NULL'
scales='fixed'
facet_level='NA'
gene='FALSE'
sampleGroup='FALSE'

while getopts "ha:A:b:B:c:C:d:D:e:E:f:F:g:G:M:N:k:i:I:R:j:J:l:L:m:n:o:O:p:P:q:Q:r:s:S:t:u:v:V:w:W:x:y:z:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		m)
			melted=$OPTARG
			;;
		a)
			xvariable=$OPTARG
			;;
		d)
			value=$OPTARG
			;;
		F)
			variable=$OPTARG
			;;
		I)
			ID_var=$OPTARG
			;;
		j)
			jitter=$OPTARG
			;;
		J)
			jitter_bp=$OPTARG
			;;
		b)
			xtics_angle=$OPTARG
			;;
		G)
			facet=$OPTARG
			;;
		g)
			facet_level=$OPTARG
			;;
		q)
			gene=$OPTARG
			;;
		Q)
			sampleGroup=$OPTARG
			;;
		M)
			nrow=$OPTARG
			;;
		N)
			ncol=$OPTARG
			;;
		k)
			scales=$OPTARG
			;;
		t)
			title=$OPTARG
			;;
		x)
			xlab=$OPTARG
			;;
		l)
			level=$OPTARG
			;;
		P)
			legend_pos=$OPTARG
			;;
		B)
			x_cut=$OPTARG
			;;
		D)
			legend_cut=$OPTARG
			;;
		L)
			x_level=$OPTARG
			;;
		n)
			notch=$OPTARG
			;;
		V)
			violin=$OPTARG
			;;
		W)
			violin_nb=$OPTARG
			;;
		A)
			scale_violin=$OPTARG
			;;
		p)
			par=$OPTARG
			;;
		y)
			ylab=$OPTARG
			;;
		w)
			uwid=$OPTARG
			;;
		u)
			vhig=$OPTARG
			;;
		R)
			rotate_plot=$OPTARG
			;;
		r)
			res=$OPTARG
			;;
		E)
			ext=$OPTARG
			;;
		o)
			outlier=$OPTARG
			;;
		O)
			out_scale=$OPTARG
			;;
		s)
			scaleY=$OPTARG
			;;
		S)
			y_add=$OPTARG
			;;
		c)
			color=$OPTARG
			;;
		C)
			color_v=$OPTARG
			;;
		v)
			scaleY_x=$OPTARG
			;;
		z)
			header=$OPTARG
			;;
		e)
			execute=$OPTARG
			;;
		i)
			ist=$OPTARG
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

mid='.boxplot'

if test "${melted}" == "FALSE" && test "${sampleGroup}" == "FALSE"; then
	if test "${value}" != "value" || test "${variable}" != "variable"; then
		value="value"
		variable="variable"
		echo "Warning, there is no need to set -d and -F for unmelted \
files. We will ignore this setting and not affect the result."
	fi
fi

if test "${xvariable}" == ""; then
	xvariable=${variable}
fi

if test "${gene}" != "FALSE"; then
	value="${gene}"
fi

if test "${value}" != "value" || test "${variable}" != "variable"; then
	mid=${mid}'.'${value}_${variable}
fi
	
#if test "${ID_var}" != ""; then
#	ID_var=${ID_var}
#fi


if test "${outlier}" == "TRUE"; then
	mid=${mid}'.noOutlier'
fi

if test ${y_add} -ne 0; then
	scaleY="TRUE"
fi

#if test "${gene}" != "FALSE"; then
#	mid=${mid}".${gene}"
#fi

if test "${scaleY}" == "TRUE"; then
	mid=${mid}'.scaleY'
fi

if test "${violin}" == "TRUE"; then
	mid=${mid}'.violin'
fi
if test "${violin_nb}" == "TRUE"; then
	mid=${mid}'.violin_nb'
fi

if test "${jitter}" == "TRUE"; then
	mid=${mid}'.jitter'
fi

if test "${jitter_bp}" == "TRUE"; then
	mid=${mid}'.jitter_bp'
fi

. `dirname $0`/sp_configure.sh

cat <<END >${file}${mid}.r

if ($ist){
	install.packages("ggplot2", repo="http://cran.us.r-project.org")
	install.packages("reshape2", repo="http://cran.us.r-project.org")
	install.packages("scales", repo="http://cran.us.r-project.org")
	if(${jitter_bp}){
		install.packages("ggbeeswarm", repo="http://cran.us.r-project.org")
	}
}

if(${jitter_bp}){
	library(ggbeeswarm)
}else if(${jitter}){
	library(ggbeeswarm)
}

library(ggplot2)
library(reshape2)
library(scales)

if(! $melted){
	ID_var <- c("${ID_var}")
	ID_var <- ID_var[ID_var!=""]
	data <- read.table(file="${file}", sep="\t", header=$header,
	row.names=1, quote="", check.names=F)
	if ("${gene}" != "FALSE") {
		data_m <- as.data.frame(t(data["${gene}", ]))
		data_m\$sample = rownames(data_m)
		if ("${sampleGroup}" != "FALSE"){
			sampleGroup <- read.table("${sampleGroup}",sep="\t",header=1,check.names=F,row.names=1)
			data_m <- merge(data_m, sampleGroup, by="row.names")
		}
	} else {
		if ("$xvariable" != "${variable}"){
			if (length(ID_var) > 0){
				ID_var <- c(ID_var, "${xvariable}")
			} else {
				ID_var <- c("${xvariable}")
			}
			data_m <- melt(data, id.vars=ID_var)
		} else {
			if (length(ID_var) > 0){
				data_m <- melt(data, id.vars=ID_var)
			} else {
				data_m <- melt(data)
			}
		}
	}
} else {
	data_m <- read.table(file="$file", sep="\t",
	header=$header, quote="")
}

if (${y_add} != 0){
	data_m\$${value} <- data_m\$${value} + ${y_add}
}

level <- c(${level})

if ("${legend_cut}" != ""){
	data_m\$${variable} <- cut(data_m\$${variable}, ${legend_cut})
} else if (length(level)>1){
	level_i <- level
	data_m\$${variable} <- factor(data_m\$${variable}, levels=level_i)
}

x_level <- c(${x_level})

if ("${x_cut}" != ""){
	data_m\$${xvariable} <- cut(data_m\$${xvariable},${x_cut})
}else if (length(x_level)){
	data_m\$${xvariable} <- factor(data_m\$${xvariable},levels=x_level)
}

facet_level <- c(${facet_level})
if (length(facet_level)>1) {
	data_m\$${facet} <- factor(data_m\$${facet},
        levels=facet_level, ordered=T)
}


p <- ggplot(data_m, aes(factor($xvariable), ${value})) + xlab("$xlab") +
ylab("$ylab") + labs(title="$title")


if (${violin}){
	p <- p + geom_violin(aes(fill=factor(${variable})), 
	stat = "ydensity", position = "dodge", trim = TRUE,  
	scale = "${scale_violin}") + 
	geom_boxplot(aes(fill=factor(${variable})), alpha=.25, width=0.15, 
	position = position_dodge(width = .9), outlier.colour='NA') + 
	stat_summary(aes(group=${variable}), fun.y=mean,  
	geom="point", fill="black", shape=19, size=1,
	position = position_dodge(width = .9))
   	
	#+ geom_jitter(height = 0)
} else if (${violin_nb}){
	p <- p + geom_violin(aes(fill=factor(${variable})), 
	stat = "ydensity", position = "dodge", trim = TRUE,  
	scale = "${scale_violin}") 
} else if (${jitter}){
	p <- p + geom_quasirandom(aes(colour=factor(${variable})))
	p <- p + stat_summary(fun.y = "mean", geom = "text", label="----", size= 10, color= "black")
	#p <- p + geom_jitter(aes(colour=factor(${variable})))
} else {
	if (${notch}){
		if (${outlier}){
		p <- p + geom_boxplot(aes(fill=factor(${variable})), notch=TRUE,
			notchwidth=0.3, outlier.colour='NA')
		}else{
		p <- p + geom_boxplot(aes(fill=factor(${variable})), notch=TRUE,
			notchwidth=0.3)
		}
	}else {
		if (${outlier}){
			p <- p + geom_boxplot(aes(fill=factor(${variable})),
			outlier.colour='NA')
		}else{
			p <- p + geom_boxplot(aes(fill=factor(${variable})))
		}
	}
}

if (${jitter_bp}){
	#p <- p + geom_jitter(aes(colour=factor(${variable})))
	#p <- p + geom_jitter()
	p <- p + geom_quasirandom()
}

if($scaleY){
	p <- p + $scaleY_x
	p <- p + stat_summary(fun.y = "mean", geom = "text", label="----", size= 10, color= "black")
}

if(${outlier}){
	#ylim_zoomin <- boxplot.stats(data_m\$${value})\$stats[c(1,5)]
	stats <- boxplot.stats(data_m\$${value})\$stats
	ylim_zoomin <- c(stats[1]/${out_scale}, stats[5]*${out_scale})
	p <- p + coord_cartesian(ylim = ylim_zoomin)
}

if($color){
	p <- p + scale_fill_manual(values=c(${color_v}))
}

if(${rotate_plot}){
	p <- p + coord_flip()	
}

if ("${facet}" != "NoMeAnInGTh_I_n_G_s"){
	p <- p + facet_wrap( ~ ${facet}, nrow=${nrow}, ncol=${ncol},
	scale="${scales}")
}


END


`ggplot2_configure`

##cat <<END >>${file}${mid}.r
##
##p <- p + theme_bw() + theme(legend.title=element_blank(),
##	panel.grid.major = element_blank(), 
##	panel.grid.minor = element_blank(),
##	legend.key=element_blank(),
##	axis.text.x=element_text(angle=${xtics_angle},hjust=1))
##
##top='top'
##botttom='bottom'
##left='left'
##right='right'
##none='none'
##legend_pos_par <- ${legend_pos}
##
###if ("${legend_pos}" != "right"){
##p <- p + theme(legend.position=legend_pos_par)
###}
##
##
##p <- p${par}
##
##
##ggsave(p, filename="${file}${mid}.${ext}", dpi=$res, width=$uwid,
##height=$vhig, units=c("cm"))
##
###png(filename="${file}${mid}.png", width=$uwid, height=$vhig,
###res=$res)
###p
###dev.off()
##END

if [ "$execute" == "TRUE" ]; then
	Rscript ${file}${mid}.r
#if [ "$?" == "0" ]; then /bin/rm -f ${file}${mid}.r; fi
fi

