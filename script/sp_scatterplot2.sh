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

This script is used to do full functional scatter plot using ggplot2. 


The parameters for logical variable are either TRUE or FALSE.

Input file (terms only exist in one or a few samples are suitable):

Samp	X_val	Y_val	Color	Size	Shape
a	1	1	grp1	10	cluster1
b	2	2	grp1	10	cluster1
c	1	3	grp1	10	cluster1
d	3	1	grp2	15	cluster2
e	2	2	grp2	15	cluster2
f	3	3	grp3	5	cluster2
g	2	1	grp3	5	cluster2



**********************A potential bug******************************
If -c column have only 1 value, program will be aborted by no reasons.


${txtbld}OPTIONS${txtrst}:
	-f	Data file (with header line, the first column is not the
 		rowname, tab seperated)${bldred}[NECESSARY]${txtrst}
	-T	A tag to define the plot. Normally used when one file is used
		to do multiple plots with different parameters to discriminate
		different output file.
		[${txtred}Default empty; Optional${txtrst}]
	-t	Title of picture[${txtred}Default empty title${txtrst}]
		[Scatter plot of horizontal and vertical variable]
	-x	xlab of picture[${txtred}Default empty xlab${txtrst}]
		[The description for horizontal variable]
	-y	ylab of picture[${txtred}Default empty ylab${txtrst}]
		[The description for vertical variable]
	-p	Legend position [Lowercase p]
		[${txtred}Default right. Accept
		top, bottom, left, none,  or c(0.08, 0.8).${txtrst}]
	-R	Rotation angle for x-axis value(anti clockwise)
		[Default 0]
	-X	The variable for horizontal axis. [Uppercase X]
		${bldred}[NECESSARY, such X_val, both text and number works]${txtrst}
	-O	The order for horizontal axis. [Uppercease O]
		${bldred}[Default alphabetical order, accept a string like
		"'K562','hESC','GM12878','HUVEC','NHEK','IMR90','HMEC'"]
		${txtrst}
	-Y	The variable for vertical axis.
		${bldred}[NECESSARY, such as Y_val, both text and number works]${txtrst}
	-B 	The order for vertical axis.
		${bldred}[Default alphabetical order, accept a string like
		"'K562','hESC','GM12878','HUVEC','NHEK','IMR90','HMEC'"]
		${txtrst}
	-c	The variable for point color.
		${bldred}[Optional, such as color]${txtrst}
	-I	The order for color variable.
		${bldred}[Default alphabetical order, accept a string like
		"'K562','hESC','GM12878','HUVEC','NHEK','IMR90','HMEC'"]
		${txtrst}
	-s	The variable for point size.
		${bldred}[Optional, such as a number or 
		a variable like count, normally should be number column]${txtrst}
	-S	The variable for point shape.
		${bldred}[Optional, such as shape]${txtrst}
	-K	The order for shape variable.
		${bldred}[Default alphabetical order, accept a string like
		"'K562','hESC','GM12878','HUVEC','NHEK','IMR90','HMEC'"]
		${txtrst}
	-C	Manually specified colors.
		${bldred}[Default system default. 
		Accept string in format like <'"green", "red"'> (both types of quotes needed).]
		${txtrst}
	-a	Type of color variable. 
		<num> represents color variables are numbers. At this time, 
		two colors should be specified in <-C> represents low and high
		color.
		<factor> represents color variables are strings. At this time, 
		same number of colors as number of samples should be given.	
	-A	Transparency value for points.
		${bldred}[Optional, such as a number or 
		a variable indicating one data column, 
		normally should be number column]${txtrst}
	-l	Get log-transformed data for given variable. [Lowercase l]
		[${txtred}Default nolog, means no log10 transform. Accept a variable
		like color to get (-1) * log10(color).${txtrst}]	
	-L	Label points.
		${bldred}[Default no-label, accept a string like <Samp> here 
		to label Samp column text to points.]${txtrst}
	-Z	Label points using <geom_text_repel> which will generate
   		non-overlap label texts by adding arrows when necessary.
		If there is sth wrong labeled especially when -J is TRUE, 
		please specify FALSE and use default
		<geom_text> to label text.
		${bldred}[Default TRUE.]${txtrst}
	-N	Label font size.
		${bldred}[Default system default. Accept a number.]${txtrst}
	-M	Points check_overlap.
		${bldred}[Optional, such as shape]${txtrst}
	-Q	Point hjust.
		[Default 0,  accept a positive (at left) and negative value (at right)]
	-J	Jitter points. Normally used when x and y axis variable is in text format 
		or represents group information to avoid point overlaps.
		${bldred}[Default FALSE]${txtrst}
	-D	Scale y axis.
		${bldred}[Default FALSE]${txtrst}
	-F	The way to scale Y-axis.
		${bldred}[scale_y_log10, coord_trans(y="log10"), 
		scale_y_continuous(trans="log2")(default), coord_trans(y="log2")
		]${txtrst}
	-w	The width of output picture.[${txtred}Default 20${txtrst}]
	-u	The height of output picture.[${txtred}Default 20${txtrst}] 
	-E	The type of output figures.[${txtred}Default pdf, accept
		eps/ps, tex (pictex), png, jpeg, tiff, bmp, svg and wmf)${txtrst}]
	-r	The resolution of output picture.[${txtred}Default 300 ppi${txtrst}]
	-b	The formula for facets.[${bldred}Default no facets, 
		+facet_grid(level ~ .) means divide by levels of 'level' vertcally.
		+facet_grid(. ~ level) means divide by levels of 'level' horizontally.
		+facet_grid(lev1 ~ lev2) means divide by lev1 vertically and lev2
		horizontally.
		+facet_wrap(~level, ncol=2) means wrap horizontally with 2
		columns.
		Example: +facet_wrap(~Size,ncol=6,scale='free')
		${txtrst}]
	-d	If you may want to specifize the order of
		other variables (default alphabetically), please supply a string like below.
		[${txtred}Accept sth like 
		(one level one sentence, separate by';') 
		data\$size <- factor(data\$size, levels=c("l1",
		"l2",...,"l10"), ordered=T) ${txtrst}]
	-z	Other parameters in ggplot format.
		[${bldred}optional${txtrst}]
	-e	Execute or not[${bldred}Default TRUE${txtrst}]
	-i	Install the required packages[${bldred}Default FALSE${txtrst}]


s-plot scatterplot2 -f file -X x_val -Y y_val -c color -s size -S shape -L Samp

EOF
}

file=''
title=''
xlab=''
ylab=''
xval=''
xval_order=''
yval_order=''
yval=''
execute='TRUE'
ist='FALSE'
color='c_t_c_t0304'
color_v=''
color_t=''
log='nolog'
uwid=20
vhig=20
res=300
ext='pdf'
facet=''
tag=''
size=''
geom_text_repel='TRUE'
shape='c_t_c_t0304'
par=''
variable_order=''
color_order=''
shape_order=''
legend_pos='right'
xtics_angle=0
hjust=0.5
point_hjust=0
vjust=1
alpha=1
jitter='FALSE'
label=''
check_overlap="FALSE"
colormodel='srgb'
label_font_size=0
scale_y='FALSE'
scale_y_way='scale_y_continuous(trans="log2")'

while getopts "hf:t:a:x:y:p:X:O:R:Y:Z:B:H:V:L:T:I:K:v:c:C:A:l:D:F:N:L:M:J:w:u:r:E:s:S:b:d:z:e:i:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		t)
			title=$OPTARG
			;;
		T)
			tag=$OPTARG
			;;
		x)
			xlab=$OPTARG
			;;
		y)
			ylab=$OPTARG
			;;
		p)
			legend_pos=$OPTARG
			;;
		R)
			xtics_angle=$OPTARG
			;;
		X)
			xval=$OPTARG
			;;
		O)
			xval_order=$OPTARG
			;;
		Y)
			yval=$OPTARG
			;;
		B)
			yval_order=$OPTARG
			;;
		c)
			color=$OPTARG
			;;
		C)
			color_v=$OPTARG
			;;
		a)
			color_t=$OPTARG
			;;
		A)
			alpha=$OPTARG
			;;
		l)
			log=$OPTARG
			;;
		L)
			label=$OPTARG
			;;
		M)
			check_overlap=$OPTARG
			;;
		Z)
			geom_text_repel=$OPTARG
			;;
		N)
			label_font_size=$OPTARG
			;;
		I)
			color_order=$OPTARG
			;;
		K)
			shape_order=$OPTARG
			;;
		D)
			scale_y=$OPTARG
			;;
		F)
			scale_y_way=$OPTARG
			;;
		w)
			uwid=$OPTARG
			;;
		u)
			vhig=$OPTARG
			;;
		r)
			res=$OPTARG
			;;
		E)
			ext=$OPTARG
			;;
		b)
			facet=$OPTARG
			;;
		d)
			variable_order=$OPTARG
			;;
		s)
			size=$OPTARG
			;;
		J)
			jitter=$OPTARG
			;;
		L)
			point_hjust=$OPTARG
			;;
		S)
			shape=$OPTARG
			;;
		z)
			par=$OPTARG
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

mid=".scatterplot2"${tag}
if [ -z $file ] || [ -z $xval ] || [ -z $yval ] ; then
	echo 1>&2 "Please give filename, xval and yval."
	usage
	exit 1
fi

. `dirname $0`/sp_configure.sh

cat <<END >${file}${mid}.r

if ($ist){
	install.packages("ggplot2", repo="http://cran.us.r-project.org")
	install.packages("ggbeeswarm", repo="http://cran.us.r-project.org")
	install.packages("ggrepel", repo="http://cran.us.r-project.org")
}
library(plyr)
library(ggplot2)
library(grid)

if (${jitter}) {
	library(ggbeeswarm)
}

if (("${label}" != "") & (${geom_text_repel}) ) {
	library('ggrepel')
}
data <- read.table(file="$file", sep="\t", quote="", comment="", header=T)

yval_order <- c(${yval_order})

if (length(yval_order) > 1) {
	data\$${yval} <- factor(data\$${yval}, levels=yval_order, ordered=T)
}
xval_order <- c(${xval_order})

if (length(xval_order) > 1) {
	data\$${xval} <- factor(data\$${xval}, levels=xval_order, ordered=T)
}

shape_order <- c(${shape_order})

if (length(shape_order) > 1) {
	data\$${shape} <- factor(data\$${shape}, levels=shape_order, ordered=T)
}

color_order <- c(${color_order})

if (length(color_order) > 1) {
	data\$${color} <- factor(data\$${color}, levels=color_order, ordered=T)
}


if ("${log}" != "nolog"){
	data\$${log} <- log10(data\$${log}) * (-1)
}

$variable_order

color_v <- c(${color_v})

if ("${shape}" != "c_t_c_t0304") {
	shape_level <- length(unique(data\$${shape}))
	shapes = (1:shape_level)%%30
}

p <- ggplot(data, aes(x=${xval},y=${yval}))

if ("$xlab" == "") {
	xlab = "${xval}"
} else {
	xlab = "$xlab"
}

if ("$ylab" == "") {
	ylab = "${yval}"
} else {
	ylab = "$ylab"
}

p <- p + labs(x=xlab, y=ylab) + labs(title="$title")


if (${jitter}) {
	if (("${size}" != "") && ("${color}" != "c_t_c_t0304") && ("${shape}" != "c_t_c_t0304")) {
		p <- p + geom_quasirandom(aes(size=${size}, color=${color}, shape=${shape}), alpha=${alpha}) 
	} else if (("${color}" != "c_t_c_t0304") && ("${shape}" != "c_t_c_t0304")) {
		p <- p + geom_quasirandom(aes(color=${color}, shape=${shape}), alpha=${alpha}) 
	} else if (("${size}" != "") && ("${shape}" != "c_t_c_t0304")) {
		p <- p + geom_quasirandom(aes(size=${size}, shape=${shape})) 
	} else if (("${size}" != "") && ("${color}" != "c_t_c_t0304")) {
		p <- p + geom_quasirandom(aes(size=${size}, color=${color}), alpha=${alpha}) 
	} else if ("${size}" != "") {
		p <- p + geom_quasirandom(aes(size=${size}))
	} else if ("${color}" != "c_t_c_t0304") {
		p <- p + geom_quasirandom(aes(color=${color}), alpha=${alpha})
	} else if ("${shape}" != "c_t_c_t0304") {
		p <- p + geom_quasirandom(aes(shape=${shape}))
	} else {
		p <- p + geom_quasirandom()
	}
} else {
	if (("${size}" != "") && ("${color}" != "c_t_c_t0304") && ("${shape}" != "c_t_c_t0304")) {
		p <- p + geom_point(aes(size=${size}, color=${color}, shape=${shape}), alpha=${alpha}) 
	} else if (("${color}" != "c_t_c_t0304") && ("${shape}" != "c_t_c_t0304")) {
		p <- p + geom_point(aes(color=${color}, shape=${shape}), alpha=${alpha}) 
	} else if (("${size}" != "") && ("${shape}" != "c_t_c_t0304")) {
		p <- p + geom_point(aes(size=${size}, shape=${shape})) 
	} else if (("${size}" != "") && ("${color}" != "c_t_c_t0304")) {
		p <- p + geom_point(aes(size=${size}, color=${color}), alpha=${alpha}) 
	} else if ("${size}" != "") {
		p <- p + geom_point(aes(size=${size}))
	} else if ("${color}" != "c_t_c_t0304") {
		p <- p + geom_point(aes(color=${color}), alpha=${alpha})
	} else if ("${shape}" != "c_t_c_t0304") {
		p <- p + geom_point(aes(shape=${shape}))
	} else {
		p <- p + geom_point()
	}
}


if (("${color}" != "c_t_c_t0304") && "${color_t}" == "factor" && length(color_v) > 1) {
	p <- p + scale_color_manual(values=color_v)
} 

if (("${color}" != "c_t_c_t0304") && "${color_t}" == "num" && length(color_v) > 1) {
	p <- p + scale_colour_gradient(low=color_v[1], high=color_v[2], name="${color}")
}



if (("${shape}" != "c_t_c_t0304") && shape_level > 6) {
	p <- p + scale_shape_manual(values=shapes)
}


if ("${label}" != "") {
	if (${geom_text_repel}) {
		if (${label_font_size} != 0) {
			p <- p + geom_text_repel(aes(label=${label}), size=${label_font_size})
		} else {
			p <- p + geom_text_repel(aes(label=${label}))
		}
	} else {
		if (${jitter}) {
			if (${label_font_size} != 0) {
				p <- p + geom_text(aes(label=${label}), 
				position=position_quasirandom(), 
				hjust=${point_hjust}, size=${label_font_size}, 
				check_overlap=${check_overlap})
			} else {
				p <- p + geom_text(aes(label=${label}), 
				position=position_quasirandom(), 
				hjust=${point_hjust}, check_overlap=${check_overlap})
			}
		} else {
			if (${label_font_size} != 0) {
				p <- p + geom_text(aes(label=${label}), 
				position="identity",
				hjust=${point_hjust}, size=${label_font_size}, 
				check_overlap=${check_overlap})
			} else {
				p <- p + geom_text(aes(label=${label}), 
				position="identity",
				hjust=${point_hjust}, size=${label_font_size}, 
				check_overlap=${check_overlap})
			}
		}
	}
}

if (${scale_y}) {
	p <- p + ${scale_y_way}
}

p <- p ${facet}

END

`ggplot2_configure`

if [ "$execute" == "TRUE" ]; then
	Rscript ${file}${mid}.r
#if [ "$?" == "0" ]; then /bin/rm -f ${file}${mid}.r; fi
fi

if [ ! -z "$log" ]; then
	log=', trans=\"'$log'\"'
fi

#convert -density 200 -flatten ${file}${mid}.eps ${first}${mid}.png
