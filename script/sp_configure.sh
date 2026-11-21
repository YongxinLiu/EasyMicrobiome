#!/bin/bash

function ggplot2_configure {

cat <<END >>${file}${mid}.r

#Configure the canvas
#legend.title=element_blank(),
p <- p + theme_bw() + theme(
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(),
	legend.key=element_blank())

if (${xtics_angle} != 0){
	if (${xtics_angle} == 90){
		p <- p + theme(axis.text.x=
		  element_text(angle=${xtics_angle},hjust=1, vjust=0.5))
	}else if (${xtics_angle} == 45){
		p <- p + theme(axis.text.x=
		  element_text(angle=${xtics_angle},hjust=0.5, vjust=0.5))
	} else {
		p <- p + theme(axis.text.x=
		  element_text(angle=${xtics_angle},hjust=0.5, vjust=0.5))
	}
}

#Set the position of legend
top='top'
botttom='bottom'
left='left'
right='right'
none='none'
legend_pos_par <- ${legend_pos}

p <- p + theme(legend.position=legend_pos_par)

#add additional ggplot2 supported commands

p <- p${par}


# output pictures

if ("${ext}" == "pdf") {
	ggsave(p, filename="${file}${mid}.${ext}", dpi=$res, width=$uwid,
	height=$vhig, units=c("cm"),colormodel="${colormodel}")
} else {
	ggsave(p, filename="${file}${mid}.${ext}", dpi=$res, width=$uwid,
	height=$vhig, units=c("cm"))
}
#png(filename="${file}${mid}.png", width=$uwid, height=$vhig,
#res=$res)
#p
#dev.off()
END
}
