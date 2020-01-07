#!/usr/bin/env R
# Thu Apr  4 17:08:32 2019
# Boyuan_Li


command=matrix(c(
	"input", "i", "1", "character" 
	),byrow=T,ncol=4)


args=getopt::getopt(command)


if (is.null(args$input)) {
	cat(paste(getopt::getopt(command, usage = T), "\n"))
	q()
}

#--- it's used to extract the result run by samtools flagstat
#setwd('/DATA/work/lbyybl/wh/ruvb2/SLAM/slam_output/map')
#library(readr)
suppressMessages(suppressWarnings(flag_file <- readr::read_csv(args$input,col_names = F)))
#View(flag_file)
a<-flag_file[5,1]

print(as.numeric(a))