#!/usr/bin/env R
# Fri Feb 22 15:05:56 2019
# Boyuan_Li

# this sript is used to merge the wt loop bedpe and ko bedpe
# defult if the distance between tow loci is short than 10000
# they will be merge to one
reso=10000
command=matrix(c(
	"Input_file", "i", "1", "character", "input filename",
	"Output_file", "o", "1", "character", "output filename",
	"Resolution", "r", "2", "integer" ,"defult is 10000"
	),byrow=T,ncol=5)


args=getopt::getopt(command)


if (is.null(args$Input_file) || is.null(args$Output_file) || is.null(args$Resolution)) {
	cat(paste(getopt::getopt(command, usage = T), "\n"))
	q()
}
file <- args$Input_file
data <- data.table::fread(file)
reso=args$Resolution
library(dplyr)

unique_data <- function(data,reso){
  data <- arrange(data,V1,V2,V3,V4,V5,V6) 
  for (i in 1:(nrow(data)-1)){
    if (data[i,2] >= (data[i+1,2]-reso) & (data[i,1]==data[i+1,1])){
      if (data[i,5] >= (data[i+1,5]-reso) & (data[i,4]==data[i+1,4])){
        data[i,2] <- min(data[i,2],data[i+1,2])
        data[i+1,2] <- min(data[i,2],data[i+1,2])
        data[i,3] <- max(data[i,3],data[i+1,3])
        data[i+1,3] <- max(data[i,3],data[i+1,3])
        data[i,5] <- min(data[i,5],data[i+1,5])
        data[i+1,5] <- min(data[i,5],data[i+1,5])
        data[i,6] <- max(data[i,6],data[i+1,6])
        data[i+1,6] <- max(data[i,6],data[i+1,6])
      }
    }
  }
  data <- unique(data)
  return(data)
}
data <- arrange(data,V1,V2,V3,V4,V5,V6)
data <- unique(data)
output_filename <- args$Output_file
data.table::fwrite(data,output_filename,sep = "\t",col.names = F)


