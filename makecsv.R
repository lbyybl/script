#!/usr/bin/R

# Boyuan-Li

# 2018-7-1

# This script is used to trans txt into csv

library(getopt)
command=matrix(c("Input","i",1,"character",

                 "Output","o",1,"character",

                 "help","h",0,"logical"),byrow=T,ncol=4)

args=getopt(command)

if (!is.null(args$help) || is.null(args$Input) || is.null(args$Output)) {

    cat(paste(getopt(command, usage = T), "\n"))

    q()

} 

#setwd(args$Input)
#data<-read.table(args$Input)
data<-read.table(args$Input,header = TRUE, stringsAsFactors = FALSE)

threadhood<-nrow(data) + 1

#data <- data[,2:ncol(data)]
#fix(data)

data <- subset(data, select = -sample_name.1)
data <- subset(data, select = -sample_name.2)

i=2
while (i <= threadhood){
  data[i,2:ncol(data)]<-as.numeric(as.character(data[(i-1),2:ncol(data)]))/as.numeric(as.character(data[(i-1),2]))
  i=i+2
}


write.csv(data,args$Output)
