#!/usr/bin/env R
# Sun Jan 27 18:45:17 2019
# Boyuan_Li

# due to the select >500K non gene region contain many region NNNNNN, so I will remove the
# first region and last region for each chr
#output_dir <- './'
command=matrix(c(
	"non_gene_bed", "b", "1", "character", 'should be the non gene region >500K',
	"output_file_name", "o", "1", "character", 'output file name that will be used',
	"output_dir", "d", "1", "character" ,'output dir'
	),byrow=T,ncol=5)


args=getopt::getopt(command)


if (is.null(args$non_gene_bed) || is.null(args$output_file_name) || is.null(args$output_dir)) {
	cat(paste(getopt::getopt(command, usage = T), "\n"))
	q()
}

library(data.table)
library(dplyr)
#--- rm fist and last line for each chr
non_gene_region <- fread(args$non_gene_bed)  
non_gene_region <- non_gene_region %>%
  filter(V2!=0)
chr_max <- c()
a=0
for (i in paste0('chr',c(1:19,'X','Y'))){
  a = a+1
  chr_max[a] <- max(non_gene_region[non_gene_region$V1==i,]$V3)
  
}
non_gene_region <- non_gene_region %>%
  filter(!(V3  %in% chr_max))
setwd(args$output_dir)
fwrite(non_gene_region,args$output_file_name,sep = '\t',col.names = F)

