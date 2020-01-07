#!/usr/bin/env R
# Sun Jan 27 18:58:34 2019
# Boyuan_Li

# need the non gene reads bed file and gene reads bed file
# calculate significant according to this;

command=matrix(c(
	"non_gene_bed", "b", "1", "character", 'non gene reads bed file',
	"gene_reads_bed", "g", "1", "character", 'gene reads bed file',
	"output_file_name", "o", "1", "character", 'output file name',
	"output_dir", "d", "1", "character" ,'output dir'
	),byrow=T,ncol=5)


args=getopt::getopt(command)


if (is.null(args$non_gene_bed) || is.null(args$gene_reads_bed) || is.null(args$output_file_name) || is.null(args$output_dir)) {
	cat(paste(getopt::getopt(command, usage = T), "\n"))
	q()
}

library(data.table)
library(dplyr)
#--- calculate lamada
non_gene_reads_num <- fread(args$non_gene_bed,
                            col.names = c('chr','st','en','id','l','strand','num'))

lamada <- sum(non_gene_reads_num$num)/sum(non_gene_reads_num$en-non_gene_reads_num$st)  

#--- calculate the sig gene and sig for gene

gene_file <- args$gene_reads_bed
sig_calcu <- function(filename){
  ext_gene <- fread(filename,
                    col.names = c('chr','st','en','id','l','strand','num'))  
  ext_gene <- ext_gene %>%
    mutate(exp=round(lamada*(en-st)))
  
  ext_gene$palue <- 'nosig'
  for (i in 1:nrow(ext_gene)){
    ext_gene$palue[i] <- poisson.test(ext_gene$num[i],ext_gene$exp[i],alternative = 'greater',conf.level = 0.95)$p.value
  }
  return(ext_gene)  
}
ext_gene <- sig_calcu(gene_file)

setwd(args$output_dir)
fwrite(ext_gene,args$output_file_name,sep="\t",col.names = F)

