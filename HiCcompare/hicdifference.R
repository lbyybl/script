#!/usr/bin/env R
# 2018/11/7
# Boyuan_li
command=matrix(c("Wt_file","w",1,"character",

                 "Ko_file","k",1,"character",
				 
				 "Output_file","o",1,"character",
				 
				 "Chr","c",1,"character",

                 "help","h",0,"logical"),byrow=T,ncol=4)

args=getopt::getopt(command)

if (!is.null(args$help) || is.null(args$Wt_file) || is.null(args$Ko_file) || is.null(args$Output_file) || is.null(args$Chr)) {

    cat(paste(getopt::getopt(command, usage = T), "\n"))

    q()

}


## ----set-options, echo=FALSE, cache=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- HiCcompare不使用inter interaction的信息所以灭有办法找出染色体间的互做差异
options(width = 400)
options(stringsAsFactors = F)

#library(data.table)
library(HiCcompare)
	wtdata <- data.table::fread(args$Wt_file)
	kodata <- data.table::fread(args$Ko_file)
	blaklist <- "/DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed"
	exclude <- read.table(blaklist, col.names=c("chr","start","end"))

	hic.table <- create.hic.table(wtdata, kodata, chr = args$Chr, exclude.regions = exclude, exclude.overlap = 0.2)
#	hic.table <- create.hic.table(wtdata, kodata, chr = args$Chr)

	hic.table <- hic_loess(hic.table, Plot = FALSE, Plot.smooth = FALSE)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- 差异分析
hic.table <- hic_compare(hic.table, A.min = 15, adjust.dist = TRUE, p.method = 'fdr', Plot = TRUE)
hic.tablesignifent <- hic.table[hic.table$p.adj < 0.05,]
data.table::fwrite(hic.table,file = args$Output_file, sep = "\t")
data.table::fwrite(hic.tablesignifent,file = paste0("signifent",args$Output_file), sep = "\t")
rm(list=ls())
