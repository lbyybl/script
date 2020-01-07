#!/usr/bin/env R
# Tue Mar 26 11:06:12 2019
# Boyuan_Li


command=matrix(c(
	"path", "p", "1", "character", 
	"outputname", "o", "1", "character", 
	"sample_name", "n", "1", "character" 
	),byrow=T,ncol=4)


args=getopt::getopt(command)


if (is.null(args$path) || is.null(args$outputname) || is.null(args$sample_name)) {
	cat(paste(getopt::getopt(command, usage = T), "\n"))
	q()
}
#--- stastic the output of HiC pro
path <- args$path
outputname <- args$outputname
samplename <- args$sample_name
suppressMessages(library(data.table))
setwd(path)
sts <- read.table('tatal.stastic',header = F)
class(t(sts))
sts_t <- t(sts)
colnames(sts_t) <- t(sts)[1,]

for (i in 1:ncol(sts_t)){
  sts_t[1,i] <- as.numeric(sts_t[2,i])/as.numeric(sts_t[2,"Total_pairs_processed"])
}
sts_t <- data.frame(sts_t)
sts_t <- sts_t[,1:ncol(sts_t)]
sts_t$sample_name <- samplename
sts_t2 <- sts_t
sts_t2[1,] <- sts_t[2,]
sts_t2[2,] <- sts_t[1,]
sts_t <- sts_t2
sts_t <-  dplyr::select(sts_t,sample_name,Total_pairs_processed,Unmapped_pairs,Pairs_with_singleton,Multiple_pairs_alignments,Unique_paired_alignments,
                        Dangling_end_pairs,Self_Cycle_pairs,Dumped_pairs,Valid_interaction_pairs,
                        valid_interaction_rmdup,cis_shortRange,cis_longRange,trans_interaction)
fwrite(sts_t,args$outputname,col.names = T,sep = ",")


