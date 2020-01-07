#!/usr/bin/env R
# Sat Nov 16 23:27:30 2019
# Boyuan_Li


command=matrix(c(
	"input_dir", "d", "1", "character", 
	"dup_bam", "b", "1", "character", 
	"unique_bam", "u", "1", "character" 
	),byrow=T,ncol=4)


args=getopt::getopt(command)


if (is.null(args$input_dir) || is.null(args$dup_bam) || is.null(args$unique_bam)) {
	cat(paste(getopt::getopt(command, usage = T), "\n"))
	q()
}

## ---- echo=FALSE, results="hide", warning=FALSE, message=FALSE-------------
suppressPackageStartupMessages({
  library(ATACseqQC)
  library(data.table)
  library(GenomicAlignments)
})

## --------------------------------------------------------------------------
setwd(args$input_dir)
# input the bamFile from the ATACseqQC package
bamfile <- args$dup_bam
bam_index <- paste0(bamfile,'.bai')
bamfile.labels <- gsub(".bam", "", basename(bamfile))
# 
# ## ---- eval=FALSE-----------------------------------------------------------
# #  source(system.file("extdata", "IGVSnapshot.R", package = "ATACseqQC"))
# 
# ## --------------------------------------------------------------------------
QC_list <- bamQC(bamfile,bam_index,mitochondria = 'chrM', outPath=paste0(bamfile.labels,'_clean.bam'))
con <-  file(paste0(bamfile.labels,"mtcars.txt"), open = "w")
writeLines(text = paste0("totalQNAMEs ", QC_list$totalQNAMEs), con = con )
writeLines(text = paste0("duplicateRate ", QC_list$duplicateRate), con = con )
writeLines(text = paste0("mitochondriaRate ", QC_list$mitochondriaRate), con = con )
writeLines(text = paste0("properPairRate ", QC_list$properPairRate), con = con )
writeLines(text = paste0("unmappedRate ", QC_list$unmappedRate), con = con )
writeLines(text = paste0("hasUnmappedMateRate ", QC_list$hasUnmappedMateRate), con = con )
writeLines(text = paste0("notPassingQualityControlsRate ", QC_list$notPassingQualityControlsRate), con = con )
writeLines(text = paste0("nonRedundantFraction ", QC_list$nonRedundantFraction), con = con )
writeLines(text = paste0("PCRbottleneckCoefficient_1 ", QC_list$PCRbottleneckCoefficient_1), con = con )
writeLines(text = paste0("PCRbottleneckCoefficient_2 ", QC_list$PCRbottleneckCoefficient_2), con = con )
write.table(x=QC_list$MAPQ, file = con, quote = FALSE, sep = "\t")
write.table(x=QC_list$idxstats, file = con, quote = FALSE, sep = "\t")
close(con)

duplicate_freq <- readsDupFreq(bamfile,index = bam_index)
fwrite(as.data.frame(duplicate_freq),paste0(bamfile.labels,'duplicate_frep.txt'),sep = '\t')
pdf(paste0(bamfile.labels,'estimate_complexity.pdf'))
plot <- estimateLibComplexity(duplicate_freq)
print(plot)
fwrite(plot,paste0(bamfile.labels,'estimate_complexity.txt'),sep = '\t')
dev.off()

## --------------------------------------------------------------------------
## generate fragement size distribution
pdf(paste0(bamfile.labels,'_fragement_size.pdf'))
fragSize <- fragSizeDist(bamfile,index = bam_index, bamFiles.labels = bamfile.labels)
print(fragSize)
fragSize2 <- as.data.frame(fragsize)
colnames(fragSize2) <- c('len','num')
fwrite(fragSize2,paste0(bamfile.labels,'fragment_length.txt'),sep = '\t')
dev.off()
FL<-sum(as.numeric(fragSize2$len)*as.numeric(fragSize2$num))/sum(as.numeric(fragSize2$num))
QC_summary <- data.frame('sample_name'=bamfile.labels,
                         'NRF'=QC_list$nonRedundantFraction,
                         'PBC1'=QC_list$PCRbottleneckCoefficient_1,
                         'PBC2'=QC_list$PCRbottleneckCoefficient_2,
                         'frag len'=FL)
fwrite(QC_summary,paste0(bamfile.labels,'QC_summary.txt'),sep = '\t')