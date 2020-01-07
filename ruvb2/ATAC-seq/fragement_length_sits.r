#!/usr/bin/env R
# Tue Jun 25 15:57:16 2019
# Boyuan_Li


command=matrix(c(
	"bam_file", "b", "1", "character", "bamfile",
	"bam_index_file", "d", "1", "character", "bamindex file",
	"output_pdf_name", "o", "1", "character" ,"output file name such as 'haha.pdf'"
	),byrow=T,ncol=5)


args=getopt::getopt(command)


if (is.null(args$bam_file) || is.null(args$bam_index_file) || is.null(args$output_pdf_name)) {
	cat(paste(getopt::getopt(command, usage = T), "\n"))
	q()
}


suppressPackageStartupMessages(library(ATACseqQC))

bamfile <- args$bam_file
bamindex <- gsub(".bai", "", args$bam_index_file)

output_file_name <- args$output_pdf_name
bamfile.labels <- gsub(".pdf","", args$output_pdf_name)
pdf(output_file_name,width = 16,height = 6)
fragSize <- fragSizeDist(bamfile, bamfile.labels,index = bamindex)
dev.off()
