#!/usr/bin/env R
# Thu Mar  7 21:37:00 2019
# Boyuan_Li


command=matrix(c(
	"bowtie_mpfile", "b", "1", "character", 
	"stastic_file", "s", "1", "character", 
	"output_filename", "o", "1", "character" 
	),byrow=T,ncol=4)


args=getopt::getopt(command)


if (is.null(args$bowtie_mpfile) || is.null(args$stastic_file) || is.null(args$output_filename)) {
	cat(paste(getopt::getopt(command, usage = T), "\n"))
	q()
}
library(dplyr)
library(data.table)


mapping_stas <- args$bowtie_mpfile

mapping_file <- suppressMessages(suppressWarnings(readr::read_csv(mapping_stas,col_names = F)))
#mapping_file <- mapping_file[is.na(mapping_file)]
mapping_file <- as.data.frame(mapping_file)
mapping_file$X1 <- as.numeric(mapping_file$X1)
mapping_ratio <- (2*(mapping_file[4,1] + mapping_file[5,1] + mapping_file[8,1])+mapping_file[13,1] + mapping_file[14,1])/(2*mapping_file[1,1])
mapping_reads <- 2*(mapping_file[4,1] + mapping_file[5,1] + mapping_file[8,1])+mapping_file[13,1] + mapping_file[14,1]
unique_mapping_reads <- 2*(mapping_file[4,1] + mapping_file[8,1])+mapping_file[13,1] 
unique_mapping_ratio <- (2*(mapping_file[4,1] + mapping_file[8,1])+mapping_file[13,1] )/(2*mapping_file[1,1])
reads_after_cut <- mapping_file[1,1]
# (2*(379520+115980+125180)+204075+130277)/2000000
# (2*(379520+125180)+204075)/2000000
stas_file <- args$stastic_file
stas_flie2 <- read.csv(stas_file,header = T,stringsAsFactors = F)
bowtie_sts <- data.frame('unique mapping reads'=unique_mapping_reads,
                         'mapping readsb'=mapping_reads,
                         'cutadaptor'=reads_after_cut)
all_stas <- cbind(stas_flie2,bowtie_sts)
rownames(all_stas) <- 'counts'
all_stas <- all_stas %>%
  select(sample_name,total_resds,cutadapt,mapping.reads,unique.mapping.reads,rm.duplicate,final.unique)
all_stas2 <- data.frame('sample_name'=c(all_stas$sample_name,''),
                        'total_reads'=c(all_stas$total_resds,all_stas$total_resds/all_stas$total_resds),
                        'cutadapt'=c(all_stas$cutadapt,all_stas$cutadapt/all_stas$total_resds),
                        'mapping.reads'=c(all_stas$mapping.reads,all_stas$mapping.reads/(all_stas$total_resds*2)),
                        'rm.duplicate'=c(all_stas$rm.duplicate,all_stas$rm.duplicate/(all_stas$total_resds*2)),
                        'unique.mapping.reads'=c(all_stas$unique.mapping.reads,all_stas$unique.mapping.reads/(all_stas$total_resds*2)),
                        'final.unique'=c(all_stas$final.unique,all_stas$final.unique/(all_stas$total_resds*2)))

#setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/test')
filename <- args$output_filename
fwrite(all_stas2,filename,quote = F)

