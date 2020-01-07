#!/usr/bin/env R
# Sun Jan 27 19:11:56 2019
# Boyuan_Li

# need non gene reads bed to calculate lamada; gene(rm TSS) bed and promoter reads bed 
# to calcu paused gene

command=matrix(c(
	"promoter_bed", "b", "1", "character", 'prmoter bed file with reads number in each 50bp region',
	"gene_reads_bed", "g", "1", "character", 'gene bed file with reads number in each region',
	"output_dir", "d", "1", "character", 'output dir',
	"prefix", "p", "1", "character", 'the output file prefix',
	'non_gene_bed','n',1,'character','non gene bed file with reads num to calcu lamada',
	'gene_bed_ext','s',1,'character','the bed file with reads for extend gene',
	'unique_all_gene_bed','u',1,'character','unique all gene bed'
	),byrow=T,ncol=5)


args=getopt::getopt(command)


if (is.null(args$promoter_bed) || is.null(args$gene_reads_bed) || 
    is.null(args$gene_bed_ext) || is.null(args$output_dir) || 
    is.null(args$prefix) || is.null(args$unique_all_gene_bed) || is.null(args$non_gene_bed)) {
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
gene_file <- args$gene_bed_ext
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
#--- calculate paused promoter
promoter_file <- args$promoter_bed
prmoter <- fread(promoter_file,
                 col.names = c('pchr','pst','pen','id','pl','pstrand','pnum'))  
#--- gene short 1k 
gene_body_file <- args$gene_reads_bed
gene_body <- fread(gene_body_file,
                   col.names = c('chr','st','en','id','l','strand','num'))  
merge_gene_promoter <- merge(prmoter,gene_body,by='id')
merge_gene_promoter <- merge_gene_promoter %>%
  mutate(pnum_exp = ceiling((pen-pst)/(pen-pst+en-st)*(pnum+num)),
         num_exp = ceiling((en-st)/(pen-pst+en-st)*(pnum+num)))
merge_gene_promoter$fisher <- '1'
for (i in 1:nrow(merge_gene_promoter)){
  compare <- matrix(c(merge_gene_promoter$pnum[i],merge_gene_promoter$pnum_exp[i],
                      merge_gene_promoter$num[i],merge_gene_promoter$num_exp[i]),nrow=2)
  merge_gene_promoter$fisher[i] <- fisher.test(compare, alternative = "greater")$p.value
}

#--- so you can select the gene you want
# 1. paused gene 1.2 unpaused gene 2. active gene 2.2 unactive gene
all_gene_file <- args$unique_all_gene_bed
all_gene <- fread(all_gene_file,
                  col.names = c('chr','st','en','id','l','strand'))
all_gene <- unique(all_gene)
#-- 1. paused gene
paused_gene <- merge_gene_promoter %>%
  filter(fisher < 0.01) %>%
  select(chr,st,en,id,l,strand,fisher)
paused_gene <- all_gene %>%
  filter(id %in% paused_gene$id)

#-- 2. active gene
active_gene <- ext_gene %>%
  filter(palue < 0.01) %>%
  select(chr,st,en,id,l,strand,palue)
active_gene <- all_gene %>%
  filter(id %in% active_gene$id)
#-- 3. unpaused gene
unpaused_gene <- all_gene %>%
  filter(!(id %in% paused_gene$id))
#-- 4. unactive gene
unactive_gene <- all_gene %>%
  filter(!(id %in% active_gene$id))
#--- class four class 1. paused active 2. paused unactive 
# 3. unpaused active 4. unpaused un active
merge_data <- function(data1,data2){
  merge_d <- merge(data1,data2,by='id')
  merge_d <- merge_d %>%
    mutate(chr=chr.x,st=st.x,en=en.x,l=l.x, strand=strand.x) %>%
    select(chr,st,en,id,l,strand)
  return(merge_d)
}

#--1. paused active
paused_active <- merge_data(paused_gene,active_gene)
#---2. paused unactive
paused_unactive <- merge_data(paused_gene,unactive_gene)
#---3. unpaused active
unpaused_active <- merge_data(unpaused_gene,active_gene)
#---4. unpaused unactive
unpaused_unactive <- merge_data(unpaused_gene,unactive_gene)

save_file <- function(data,name){
  fwrite(data,name,sep = "\t",col.names = F)
}
output_dir<-args$output_dir
setwd(output_dir)
save_file(paused_active,paste0(args$prefix,'paused_active.bed'))
save_file(paused_unactive,paste0(args$prefix,'paused_unactive.bed'))
save_file(unpaused_active,paste0(args$prefix,'unpaused_active.bed'))
save_file(unpaused_unactive,paste0(args$prefix,'unpaused_unactive.bed'))

