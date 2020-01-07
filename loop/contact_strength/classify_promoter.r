#--- read bam file
library(GenomicRanges)
library(GenomicAlignments)
setwd('/DATA2/work/lbyybl/pol2_pro_seq2015/mapping')
plus_gro_seq <- readGAlignments('SRR1991266_plus_sort.bam')
minus_gro_seq <- readGAlignments('SRR1991266_minus_sort.bam')
frame_gro_plus <- as.data.frame(plus_gro_seq)
frame_gro_minus <- as.data.frame(minus_gro_seq)
frame_gro_plus <- frame_gro_plus %>%
  select(seqnames,start,end,width,njunc,strand)
frame_gro_minus <- frame_gro_minus %>%
  select(seqnames,start,end,width,njunc,strand)
frame_gro <- bind_rows(frame_gro_minus,frame_gro_plus)
fwrite(frame_gro,'gro_seq_mapping.bed', sep = "\t",col.names = F)

#---------------------------------------------------
# write a script to cauculate the activety gene and paused gene
#---------------------------------------------------
# gene activity, need lamada and N/L
# lamada need an no duplicate gene_bed file to calculate the gene_length
# or you can used the merged gene_bed file to calculate it;
# N/L you need a file to contain how many reads fall into one gene region;
# for this file you can used the bedtools intersect to get the reads in 
# every region, for speed you can used +- seperately;
# and then make a possion test 

#-------for paused gene
# stastic the reads num in promoter and in gene body and cautulate the 
# thory num in promoter and gene body according to the length;
# and give fisher stastic; but used the promoter proximal peaks is very
# difficulty, so i want used the all promoter to do this and ave(p) >
# ave(g) and significant is called paused
# summary, 1. a bed file contain promoter and gene body;
# 2. a bed file contail merge gene to calculate lambda (notice extend 1k up);
# 3. a 
#--- produce promoter file step 5 bp 
setwd('/DATA2/work/lbyybl/pol2_pro_seq2015/mapping/paused_gene/bash')
# promoter <- fread('promoter_final_1000.bed',
#                   col.names = c('chr','st','en','id','score','strand'))
# nrow(unique(promoter))
# 
# for (i in 1:nrow(promoter)){
#   if (sum(promoter[i,c(1:3,6)]==promoter[i+1,c(1:3,6)])==4){
#     promoter[i+1,4] <- promoter[i,4] 
#   }
# }
# promoter <- unique(promoter)
promoter <- fread('gene_len3k.bed',
                  col.names = c('chr','st','en','id','score','strand'))
promoter <- promoter %>% 
  arrange(chr,st,en)

n <- nrow(promoter)
for (i in 1:n){
  if (i+1 <= n){
    if (sum(promoter[i,c(1:3,6)]==promoter[i+1,c(1:3,6)])==4){
      promoter[i+1,4] <- promoter[i,4]

    }
  }

}

uniq_progene <- unique(promoter)
k <- nrow(uniq_progene)
uniq_progene <- uniq_progene %>%
  arrange(chr,st,en)

for (i in 1:k){
  if (i+1 <= k){
    if (sum(uniq_progene[i,c(1,2,6)]==uniq_progene[i+1,c(1,2,6)])==3){
      if (uniq_progene[i,3]<=uniq_progene[i+1,3]){
        uniq_progene[i+1,c(3,4)]<-uniq_progene[i,c(3,4)]
      }
      else {
        uniq_progene[i,c(3,4)]<-uniq_progene[i+1,c(3,4)] 
      }
      
    }
  }
}  
  
uniq_progene <- unique(uniq_progene)
k <- nrow(uniq_progene)
uniq_progene <- uniq_progene %>%
  arrange(chr,st)
  for (i in 1:k){
    if (i+1 <= k){
      if (sum(uniq_progene[i,c(1,3,6)]==uniq_progene[i+1,c(1,3,6)])==3){
        if (uniq_progene[i,2]>=uniq_progene[i+1,2]){
          uniq_progene[i+1,c(2,4)]<-uniq_progene[i,c(2,4)]
        }
        else {
          uniq_progene[i,c(2,4)]<-uniq_progene[i+1,c(2,4)] 
        }
        
      }
    }
  }
uniq_progene <- unique(uniq_progene)
fwrite(uniq_progene,'unique_short_gene.bed',sep = "\t")

#--- calculate mappable length
setwd('/DATA2/work/lbyybl/pol2_pro_seq2015/mapping/paused_gene/result/reads')
mappalbe_bed <- fread('mappable_base.bed',col.names = c('chr','st','en','strand'))

sum(mappalbe_bed$en-mappalbe_bed$st)
map_reads <- fread('../gro_seq_mapping.bed')
sum(map_reads$V4)

#--- calculate lamada
  #--- rm fist and last line for each chr
  setwd("/DATA2/work/lbyybl/pol2_pro_seq2015/mapping/paused_gene/result/backgroung")
  non_gene_region <- fread('non_gene_region.bed')  
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
  fwrite(non_gene_region,'non_gene_region.bed',sep = '\t',col.names = F)
  #--- calculate lamada
  non_gene_reads_num <- fread('non_gene.bed',
                              col.names = c('chr','st','en','id','l','strand','num'))
  
  lamada <- sum(non_gene_reads_num$num)/sum(non_gene_reads_num$en-non_gene_reads_num$st)  
  
#--- calculate the sig gene and sig for gene
  setwd('/DATA2/work/lbyybl/pol2_pro_seq2015/mapping/paused_gene/result/reads_num')
  gene_file <- 'genebody_final_-1000_reads_num.bed'
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
  setwd('/DATA2/work/lbyybl/pol2_pro_seq2015/mapping/paused_gene/result/reads')
  promoter_file <- 'promoter_final_1000_stis.bed'
  #prmoter_sti <- sig_calcu(promoter_file)
  prmoter <- fread(promoter_file,
                   col.names = c('pchr','pst','pen','id','pl','pstrand','pnum'))  
  #--- gene short 1k 
  gene_body_file <- '/DATA2/work/lbyybl/pol2_pro_seq2015/mapping/paused_gene/result/reads_num/genebody_final_1000_reads_num.bed'
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
  all_gene_file <- '../../bash/unique_short_gene.bed'
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
  output_dir<-'/DATA2/work/lbyybl/pol2_pro_seq2015/mapping/paused_gene/result/class'
  setwd(output_dir)
  save_file(paused_active,'paused_active.bed')
  save_file(paused_unactive,'paused_unactive.bed')
  save_file(unpaused_active,'unpaused_active.bed')
  save_file(unpaused_unactive,'unpaused_unactive.bed')
  # paused_unactive <- merge(paused_gene,unactive_gene,by='id')
  # paused_unactive <- paused_unactive %>%
  #   mutate(chr=chr.x,st=st.x,en=en.x,l=l.x, strand=strand.x) %>%
  #   select(chr,st,en,id,l,strand)
  # 