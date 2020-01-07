setwd('/WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl/Pol2/graph/travel')
R1 <- fread('R1_tss_reads.bed')
R2 <- fread('R2_tss_reads.bed')
merge <- merge(R1,R2,by=c('V1','V2','V3','V4','V5','V6'))
merge <- merge %>%
  mutate(sum=V7.x+V7.y) %>%
  dplyr::select(V1,V2,V3,V4,V5,V6,sum)
fwrite(merge,'Ruvbl_tss_counts.bed',sep = '\t',col.names = F)
merge <- merge %>%
  dplyr::filter(V1 %in% paste0('chr',c(1:19,'X','Y')))
tr <- fread('Pol2_tr.bed',col.names = c('chr','st','en','ratiodaox','ratio05','ratio1','gene','strand'))
tr_merge <- merge(tr,merge,by.x='gene',by.y='V4')
tr_merge <- tr_merge %>%
  dplyr::select(chr,st,en,ratiodaox,ratio05,ratio1,sum) %>%
  arrange(-sum)
nrow(unique(tr_merge))
top1000 <- tr_merge[1:1000,1:6]
tail1000 <- tr_merge[(nrow(tr_merge)-999):nrow(tr_merge),1:6]
fwrite(top1000,'Ruvbl_binding_gene_tr.bed',sep = '\t',col.names = F)
fwrite(tail1000,'Ruvbl_unbinding_gene_tr.bed',sep = '\t',col.names = F)

#--- select prmoter for meta analysis
setwd('/WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl/Pol2/graph/travel')
R1 <- fread('R1_tss_reads.bed')
R2 <- fread('R2_tss_reads.bed')
merge <- merge(R1,R2,by=c('V1','V2','V3','V4','V5','V6'))
merge <- merge %>%
  mutate(sum=V7.x+V7.y) %>%
  dplyr::select(V1,V2,V3,V4,V5,V6,sum)
#fwrite(merge,'Ruvbl_tss_counts.bed',sep = '\t',col.names = F)
merge <- merge %>%
  dplyr::filter(V1 %in% paste0('chr',c(1:19,'X','Y')))
tr <- fread('Pol2_tr.bed',col.names = c('chr','st','en','ratiodaox','ratio05','ratio1','gene','strand'))
setwd('/WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl/Pol2/graph/TSS_profile')
tr_merge <- merge(tr,merge,by.x='gene',by.y='V4')
tr_merge <- tr_merge %>%
  dplyr::select(chr,V2,V3,gene,V5,strand,sum) %>%
  arrange(-sum)
tr_merge <- tr_merge %>%
  dplyr::filter(chr %in% paste0('chr',c(1:19,'X','Y')))
nrow(unique(tr_merge))
top1000 <- tr_merge[1:1000,1:6]
tail1000 <- tr_merge[(nrow(tr_merge)-999):nrow(tr_merge),1:6]
fwrite(tr_merge,'Ruvbl_tss_rank.bed',sep = '\t',col.names = F)
fwrite(top1000,'Ruvbl_binding_TSS.bed',sep = '\t',col.names = F)
fwrite(tail1000,'Ruvbl_unbinding_TSS.bed',sep = '\t',col.names = F)
