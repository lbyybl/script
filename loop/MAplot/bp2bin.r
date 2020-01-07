setwd("/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/MAplot")
anchor_pair <- read.table("pol2_10pet_anchor_uniq.bed",header = FALSE)
#View(anchor_pair)
colnames(anchor_pair) <- c("chr1","start1","end1","chr2","start2","end2")
reso<-5000
around_bin <-0
anchor_pair$bin1_st <- floor((anchor_pair$start1+anchor_pair$end1)/(2*reso)-around_bin)*reso
anchor_pair$bin1_en <- anchor_pair$bin1_st
#anchor_pair$bin1_en <- ceiling((anchor_pair$start1+anchor_pair$end1)/(2*reso)+around_bin)*reso
anchor_pair$bin2_st <- floor((anchor_pair$start2+anchor_pair$end2)/(2*reso)-around_bin)*reso
anchor_pair$bin2_en <- anchor_pair$bin2_st
#anchor_pair$bin2_en <- ceiling((anchor_pair$start2+anchor_pair$end2)/(2*reso)+around_bin)*reso
#library(dplyr)
anchor_pair_bin <- anchor_pair[,c(1,7,8,4,9,10)]
library(data.table)
fwrite(anchor_pair_bin,file="loop_pair_bin.bedpe",sep = "\t",col.names = FALSE)

