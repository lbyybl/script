library(ChIPpeakAnno)
library(UpSetR)
setwd('/WORK/lbyybl/WH/rvb/chip-seq/PolII/unique/diff_peak')


# R105_05 <- 'R105_05_high_deseq2.bed'
# R105_05 <- toGRanges(R105_05,format='BED',header=F)
# 
# R105_1 <- '/DATA/work/lbyybl/wh/ruvb2/chip-seq/Ruvbl1/sample20190509/unique/peak/consistInput_GFP.bed'
# Ruvbl1_peak <- toGRanges(Ruvbl1_bed,format='BED',header=F)
# 
# # h2az <- '/DATA2/work/lbyybl/H2AZ/unique/peak/H2AZ_peaks.narrowPeak'
# # h2az_peak <- toGRanges(h2az,format='BED',header=F)
# 
# znhit1 <- '/DATA2/work/lbyybl/znhit1/unique/peak/znhit1_peaks.narrowPeak'
# znhit1_peak <- toGRanges(znhit1,format='BED',header=F)
# 
# ino80 <- '/DATA2/work/lbyybl/ino80_mES/unique/peak/ino80_merge_rmdup_peaks.narrowPeak'
# ino80_peak <- toGRanges(ino80,format='BED',header=F)
# 
# tip60 <- '/DATA2/work/lbyybl/Tip60_chip/unique/peak/TIP60_rmdup_uniqe_peaks.narrowPeak'
# tip60_peak <- toGRanges(tip60,format='BED',header=F)
# 
# tip60 <- '/DATA2/work/lbyybl/Tip60_chip/unique/peak/TIP60_rmdup_uniqe_peaks.narrowPeak'
# tip60_peak <- toGRanges(tip60,format='BED',header=F)
# 
# MCRS1 <- '/DATA2/work/lbyybl/MCRS1_chip/unique/peak/MCRS1_merge_rmdup_uniqe_peaks.narrowPeak'
# MCRS1_peak <- toGRanges(MCRS1,format='BED',header=F)
for (i in c('R105_05_high','R105_1_high','R105unchange',
           'R1dox_1_high','R1doxr1_high','R1doxunchange',
           'R205dox_05_high','R205dox_dox_high','R205doxunchange')){
  assign(i,toGRanges(paste0(i,'_deseq2.bed'),format='BED',header=F))
}
ov <- findOverlapsOfPeaks(R105_05_high,R105_1_high,R105unchange,
                          R1dox_1_high,R1doxr1_high,R1doxunchange,
                          R205dox_05_high,R205dox_dox_high,R205doxunchange)

pdf('overlap.pdf')
makeVennDiagram(ov)
dev.off()
# names(ov$uniquePeaks)
# names(ov$overlappingPeaks)
# gsub("///",'_',names(ov$peaklist)[11])
# as.data.frame(ov$peaklist$tip60_peak)
# ov$peaklist[10]

#--- for elements shared loops
#setwd('/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/Contact_strength/hiccups')

got_file_list <- function(prefix){
  unset_list<-list()
  i=1
  for (file in system('ls *_high_deseq2.bed *unchange_deseq2.bed',intern = T)){
    #print(file)
    read_data <- read.table(file, col.names = c('chr1','st1','en1'))
    read_data <- read_data %>%
      mutate(loop=str_c(chr1,st1,en1,sep = "-"))
    assign(paste0(prefix,strsplit(file,"\\.")[[1]][1]), read_data$loop)
    unset_list[[i]] <- get(paste0(prefix,strsplit(file,"\\.")[[1]][1]))
    names(unset_list)[i] <- paste0(prefix,strsplit(file,"\\.")[[1]][1])
    i=i+1
  }
  return(unset_list)
}
unset_list <- got_file_list('')
pdf('overlap.pdf',width = 12, height = 6)
upset(fromList(unset_list),nsets=9, order.by = "freq", mb.ratio = c(0.5, 0.5))
dev.off()


for (i in 1:length(names(ov$peaklist))){
  fwrite(as.data.frame(ov$peaklist[i]),gsub("///",'_',names(ov$peaklist)[i]),sep = '\t',col.names = F)
}