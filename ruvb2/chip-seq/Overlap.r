library(ChIPpeakAnno)
setwd('/WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl')
bed <- system.file("extdata", "MACS_output.bed", package="ChIPpeakAnno")
gr1 <- toGRanges(bed, format="BED", header=FALSE) 
## one can also try import from rtracklayer
gff <- system.file("extdata", "GFF_peaks.gff", package="ChIPpeakAnno")
gr2 <- toGRanges(gff, format="GFF", header=FALSE, skip=3)
## must keep the class exactly same as gr1$score, i.e., numeric.
gr2$score <- as.numeric(gr2$score) 
ol <- findOverlapsOfPeaks(gr1, gr2)
## add metadata (mean of score) to the overlapping peaks
ol <- addMetadata(ol, colNames="score", FUN=mean) 
ol$peaklist[["gr1///gr2"]][1:2]
makeVennDiagram(ol)

Ruvbl2_bed <- '/WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl/overlap/Ruvbl2_merge.bed'
Ruvbl2_peak <- toGRanges(Ruvbl2_bed,format='BED',header=F)

Ruvbl1_bed <- '/WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl/overlap/Ruvbl1_merge.bed'
Ruvbl1_peak <- toGRanges(Ruvbl1_bed,format='BED',header=F)

# h2az <- '/DATA2/work/lbyybl/H2AZ/unique/peak/H2AZ_peaks.narrowPeak'
# h2az_peak <- toGRanges(h2az,format='BED',header=F)

znhit1 <- 'znhit1_peak.bed'
znhit1_peak <- toGRanges(znhit1,format='BED',header=F)

ino80 <- 'ino80_peak.bed'
ino80_peak <- toGRanges(ino80,format='BED',header=F)

tip60 <- 'tip60_peak.bed'
tip60_peak <- toGRanges(tip60,format='BED',header=F)

# tip60 <- '/DATA2/work/lbyybl/Tip60_chip/unique/peak/TIP60_rmdup_uniqe_peaks.narrowPeak'
# tip60_peak <- toGRanges(tip60,format='BED',header=F)

# MCRS1 <- '/DATA2/work/lbyybl/MCRS1_chip/unique/peak/MCRS1_merge_rmdup_uniqe_peaks.narrowPeak'
# MCRS1_peak <- toGRanges(MCRS1,format='BED',header=F)

ov <- findOverlapsOfPeaks(Ruvbl1_peak,Ruvbl2_peak,znhit1_peak,
                          ino80_peak,tip60_peak,connectedPeaks = 'merge')
pdf('overlap.pdf')
makeVennDiagram(ov)
dev.off()
# names(ov$uniquePeaks)
# names(ov$overlappingPeaks)
# gsub("///",'_',names(ov$peaklist)[11])
# as.data.frame(ov$peaklist$tip60_peak)
# ov$peaklist[10]

for (i in 1:length(names(ov$peaklist))){
  fwrite(as.data.frame(ov$peaklist[i]),gsub("///",'_',names(ov$peaklist)[i]),sep = '\t',col.names = F)
}
