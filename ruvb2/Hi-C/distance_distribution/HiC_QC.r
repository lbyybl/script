library(HiTC)
# matrix <- '/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/merge/ko_out/hic_results/matrix/ko/iced/100000/ko_100000_iced.matrix'
# bed <- '/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/merge/ko_out/hic_results/matrix/ko/raw/100000/ko_100000_abs.bed'
setwd('/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/loop/distance_distribute/100K')

read_data <- function(prefix){
  matrix <- paste0(prefix,'_iced.matrix')
  bed <- paste0(prefix,'_abs.bed')
  hic_data <- importC(matrix,bed,bed)
  return(hic_data)
}
R2_ko_rep1 <- read_data('R2-IAA-Hi-C1_100000')

pdf('R2_ko_rep1.pdf')
par(mfrow=c(3,2))
CQC(R2_ko_rep1, winsize = 1e+06, dev.new=FALSE, hist.dist=T)
dev.off()

R2_ko_rep2 <- read_data('R2-IAA-Hi-C2_100000')

pdf('R2_ko_rep2.pdf')
par(mfrow=c(3,2))
CQC(R2_ko_rep2, winsize = 1e+06, dev.new=FALSE, hist.dist=T)
dev.off()

R2_wt_rep1 <- read_data('R2-NT-Hi-C1_100000')

pdf('R2_wt_rep1.pdf')
par(mfrow=c(3,2))
CQC(R2_wt_rep1, winsize = 1e+06, dev.new=FALSE, hist.dist=T)
dev.off()

R2_wt_rep2 <- read_data('R2-NT-Hi-C2_100000')

pdf('R2_wt_rep2.pdf')
par(mfrow=c(3,2))
CQC(R2_wt_rep2, winsize = 1e+06, dev.new=FALSE, hist.dist=T)
dev.off()

# hic_data <- importC(matrix,bed,bed)
# par(mfrow=c(2,2))
# CQC(hic_data, winsize = 1e+06, dev.new=FALSE, hist.dist=FALSE)
