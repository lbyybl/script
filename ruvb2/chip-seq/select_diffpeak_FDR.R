#-----------------------------------
# this script is used to clasified the protein to different class
setwd('/WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl/Pol2/diff_peak')
library(data.table)
library(dplyr)
#-------------------------------
# got different region for Ruvbl1/2 peak
setwd('/WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl/Pol2/diff_peak')
file <- 'R1dox_deseq2.txt'
name <- 'R1dox'
filter_fdr <- function(file,name){
  data <- fread(file,header = T)
  #data2 <- fread('RUVB2_POL2_deseq2_sig.bed',header = F)
  names(data)
  data <- as.data.frame(data)
  unchange <- data %>%
    filter(FDR >= 0.05 | (FDR < 0.05 & abs(Fold) < log(2))) %>%
    dplyr::select(seqnames,start,end)
  #unchange[,.(seqnames,start,end)]
  down <- data %>% 
    filter(FDR <0.05 & Fold > log(2)) %>%
    dplyr::select(seqnames,start,end)
  up <- data %>% 
    filter(FDR < 0.05 & Fold < -log(2)) %>%
    dplyr::select(seqnames,start,end)
  nrow(unchange)
  nrow(down)
  nrow(up)
  fwrite(unchange,paste0(name,'unchange_deseq2.txt'),sep = '\t',col.names = F)
  fwrite(down,paste0(name,'r2_high_deseq2.txt'),sep = '\t',col.names = F)
  fwrite(up,paste0(name,'r1_high_deseq2.txt'),sep = '\t',col.names = F)
  
}

filter_fdr('R05dox_deseq2.txt','R05dox')
filter_fdr('R105_deseq2.txt','R105')
filter_fdr('R1dox_deseq2.txt','R1dox')
