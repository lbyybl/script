#---分四类基因
library(data.table)
library(dplyr)
library(stringr)
active_gene_file <- '/DATA2/work/lbyybl/pol2_pro_seq2015/mapping/paused_gene/result/class/define_with_pol2chip/active.bed'
unactive_gene_file <- '/DATA2/work/lbyybl/pol2_pro_seq2015/mapping/paused_gene/result/class/define_with_pol2chip/unactive.bed'

travelratio_file <- '/DATA2/work/lbyybl/pol2_pro_seq2015/mapping/paused_gene/result/class/define_with_pol2chip/pol2_chip_TR.mat'

active_gene <- fread(active_gene_file)
unactive_gene <- fread(unactive_gene_file)
travelratio <- fread(travelratio_file)
travelratio <- travelratio[,3:4]
travelratio$TSS <- as.numeric(strsplit2(travelratio$V4, ",")[,1])
travelratio$GB <- as.numeric(strsplit2(travelratio$V4, ",")[,2])
travelratio <- travelratio %>%
  mutate(bi=(TSS+1)/(GB+1))
travelratio$class <- 'no'
active_gene$class <- 'active'
unactive_gene$class <- 'unactive'
gene_acorunac <- bind_rows(active_gene,unactive_gene)
for (i in 1:nrow(travelratio)){
  if (travelratio$bi[i] >=2){
    travelratio$class[i] <- 'paused'
  }
  else{
    travelratio$class[i] <- 'unpaused'
  }
}

tow_class_merge <- merge(gene_acorunac,travelratio,by.x='V4',by.y='V3')

tow_class_merge <- tow_class_merge %>%
  mutate(class=paste0(class.x,"-",class.y)) 
tow_class_merge <- tow_class_merge[,c(2,3,4,1,11,6,13)]

filterdata <- function(data,type){
  data2 <- data %>%
    filter(class==type)
  return(data2)
}
active_paused <- filterdata(tow_class_merge,'active-paused')
unactive_paused <- filterdata(tow_class_merge,'unactive-paused')
unactive_unpaused <- filterdata(tow_class_merge,'unactive-unpaused')
active_unpaused <- filterdata(tow_class_merge,'active-unpaused')

fwrite(active_paused,'paused_active.bed',col.names = F, quote = F,sep = "\t")
fwrite(unactive_paused,'paused_unactive.bed',col.names = F, quote = F,sep = "\t")
fwrite(unactive_unpaused,'unpaused_unactive.bed',col.names = F, quote = F,sep = "\t")
fwrite(active_unpaused,'unpaused_active.bed',col.names = F, quote = F,sep = "\t")


