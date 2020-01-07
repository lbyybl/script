#---筛选+1处都有核小体的基因；与可重复的至少有一个核小体的重复；
setwd('/WORK/lbyybl/WH/rvb/cor_uniq/MNase/graph/TSS_profile/deeptools')
data <- fread('plus1nuclesome_summary_peak.tab',skip = 1, col.names = c('chr','st','en',
                                                                        'h05rep1','h1rep1','doxrep1',
                                                                        'h05rep2','h1rep2','doxrep2'),
              stringsAsFactors = F)
data <- as.data.frame(data)
data <- data %>%
  arrange(-h1rep1)
head(data)
data[which(data[,7]==max(data[,4:9])),]
med <- median(c(data[,c(4)],data[,5],data[,6],data[,7],data[,8],data[,9]))
for (i in 4:ncol(data)){
  data[data[,i]>med,i] <- 1
  data[data[,i]<=0,i] <- 0
}

data_all <- data %>%
  filter(h05rep1+h05rep2+h1rep1+h1rep2+doxrep1+doxrep2==6) %>%
  select(chr,st,en)

data_any <- data %>%
  dplyr::filter((h05rep1==h05rep2) & (h1rep1==h1rep2) & (doxrep1==doxrep2) & 
           (h05rep1+h05rep2+h1rep1+h1rep2+doxrep1+doxrep2 < 6) &
           (h05rep1+h05rep2+h1rep1+h1rep2+doxrep1+doxrep2 > 0 )) %>%
  select(chr,st,en)

plus1nuc <- fread('plus1nuclesome.bed',
                  col.names = c('chr','st','en','name','sc','strand')) %>%
  select(chr,st,en,name)

gene <- fread('/WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl/Pol2/graph/travel/Pol2_tr.bed',
              col.names = c('chr','st','en','no1','no2','no3','name','strand')) %>%
  select(name,strand)
plus1nuc <- merge(plus1nuc,gene,by='name')
data <- fread('plus1nuclesome_summary_peak.tab',skip = 1, col.names = c('chr','st','en',
                                                                        'h05rep1','h1rep1','doxrep1',
                                                                        'h05rep2','h1rep2','doxrep2'),
              stringsAsFactors = F)

data_all <- merge(data_all,plus1nuc,by=c('chr','st','en'))
data_all <- merge(data_all,data,by=c('chr','st','en'))
data_all <- unique(data_all)
data_any <- merge(data_any,plus1nuc,by=c('chr','st','en'))
data_any <- merge(data_any,data,by=c('chr','st','en'))
data_any <- unique(data_any)
fwrite(data_all,'all_have_plus1.bed',sep = '\t',col.names = F)
fwrite(data_any,'any_have_plus1.bed',sep = '\t',col.names = F)
