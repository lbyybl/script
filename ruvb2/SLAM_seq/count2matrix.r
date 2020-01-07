#--- it's used to extract the result run by samtools flagstat
setwd('/DATA/work/lbyybl/wh/ruvb2/SLAM/slam_output/count_2TC')
library(readr)
# #suppressMessages(suppressWarnings(flag_file <- read_tsv('L13-R2-IAA-IAA-SLAM_FKDL171549253-1A_1.fq_slamdunk_mapped_filtered_tcount.tsv')))
# R2_IAA_IAA_1 <- read_tsv('L13-R2-IAA-IAA-SLAM_FKDL171549253-1A_1.fq_slamdunk_mapped_filtered_tcount.tsv',skip = 2)
# #View(R2_IAA_IAA_1)
# R2_IAA_IAA_1 <- R2_IAA_IAA_1[,c(4,13)]
# R2_IAA_IAA_2 <- read_tsv('L13-R2-IAA-IAA-SLAM_FKDL171549253-1A_2.fq_slamdunk_mapped_filtered_tcount.tsv',skip = 2)
# #View(R2_IAA_IAA_2)
# R2_IAA_IAA_2 <- R2_IAA_IAA_2[,c(4,13)]
# R2_IAA_IAA <- merge(R2_IAA_IAA_1,R2_IAA_IAA_2,by='Name',all=T)
# R2_IAA_IAA$TcReadCount <- R2_IAA_IAA$TcReadCount.x + R2_IAA_IAA$TcReadCount.y
# R2_IAA_IAA <- R2_IAA_IAA[,c(1,4)]
#View(R2_IAA_IAA)
readdata <- function(name){
  data1 <- read_tsv(paste0(name,'_1.fq_slamdunk_mapped_filtered_tcount.tsv'),skip = 2)
  data1 <- data1[,c(4,13)]
  data2 <- read_tsv(paste0(name,'_2.fq_slamdunk_mapped_filtered_tcount.tsv'),skip = 2)
  data2 <- data2[,c(4,13)]
  data_merge <- merge(data1,data2,by='Name',all=T)
  data_merge$TcReadCount <- data_merge$TcReadCount.x + data_merge$TcReadCount.y
  data_merge <- data_merge[,c(1,4)]
  print(paste0('data1 => ',sum(data1$TcReadCount),'    data2 => ',sum(data2$TcReadCount),'    data_merge => ',sum(data_merge$TcReadCount)))
  return(data_merge)
}
R2_IAAplusIAA <- readdata('L12-R2-IAA+IAA-SLAM_FKDL171549252-1A')
R2_WTplusIAA <- readdata('L10-R2-NT+IAA-SLAM_FKDL171549250-1A')
all_data_matrix <- merge(R2_WTplusIAA,R2_IAAplusIAA,by='Name',all=T)
colnames(all_data_matrix) <- c('gene','R2_WTplusIAA','R2_IAAplusIAA')

nrow(all_data_matrix)
nrow(all_data_matrix[all_data_matrix$R2_WTplusIAA + all_data_matrix$R2_IAAplusIAA>10,])
View(all_data_matrix)
sum(all_data_matrix$R2_WTplusIAA)
sum(all_data_matrix$R2_IAAplusIAA)

R2_IAAnoIAA <- readdata('L13-R2-IAA-IAA-SLAM_FKDL171549253-1A')
sum(R2_IAAnoIAA$TcReadCount)
R2_WTnoIAA <- readdata('L11-R2-NT-IAA-SLAM_FKDL171549251-1A')
sum(R2_WTnoIAA$TcReadCount)