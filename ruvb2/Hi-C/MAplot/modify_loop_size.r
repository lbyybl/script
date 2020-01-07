file <- '/DATA2/work/lkh/pol2_1226/loop/wt_loop/merged_loops.bedpe'
library(data.table)
library(dplyr)

loop_file <- fread(file)
loop_file <- loop_file[,1:6]

for (i in 1:nrow(loop_file)) {
  dis <- loop_file[i,3]-loop_file[i,2]
  if (dis != 25000){
    loop_file[i,3] <- (floor(loop_file[i,3]/25000)+1)*25000
    loop_file[i,2] <- floor(loop_file[i,2]/25000)*25000
  } 
  dis2 <- loop_file[i,6]-loop_file[i,5]
  if (dis2 != 25000){
    loop_file[i,6] <- (floor(loop_file[i,6]/25000)+1)*25000
    loop_file[i,5] <- floor(loop_file[i,5]/25000)*25000
  } 
}
colnames(loop_file) <- c('chr1','st1','en1','chr2','st2','en2')
loop_file$chr1 <- paste0('chr',loop_file$chr1)
loop_file$chr2 <- paste0('chr',loop_file$chr2)
# for (i in 1:nrow(loop_file)){
#   dis <- loop_file[i,3]-loop_file[i,2]
#   if (dis==5000){
#     loop_file[i,3] <- loop_file[i,3]+10000
#     loop_file[i,2] <- loop_file[i,2]-10000
#   }else if (dis==10000) {
#     loop_file[i,3] <- loop_file[i,3]+5000
#     loop_file[i,2] <- loop_file[i,2]-10000 
#   }
#   dis2 <- loop_file[i,6]-loop_file[i,5]
#   if (dis2==5000){
#     loop_file[i,6] <- loop_file[i,6]+10000
#     loop_file[i,5] <- loop_file[i,5]-10000
#   }else if (dis2==10000) {
#     loop_file[i,6] <- loop_file[i,6]+5000
#     loop_file[i,5] <- loop_file[i,5]-10000 
#   }
# }
setwd('/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/loop')
fwrite(loop_file,'Pol2_wt_loop.bedpe',col.names = F,sep = '\t')

unique(sort(loop_file[,3]-loop_file[,2]))