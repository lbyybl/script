# find the gain lost stable loop
library(data.table)
library(dplyr)
setwd("/DATA/work/lbyybl/ypjiang/hichip/sample0512/hichipper/pol2_ko_hichipper")

file1 <- "pol2_ko_1.filt.intra.loop_counts.bedpe"
file2 <- "pol2_ko_2.filt.intra.loop_counts.bedpe"
bedpefile1 <- fread(file1, col.names = c("chr1","start1","end1","chr2","start2","end2","name","score"))
bedpefile2 <- fread(file2, col.names = c("chr1","start1","end1","chr2","start2","end2","name","score"))
#--- get gain lost stable

file1region <- bedpefile1 %>%
  mutate(region=paste0(chr1,"-",start1,"-",end1,"-",chr2,"-",start2,"-",end2)) 
file2region <- bedpefile2 %>%
  mutate(region=paste0(chr1,"-",start1,"-",end1,"-",chr2,"-",start2,"-",end2))

file1_stable <- file1region %>%
  filter(region %in% file2region$region) %>%
  select(chr1,start1,end1,chr2,start2,end2,name,score)
file1_spe <- file1region %>%
  filter(!(region %in% file2region$region))%>%
  select(chr1,start1,end1,chr2,start2,end2,name,score)

file2_stable <- file2region %>%
  filter(region %in% file1region$region)%>%
  select(chr1,start1,end1,chr2,start2,end2,name,score)
file2_spe <- file2region %>%
  filter(!(region %in% file1region$region))%>%
  select(chr1,start1,end1,chr2,start2,end2,name,score)

#--- write to the file
lengthfile1 <- length(stringr::str_split(file1,'\\.')[[1]])
outputfilname <- function(file1, str){
  stringr::str_c(c(stringr::str_split(file1,'\\.')[[1]][-lengthfile1],str,stringr::str_split(file1,'\\.')[[1]][lengthfile1]),collapse = ".")
}

fwrite(file1_stable,file=outputfilname(file1,"stable"),sep="\t",col.names = FALSE)
fwrite(file1_spe,file=outputfilname(file1,"spe"),sep="\t",col.names = FALSE)
fwrite(file2_stable,file=outputfilname(file2,"stable"),sep="\t",col.names = FALSE)
fwrite(file2_spe,file=outputfilname(file2,"spe"),sep="\t",col.names = FALSE)

