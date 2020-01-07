#!/usr/bin/env R
# Sun Jan 27 11:16:02 2019
# Boyuan_Li

# there are many gene share the same TSS or TSS and TTS; this script is used to remove the
# gene share the same TSS and TTS; if the script share the same TSS, I will choose the shorter
# gene and remove the longer gene;
output_dir <- './'
command=matrix(c(
	"input_file", "i", "1", "character", 'need the  format chr st en id score strand',
	"output_file", "o", "1", "character", 'the output file name',
	"output_dir", "d", "2", "character" ,'defult is ./; output_file will place in it'
	),byrow=T,ncol=5)


args=getopt::getopt(command)


if (is.null(args$input_file) || is.null(args$output_file) || is.null(args$output_dir)) {
	cat(paste(getopt::getopt(command, usage = T), "\n"))
	q()
}

library(data.table)
library(dplyr)

promoter <- fread(args$input_file,
                  col.names = c('chr','st','en','id','score','strand'))
promoter <- promoter %>% 
  arrange(chr,st,en)

n <- nrow(promoter)
for (i in 1:n){
  if (i+1 <= n){
    if (sum(promoter[i,c(1:3,6)]==promoter[i+1,c(1:3,6)])==4){
      promoter[i+1,4] <- promoter[i,4]
      
    }
  }
  
}

uniq_progene <- unique(promoter)
k <- nrow(uniq_progene)
uniq_progene <- uniq_progene %>%
  arrange(chr,st,en)

for (i in 1:(k-1)){
  if (uniq_progene[i,6]=="+" & uniq_progene[i+1,6]=="+"){
    if (sum(uniq_progene[i,c(1,2,6)]==uniq_progene[i+1,c(1,2,6)])==3){
      if (uniq_progene[i,3]<=uniq_progene[i+1,3]){
        uniq_progene[i+1,c(3,4)]<-uniq_progene[i,c(3,4)]
      }
      else {
        uniq_progene[i,c(3,4)]<-uniq_progene[i+1,c(3,4)] 
      }
      
    }
  }
}  

uniq_progene <- unique(uniq_progene)
k <- nrow(uniq_progene)
uniq_progene <- uniq_progene %>%
  arrange(chr,st)
for (i in 1:(k-1)){
  if (uniq_progene[i,6]=="-" & uniq_progene[i+1,6]=="-"){
    if (sum(uniq_progene[i,c(1,3,6)]==uniq_progene[i+1,c(1,3,6)])==3){
      if (uniq_progene[i,2]>=uniq_progene[i+1,2]){
        uniq_progene[i+1,c(2,4)]<-uniq_progene[i,c(2,4)]
      }
      else {
        uniq_progene[i,c(2,4)]<-uniq_progene[i+1,c(2,4)] 
      }
      
    }
  }
}
uniq_progene <- unique(uniq_progene)
output_dir<-args$output_dir
setwd(output_dir)
fwrite(uniq_progene,args$output_file,sep = "\t",col.names = F)