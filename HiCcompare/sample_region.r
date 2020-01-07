#--- this file is used to produce the sampleed region
# 1. 使用用R写一个程序：a.该程序可以生成每条染色体相应分辨率的互做对的文件；
# b. 可以对生成的文件按照距离分布分成若干份；
# c. 根据距离分布文件提出相互作用对；
#------------------------------------------------------------------
library(dplyr)
library(data.table)
library(purrr)
options(scipen = 200)
options(stringsAsFactors = FALSE)
# chr19   61431566
#chr19_len <- 61431566
#--- produce bin according to the chr length
#chr_len <- chr19_len
#chr="chr19"
reso <- 50000

gotbin_pair <- function(chr_len,chr,reso){
  bin=(0:floor(chr_len/reso))*reso
  bin_pair=data.frame("chr"=chr,
                      "bin1"=rep(bin,each=length(bin)),
                      "bin2"=rep(bin,length(bin)))
  #--- got the distribute
  
  bin_pair$distrib <- abs((bin_pair$bin2-bin_pair$bin1)/reso)
  return(bin_pair)
}

chr_size_file <- '/home/boyuanli/tools/hic-pro/HiC-Pro_2.9.0/annotation/chrom_mm10.sizes'

chr_size <- fread(chr_size_file,col.names = c("chr","length"))
chr_size <- chr_size[-which(chr_size$chr=='chrM'),]
chr_size <- chr_size[-which(chr_size$chr=='chrMT'),]
chr_bin_pair <- data.frame("chr"=NULL,
                           "bin1"=NULL,
                           "bin2"=NULL,
                           "distrib"=NULL)
for (i in 1:nrow(chr_size)){
  assign(paste0(chr_size$chr[i],"_binpair"), gotbin_pair(chr_size$length[i],chr_size$chr[i],reso))
  chr_bin_pair <- bind_rows(chr_bin_pair,get(paste0(chr_size$chr[i],"_binpair")))
}

#fwrite(chr_bin_pair,'all_pair.bed',sep = '\t',col.names = T)


#--------------------------------------------------------------------
#
#--------------------------------------------------------------------
setwd('/home/boyuanli/bashscript/bin/HiCcompare/sample')
chr_bin_pair <- fread('all_pair.bed',header = T)
#sig_dis <- fread('signignificant_interaction.bed')

#--- got the proportion of different distance pets 

YY1 <- '/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/YY1/all_signifent_YY1_sub.allValidPairs_.diff'
CTCF <- '/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/CTCF/all_signifent_CTCF_sub.allValidPairs_.diff'
pol1 <- '/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/pol1/all_signifent_pol1_75mi_sub.allValidPairs_.diff'
#pol1vspol2 <- '/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/pol1vspol2/all_signifent_pol2_wt_75mi_sub.allValidPairs_.diff'
#pol1vspol3 <- '/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/pol1vspol3/all_signifent_pol3_wt_75mi_sub.allValidPairs_.diff'
pol2 <- '/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/pol2/all_signifent_pol2_75mi_sub.allValidPairs_.diff'
#pol2vspol3 <- '/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/pol2vspol3/all_signifent_pol3_wt_75mi_sub.allValidPairs_.diff'
pol3 <- '/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/pol3/all_signifent_pol3_75mi_sub.allValidPairs_.diff'

#file <- YY1
readfile <- function(file){
  data <- fread(file) 
  data <- dplyr::select(data,c(V1,V2,V5))
  colnames(data) <- c('chr','bin1','bin2')
  data$distrib <- abs((data$bin2-data$bin1)/reso)
  return(data)
}
sig_pair <- data.frame("chr"=NULL,
                           "bin1"=NULL,
                           "bin2"=NULL,
                           "distrib"=NULL,
                       "name"=NULL)
for (i in c('YY1','CTCF','pol1','pol2','pol3')){ #,'pol1vspol2','pol1vspol3','pol2vspol3')){
  assign(paste0(i,"_binpair"), readfile(get(i)))
  #get(paste0(i,"_binpair"))$name <- i
  sig_pair <- bind_rows(sig_pair,get(paste0(i,"_binpair")))
}

sig_dis <- sig_pair %>%
  group_by(distrib) %>%
  summarise(num=round(n()/8),pro=n()/nrow(sig_pair))
#fwrite(sig_dis,'signignificant_interaction.bed',sep = '\t')


#--- sample the pirs with the same distance distribution
sample_pair <- data.frame("chr"=NULL,
                           "bin1"=NULL,
                           "bin2"=NULL,
                           "distrib"=NULL)
# for (i in sig_dis$distrib){
#   sample <- sample_n(chr_bin_pair[chr_bin_pair$distrib==i,],as.numeric(sig_dis[sig_dis$distrib==i,2]))
#   sample_pair <- bind_rows(sample_pair,sample)
# }

sample_num <- function(num){
  return(sample_n(chr_bin_pair[chr_bin_pair$distrib==num,],as.numeric(sig_dis[sig_dis$distrib==num,2])))
  #sample_pair <- bind_rows(sample_pair,sample)
  #return(sample)
}
#sample_n(chr_bin_pair[chr_bin_pair$distrib==i,],as.numeric(sig_dis[sig_dis$distrib==i,2]))
#apply(data.frame('data'=sig_dis$distrib),sample_num)
sig_dis <- sig_dis[-which(sig_dis$num==0),]
format_data <- function(data){
  data$bin1 <- format(data$bin1, scientific = F)
  data$bin2 <- format(data$bin2, scientific = F)
  return(data)
}
sample_one_pair <- function(name){
  assign(paste0(name,'sample_pair'),sig_dis$distrib %>%
           map_dfr(sample_num))
  assign(paste0(name,'sample_pair'),format_data(get(paste0(name,'sample_pair'))))
  fwrite(get(paste0(name,'sample_pair')),paste0(name,'sample_pair'),sep='\t',col.names = F)
}

# paste0('sample',1:100) %>%
#   map_chr(sample_one_pair)

for (i in paste0('sample',1:100)){
  sample_one_pair(i)
}
sample_pair <- sig_dis$distrib %>%
              map_dfr(sample_num)

#--- random region
setwd('/home/boyuanli/bashscript/bin/HiCcompare/sample/sample_pair_result')

sample100 <- 'sample100sample_pair/sample100sample_pair.txt'
data <- fread(sample100,col.names = c('GG','GN','NN'))
