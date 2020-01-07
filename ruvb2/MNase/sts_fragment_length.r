setwd('/WORK/lbyybl/WH/rvb/ATAC/sample20191014/test')
ratio <- fread('ratio.txt')
ratio_chr <- ratio %>%
  filter(seqnames %in% paste0('chr',c(1:19,'X','Y')))
sum(ratio_chr$mapped)/sum(ratio$mapped)
sum(ratio$mapped[22])/sum(ratio$mapped)
A=c(1,2,1,2,1,2)
fft(A)
plot(A)
line(A)
a <- c(1:300)
y<-abs(cos(pi/150*a))
plot(a,y,type='l')
max(y)
sum(abs(y))

setwd('/WORK/lbyybl/WH/rvb/cor_uniq/MNase/danpos')
chr_size <- fread('mm10.chrom.size',col.names = c('chr','len'))
floor(chr_size$len[2]/10+1)*10
chr_size2 <- chr_size %>%
  mutate(len=floor(len/10+1)*10)
for (i in 1:nrow(chr_size)){
  if (chr_size$len[i]==chr_size2$len[i]){
    print(i)
  }else{
    print(paste('no equ',chr_size$len[i],chr_size2$len[i]))
  }
}
fwrite(chr_size2,'mm10.chrom.size2',sep = '\t',col.names = F)

setwd('/WORK/lbyybl/WH/rvb/ATAC/sample20191014/ATAC_QC')
fl <- fread('R21hATAC-1_FKDL190765001-1a_FL_count.txt',col.names = c('num','len'))
plot(fl$len,log(fl$num),type='h')
plot(fl$len,fl$num,type='h')
ggplot(fl,aes(len,num))+geom_col()+
  scale_x_continuous(limits = c(0,1000))
fl2 <- fread('../test_frag.txt',col.names = 'len')
ggplot(fl2,aes(len))+geom_freqpoly()
fl2_summ <- fl2 %>%
  group_by(len) %>%
  summarise(num=n())
undebug(bamQC)
bamfile<-'/WORK/lbyybl/WH/rvb/ATAC/sample20191014/test/R205hATAC-1.bam'
bam_index<-'/WORK/lbyybl/WH/rvb/ATAC/sample20191014/test/R205hATAC-1.bam.bai'
bamQC(bamfile,bam_index,mitochondria = 'chrM', outPath=NULL)

tibble(x = 1:10, y = c(-5:-1,1:5)) %>% ggplot(aes(x = x)) + 
  geom_col(aes(y = y)) + 
  geom_text(aes(y = 0, label = x, hjust = case_when(y > 0 ~ 1, y < 0 ~ 0, TRUE ~ .5))) + 
  coord_flip()

setwd('/WORK/lbyybl/WH/rvb/MNase/sample20190930/test/test_gene')
merge <- fread('join_trans.bed',header = F,stringsAsFactors = F,
              col.names = c('symbol','count','chr','st','en','score','strand'))
merge <- unique(merge)
merge <- merge %>%
  arrange(symbol)
#merge$length <- abs(merge$en-merge$st)
#merge <- unique(merge)
for (i in 1:(nrow(merge)-1)){
  if (merge$symbol[i]==merge$symbol[i+1]){
    if (merge$count[i]==merge$count[i+1]){
      if (merge$length[i]==merge$length[i+1]){
        merge[i+1,1:8]<- merge[i,1:8]
      }else if (merge$length[i]>merge$length[i+1]){
        merge[i,c(1:8)]<- merge[i+1,c(1:8)]
      }else if (merge$length[i]<merge$length[i+1]){
        merge[i+1,c(1:8)]<- merge[i,c(1:8)]
      }
    }else if (merge$count[i]>merge$count[i+1]){
      merge[i,c(1:8)]<- merge[i+1,c(1:8)]
    }else if (merge$count[i]<merge$count[i+1]){
      merge[i,c(1:8)]<- merge[i+1,c(1:8)]
    }
  }
}
nrow(unique(merge))
merge2 <- merge[,c(3:5,1,6:7)]
fwrite(merge2,'final_gene.bed',col.names = F,sep = '\t')
merge3 <- merge %>%
  arrange(-count)
merge3 <- merge3[1:4500,c(3:5,1,6:7)]
fwrite(merge3,'high_gene.bed',col.names = F,sep = '\t')
#--------------------------------------------------------
setwd('/WORK/lbyybl/WH/rvb/MNase/sample20190930/test/test_gene')
file <- fread('frg_len.txt',col.names = c('len'))
filecount <- data.frame('count'=seq(10,300,10),
                        'num'=rep(0,30))

for (i in seq(10,300,10)){
  assign(paste0('sts',i),file[file$len<i & file$len>(i-10),])
}

for (i in seq(10,300,10)){
  data <- get(paste0('sts',i))
  print(paste0(paste0('sts',i),'=>',nrow(data)))
  filecount[filecount$count==i,2] <- nrow(data)
  
}

ggplot(filecount,aes(x=count,y=num)) + geom_col()
mean(sum(filecount$num*filecount$count)/sum(filecount$num))
filecount %>%
  dplyr::filter(count<181 & count > 151) %>%
  summarise(sum(num))
#---BinBinLai
setwd('/DATA3/lbyybl/download/MNase/BinBinLai')
file2 <- fread('frg_len.txt',col.names = c('len'))
filecount2 <- data.frame('count'=seq(10,300,10),
                        'num'=rep(0,30))

for (i in seq(10,300,10)){
  assign(paste0('sts2',i),file2[file2$len<i & file2$len>(i-10),])
}

for (i in seq(10,300,10)){
  data <- get(paste0('sts2',i))
  print(paste0(paste0('sts2',i),'=>',nrow(data)))
  filecount2[filecount2$count==i,2] <- nrow(data)
  
}

ggplot(filecount2,aes(x=count,y=num)) + geom_col()
mean(sum(filecount2$num*filecount2$count)/sum(filecount2$num))
filecount2 %>%
  dplyr::filter(count<181 & count >130) %>%
  summarise(sum(num))

fwrite(filecount,'mydata.txt')
fwrite(filecount2,'BinBindata.txt')
filecount <- fread('mydata.txt',stringsAsFactors = F)

library(nucleR)
library(dplyr)
library(IRanges)
library(GenomicRanges)
library(ShortRead)
library(doParallel)
library(ggplot2)
library(magrittr)
library(NucDyn)
data(readsG2_chrII)
data(readsM_chrII)
data(nuc_chrII)
dyn <- nucleosomeDynamics(setA=readsG2_chrII, setB=readsM_chrII)
findHotspots(dyn, nuc_chrII)

#---
library(data.table)
library(dplyr)
freg <- fread('h05rep1-2.txt')
names(freg) <- 'len'
freg_freq <- freg %>%
  dplyr::group_by(len) %>%
  summarise(n=n())
freg_freq$freq <- freg_freq$n/nrow(freg)
plot(freg_freq$len,freg_freq$n,type='l',ylab='frequency',xlab='frengment length')
#axis(1,0:1000)
label<-data.frame(x=c(130,290,460,200,390,567),y=rep(0,6))
ggplot(freg_freq,aes(x=len,y=freq)) + geom_line() + scale_x_continuous(limits=c(0,1000))+
  geom_vline(xintercept =c(130,290,460),color='blue')+
  geom_vline(xintercept = c(200,390,567),color='red')+xlab('frengment length')+
  ylab('frequency')+ geom_text(data=label,aes(x,y,label=x),size=2)+
  theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(), panel.background=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2))
ggsave('iswi_freg.pdf')
abline(v=130,lwd=1,col="blue") # 134
abline(v=290,lwd=1,col="blue") # 308
abline(v=460,lwd=1,col="blue") # 447
abline(v=200,lwd=1,col="red") # 177
abline(v=390,lwd=1,col="red") # 367
abline(v=567,lwd=1,col="red") 
# y.predict <- predict(loess(freg_freq$freq ~ freg_freq$len), se = TRUE)
# plot(freg_freq$len,y.predict$fit)
jz <- freg_freq %>%
  filter(len >558 & len <615)
n<-which(jz$freq==max(jz$freq))
jz[n,]
# library(ggplot2)
# g <- ggplot(freg,aes(x=len,y=..count../nrow(freg))) + geom_freqpoly(binwidth = 1)+
#   scale_x_continuous(limits = c(0,1000)) + geom_vline(xintercept = 100) +
#   geom_vline(xintercept = 180) + geom_vline(xintercept = 247) +
#   geom_vline(xintercept = 315) + geom_vline(xintercept = 473) +
#   geom_vline(xintercept = 558) + geom_vline(xintercept = 615)
# g1 <- ggplotGrob(g)
# g + annotation_custom(grob = g1)

freg2 <- fread('R2DoxATAC_merge_frag.txt')
names(freg2) <- 'len'
freg_freq <- freg2 %>%
  dplyr::group_by(len) %>%
  summarise(n=n())
freg_freq$freq <- freg_freq$n/nrow(freg)
plot(freg_freq$len,freg_freq$freq,type='l',ylab='frequency',xlab='frengment length')
#axis(1,0:1000)
label<-data.frame(x=c(150,380,30,210,80),y=rep(0,5))
# jz <- freg_freq %>%
#   filter(len >558 & len <615)
# n<-which(jz$freq==max(jz$freq))
# jz[n,]
ggplot(freg_freq,aes(x=len,y=freq)) + geom_line() + scale_x_continuous(limits=c(0,1000))+
  geom_vline(xintercept =c(150,380,30),color='blue')+
  geom_vline(xintercept = c(210,80),color='red')+xlab('frengment length')+
  ylab('frequency')+ geom_text(data=label,aes(x,y,label=x),size=2)+
  theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(), panel.background=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2))
ggsave('dox_freg.pdf')
ggplot(freg2,aes(x=len,y=..count../nrow(freg))) + geom_freqpoly(binwidth = 1)+
  scale_x_continuous(limits = c(0,1000)) + geom_vline(xintercept = 100) +
  geom_vline(xintercept = 180) + geom_vline(xintercept = 247) +
  geom_vline(xintercept = 315) + geom_vline(xintercept = 473) +
  geom_vline(xintercept = 558) + geom_vline(xintercept = 615)

