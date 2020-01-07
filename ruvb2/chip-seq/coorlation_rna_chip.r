#--see the coorlation of RNA-seq and chip-seq
library(data.table)
library(dplyr)
library(clusterProfiler)
library(RColorBrewer)
library(gplots)
library(reshape2)
library(tidyr)
library(ggplot2)
setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/bw_bam/Ruvbl2/merge/bam/express_coor')
read_file <- function(file){
  data <- fread(file,col.names = c('chr','st','en','id','score','strand','num'))
  data_total <- sum(data$num)
  data <- data %>%
    mutate(rpkm=num*1000000/(data_total*as.numeric(en-st)))
  data <- as.data.frame(data)
  return(data)
} 
rna_seq_file <- '/DATA/work/lbyybl/wh/ruvb2/chip-seq/bw_bam/Pol2/merge/bam/express_coor/RNA_seq.count'
rna_seq <- read_file(rna_seq_file)
h3k27 <- read_file('H3k27ac.count')
ruvbl2 <- read_file('Ruvb2.count')
# rna_seq <- fread(rna_seq_file,col.names = c('chr','st','en','id','score','strand','num'))
# h3k27 <- fread('H3k27ac.count',col.names = c('chr','st','en','id','score','strand','num'))
# rna_seq_total <- sum(rna_seq$num)
# h3k27_total <- sum(h3k27$num)
# rna_seq <- rna_seq %>%
#   mutate(rpkm=num*1000000/(rna_seq_total*as.numeric(en-st)))
# h3k27 <- h3k27 %>%
#   mutate(rpkm=num*1000000/(h3k27_total*as.numeric(en-st)))
merge <- merge(h3k27,rna_seq,by='id')
merge <- merge[,c(1,8,15)]
colnames(merge) <- c('id','h3k27','rna')
merge <- merge(ruvbl2,merge,by='id')
merge <- merge[,c(1,8,9,10)]
colnames(merge)[2] <- 'ruvbl2'
merge <- merge[order(merge$rna),]
head(merge)
# merge$h3k27 <- (merge$h3k27-mean(merge$h3k27))/var(merge$h3k27)
# merge$rna <- (merge$rna-mean(merge$rna))/var(merge$rna)
merge3 <- merge %>% arrange(rna)
merge3$id <- 1:nrow(merge3)
merge3 <- gather(merge3,class,num,c('ruvbl2','h3k27','rna'))
head(merge3)
merge3_rna <- merge3 %>% filter(class=='rna')
merge3_h3k27 <- merge3 %>% filter(class=='h3k27')
merge3_ruvbl2 <- merge3 %>% filter(class=='ruvbl2')
ggplot(merge3_rna,aes(x=id,y=class,fill=num))+geom_tile()+
  scale_fill_distiller(palette = "RdYlBu",values=c(0,0.5,0.50001,1),limits=c(-1.6,1.6))+
  theme(plot.title=element_text(size=40,hjust=0.5), 
        panel.border = element_rect(colour="black",  fill=NA, size=1),
        axis.title=element_text(size=30),
        axis.text=element_text(size=28),
        legend.title=element_text(size=25),
        panel.background=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.text=element_text(size=22))
ggsave('rna_rank.pdf',width = 12,height = 1.5)
ggplot(merge3_h3k27,aes(x=id,y=class,fill=num))+geom_tile()+
  scale_fill_distiller(palette = "RdYlBu",values=c(0,0.5,0.52,1),limits=c(-0.042,0.042))+
  theme(plot.title=element_text(size=40,hjust=0.5), 
        panel.border = element_rect(colour="black",  fill=NA, size=1),
        axis.title=element_text(size=30),
        axis.text=element_text(size=28),
        legend.title=element_text(size=25),
        panel.background=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.text=element_text(size=22))
ggsave('h3k27.pdf',width = 12,height = 1.5)
ggplot(merge3_ruvbl2,aes(x=id,y=class,fill=num))+geom_tile()+
  scale_fill_distiller(palette = "RdYlBu",values=c(0,0.5,0.51,1),limits=c(-0.085,0.085))+
  theme(plot.title=element_text(size=40,hjust=0.5), 
        panel.border = element_rect(colour="black",  fill=NA, size=1),
        axis.title=element_text(size=30),
        axis.text=element_text(size=28),
        legend.title=element_text(size=25),
        panel.background=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.text=element_text(size=22))
ggsave('ruvbl2.pdf',width = 12,height = 1.5)
