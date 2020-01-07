#--see the coorlation of RNA-seq and chip-seq
library(data.table)
library(dplyr)
library(clusterProfiler)
library(RColorBrewer)
library(gplots)
library(reshape2)
library(tidyr)
library(ggplot2)
setwd('/DATA2/work/lbyybl/coorlaborate/YB/coor_4c_Rna')
read_file <- function(file){
  data <- fread(file,col.names = c('chr','st','en','id','score','strand','num'))
  data_total <- sum(data$num)
  data <- data %>%
    mutate(rpkm=num*1000000/(data_total*as.numeric(en-st)))
  data <- as.data.frame(data)
  return(data)
} 

Ac12_rna <- read_file('Ac12_rna.count')
Ad38nodox6_rna <- read_file('Ad38nodox_rna.count')
Ac12_4c <- read_file('AC12_4c.count')
Ad38nodox6_4c <- read_file('AD38noDOX6d_4c.count')


merge <- merge(Ac12_rna,Ad38nodox6_rna,by='id')
merge <- merge[,c(1,8,15)]
colnames(merge) <- c('id','Ac12_rna','Ad38nodox6_rna')
merge <- merge(Ac12_4c,merge,by='id')
merge <- merge[,c(1,8,9,10)]
colnames(merge)[2] <- 'Ac12_4c'
merge <- merge(Ad38nodox6_4c,merge,by='id')
merge <- merge[,c(1,8,9,10,11)]
colnames(merge)[2] <- 'Ad38nodox6_4c'

#merge <- merge[order(merge$rna),]
head(merge)
merge <- merge %>%
  filter(Ac12_rna+Ad38nodox6_rna+Ac12_4c+Ad38nodox6_4c>0)
# merge$h3k27 <- (merge$h3k27-mean(merge$h3k27))/var(merge$h3k27)
# merge$rna <- (merge$rna-mean(merge$rna))/var(merge$rna)
merge3 <- merge %>% arrange(Ac12_rna)
merge3$id <- 1:nrow(merge3)
merge3 <- gather(merge3,class,num,c('Ac12_rna','Ad38nodox6_rna','Ad38nodox6_4c','Ac12_4c'))
head(merge3)
merge3_Ac12_4c <- merge3 %>% filter(class=='Ac12_4c')
merge3_Ad38nodox6_rna <- merge3 %>% filter(class=='Ad38nodox6_rna')
merge3_Ac12_rna <- merge3 %>% filter(class=='Ac12_rna')
ggplot(merge3_Ac12_rna,aes(x=id,y=class,fill=num))+geom_tile()+
  scale_fill_distiller(palette = "RdYlBu",values=c(0,0.5,0.50001,0.5001,0.501,0.51,1),limits=c(-2.7,2.7))+
  theme(plot.title=element_text(size=40,hjust=0.5), 
        panel.border = element_rect(colour="black",  fill=NA, size=1),
        axis.title=element_text(size=30),
        axis.text=element_text(size=28),
        legend.title=element_text(size=25),
        panel.background=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.text=element_text(size=22))
ggsave('Ac12_rna_rank.pdf',width = 12,height = 1.5)
ggplot(merge3_Ac12_4c,aes(x=id,y=class,fill=num))+geom_tile()+
  scale_fill_distiller(palette = "RdYlBu",values=c(0,0.5,0.51,1),limits=c(-0.041,0.041))+
  theme(plot.title=element_text(size=40,hjust=0.5), 
        panel.border = element_rect(colour="black",  fill=NA, size=1),
        axis.title=element_text(size=30),
        axis.text=element_text(size=28),
        legend.title=element_text(size=25),
        panel.background=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.text=element_text(size=22))
ggsave('Ac12_4c.pdf',width = 12,height = 1.5)

merge3 <- merge %>% arrange(Ad38nodox6_rna)
merge3$id <- 1:nrow(merge3)
merge3 <- gather(merge3,class,num,c('Ac12_rna','Ad38nodox6_rna','Ad38nodox6_4c','Ac12_4c'))
head(merge3)
merge3_Ad38nodox6_4c <- merge3 %>% filter(class=='Ad38nodox6_4c')
merge3_Ad38nodox6_rna <- merge3 %>% filter(class=='Ad38nodox6_rna')
merge3_Ac12_rna <- merge3 %>% filter(class=='Ac12_rna')
ggplot(merge3_Ad38nodox6_rna,aes(x=id,y=class,fill=num))+geom_tile()+
  scale_fill_distiller(palette = "RdYlBu",values=c(0,0.5,0.50001,0.5001,0.501,0.51,1),limits=c(-2.7,2.7))+
  theme(plot.title=element_text(size=40,hjust=0.5), 
        panel.border = element_rect(colour="black",  fill=NA, size=1),
        axis.title=element_text(size=30),
        axis.text=element_text(size=28),
        legend.title=element_text(size=25),
        panel.background=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.text=element_text(size=22))
ggsave('Ad38nodox6_rna_rank.pdf',width = 12,height = 1.5)
ggplot(merge3_Ad38nodox6_4c,aes(x=id,y=class,fill=num))+geom_tile()+
  scale_fill_distiller(palette = "RdYlBu",values=c(0,0.5,0.505,1),limits=c(-0.041,0.041))+
  theme(plot.title=element_text(size=40,hjust=0.5), 
        panel.border = element_rect(colour="black",  fill=NA, size=1),
        axis.title=element_text(size=30),
        axis.text=element_text(size=28),
        legend.title=element_text(size=25),
        panel.background=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.text=element_text(size=22))
ggsave('Ad38nodox6_4c.pdf',width = 12,height = 1.5)