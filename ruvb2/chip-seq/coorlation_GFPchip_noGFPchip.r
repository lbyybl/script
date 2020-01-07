#--see the coorlation of RNA-seq and chip-seq
library(data.table)
library(dplyr)
library(clusterProfiler)
library(RColorBrewer)
library(gplots)
library(reshape2)
library(tidyr)
library(ggplot2)
setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/unGFP-chip/graph/coor/sample')
read_file <- function(file){
  data <- fread(file,col.names = c('chr','st','en','num'))
  data_total <- sum(data$num)
  data <- data %>%
    mutate(rpkm=num*1000000/(data_total*as.numeric(en-st)))
  data <- as.data.frame(data)
  data$id <- paste0(data$chr,"_",data$st,"_",data$en)
  return(data)
} 
GFP_Ruvb1 <- read_file('GFP_Ruvb1.count')
GFP_Ruvbl2 <- read_file('GFP_Ruvbl2.count')
ruvbl1 <- read_file('Ruvb1.count')
ruvbl2 <- read_file('Ruvb2.count')
h3k27ac <- read_file('H3k27ac.count')
H3K9me3 <- read_file('H3K9me3.count')

mergeG12 <- left_join(GFP_Ruvb1,GFP_Ruvbl2,by='id') 

mergeG12 <- mergeG12[,c(5,11,6)]
colnames(mergeG12) <- c('Gr1','Gr2','id')
merge1G12 <- left_join(ruvbl1,mergeG12,by='id')
merge1G12 <- merge1G12[,c(5,7,8,6)]
colnames(merge1G12)[1] <- 'ruvbl1'
merge12G12 <- left_join(ruvbl2,merge1G12,by='id')
merge12G12 <- merge12G12[,c(6,5,7,8,9)]
colnames(merge12G12)[2] <- 'ruvbl2'
merge12G12 <- left_join(h3k27ac,merge12G12,by='id')
merge12G12 <- merge12G12[,c(6,5,7,8,9,10)]
colnames(merge12G12)[2] <- 'h3k27ac'
merge12G12 <- left_join(H3K9me3,merge12G12,by='id')
merge12G12 <- merge12G12[,c(6,5,7,8,9,10,11)]
colnames(merge12G12)[2] <- 'H3K9me3'
merge12G12 <- merge12G12[order(merge12G12$Gr2),]

head(merge12G12)
# heatmap_matrix <- merge12G12[2:5]
# rownames(heatmap_matrix) <- merge12G12$id
# mycol <- colorpanel(1000,"blue","white","red")
# pheatmap(as.matrix(heatmap_matrix),color = mycol,cluster_row = FALSE,
#          cluster_col = FALSE,#scale = 'row',
#          filename = 'AD38noDOXAD38withDOX_heatmap.pdf')


merge3 <- merge12G12 %>% arrange(Gr2)
merge3 <- merge3[round(nrow(merge3)*0.03):round(nrow(merge3)-round(nrow(merge3)*0.03)),]
merge3$id <- 1:nrow(merge3)
merge3 <- gather(merge3,class,num,c('Gr1','Gr2','ruvbl1','ruvbl2','h3k27ac','H3K9me3'))
head(merge3)
merge3_Gr1 <- merge3 %>% filter(class=='Gr1')
merge3_Gr2 <- merge3 %>% filter(class=='Gr2')
merge3_ruvbl2 <- merge3 %>% filter(class=='ruvbl2')
merge3_ruvbl1 <- merge3 %>% filter(class=='ruvbl1')
merge3_h3k27 <- merge3 %>% filter(class=='h3k27ac')
merge3_h3k9 <- merge3 %>% filter(class=='H3K9me3')
ggplot(merge3_Gr1,aes(x=id,y=class,fill=num))+geom_tile()+
  scale_fill_distiller(palette = "RdYlBu",values=c(0,0.02,0.05,1),limits=c(0,0.6))+
  theme(plot.title=element_text(size=40,hjust=0.5), 
        panel.border = element_rect(colour="black",  fill=NA, size=1),
        axis.title=element_text(size=30),
        axis.text=element_text(size=28),
        legend.title=element_text(size=25),
        panel.background=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.text=element_text(size=22))
ggsave('GFP_Ruvbl1.pdf',width = 12,height = 1.5)
ggplot(merge3_Gr2,aes(x=id,y=class,fill=num))+geom_tile()+
  scale_fill_distiller(palette = "RdYlBu",values=c(0,0.018,0.05,1),limits=c(0,0.6))+
  theme(plot.title=element_text(size=40,hjust=0.5), 
        panel.border = element_rect(colour="black",  fill=NA, size=1),
        axis.title=element_text(size=30),
        axis.text=element_text(size=28),
        legend.title=element_text(size=25),
        panel.background=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.text=element_text(size=22))
ggsave('GFP_Ruvbl2.pdf',width = 12,height = 1.5)
ggplot(merge3_ruvbl1,aes(x=id,y=class,fill=num))+geom_tile()+
  scale_fill_distiller(palette = "RdYlBu",values=c(0,0.02,0.05,1),limits=c(0,0.6))+
  theme(plot.title=element_text(size=40,hjust=0.5), 
        panel.border = element_rect(colour="black",  fill=NA, size=1),
        axis.title=element_text(size=30),
        axis.text=element_text(size=28),
        legend.title=element_text(size=25),
        panel.background=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.text=element_text(size=22))
ggsave('Ruvbl1.pdf',width = 12,height = 1.5)
ggplot(merge3_ruvbl2,aes(x=id,y=class,fill=num))+geom_tile()+
  scale_fill_distiller(palette = "RdYlBu",values=c(0,0.018,0.05,1),limits=c(0,0.6))+
  theme(plot.title=element_text(size=40,hjust=0.5), 
        panel.border = element_rect(colour="black",  fill=NA, size=1),
        axis.title=element_text(size=30),
        axis.text=element_text(size=28),
        legend.title=element_text(size=25),
        panel.background=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.text=element_text(size=22))
ggsave('Ruvbl2.pdf',width = 12,height = 1.5)
ggplot(merge3_h3k27,aes(x=id,y=class,fill=num))+geom_tile()+
  scale_fill_distiller(palette = "RdYlBu",values=c(0,0.05,0.52,1),limits=c(0,0.05))+
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
ggplot(merge3_h3k9,aes(x=id,y=class,fill=num))+geom_tile()+
  scale_fill_distiller(palette = "RdYlBu",values=c(0,0.25,0.52,1),limits=c(0,0.05))+
  theme(plot.title=element_text(size=40,hjust=0.5), 
        panel.border = element_rect(colour="black",  fill=NA, size=1),
        axis.title=element_text(size=30),
        axis.text=element_text(size=28),
        legend.title=element_text(size=25),
        panel.background=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.text=element_text(size=22))
ggsave('h3k9.pdf',width = 12,height = 1.5)