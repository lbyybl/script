setwd('/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/loop/MAplot')
library(data.table)
library(dplyr)
library(ggplot2)
#library()
ko_loop <- '/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/loop/MAplot/ruvbl2_ko_loop.bed'
wt_loop <- '/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/loop/MAplot/ruvbl2_wt_loop.bed'

readfile <- function(file){
  data <- fread(file,col.names = c('loc','inte'))
}

ko_file <- readfile(ko_loop)
colnames(ko_file) <- c('loc','ko_int')
wt_file <- readfile(wt_loop)
colnames(wt_file) <- c('loc','wt_int')
merge_file <- merge(ko_file,wt_file,by='loc')
merge_file$logFC <- log2(merge_file$ko_int/merge_file$wt_int)
merge_file$logCPM <- log2(merge_file$ko_int*merge_file$wt_int)

ggplot(merge_file,aes(x=logCPM,y=logFC))+ geom_point(size=0.3,alpha=0.7,color='red') +
  geom_hline(aes(yintercept = 0),colour="black",size=0.2,linetype="dashed")+
  #geom_hline(aes(yintercept = log2(2)),color="cyan3",size=0.2,linetype="dashed")+
  #geom_hline(aes(yintercept = log2(0.5)),color="cyan3",size=0.2,linetype="dashed")+
  labs(x="log2 (Degron*untreated)/2",y="log2 Degron/untreated FC",color="") +
  scale_y_continuous(expand=c(0,0),limits = c(-3,3)) +
  #scale_color_manual(values = c("grey67","blue","blue","red","red"))+
  scale_color_manual(values = c("red","red","red","red","red"))+
  #scale_color_manual(values = c("grey67","grey67","blue","red"))+
  scale_x_continuous(limits = c(0,10))+
  theme(
    legend.key.width = unit(0.1,"line"),
    legend.key.height = unit(1,"cm"),
    legend.text = element_text(
      size = 10,
      hjust = 0,
      #face = "italic",
      colour = "black",
      angle = 0
    ),
    
    plot.title=element_text(colour="black", size=10),
    axis.line.x = element_blank(),#element_line(color="black", size = 0.5),
    axis.line.y = element_blank(), #element_line(color="black", size = 0.5),
    axis.text.x=element_text(angle=0, hjust=0.5,vjust = 0,
                             colour="black", size=12),
    axis.text.y=element_text(#face = "bold",
      colour="black", size=12),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    panel.background=element_blank()
  )

ggsave('Ruvbl2_loop_MA.pdf',width = 3, height = 3)


ko_file$sample <- 'ko'
colnames(ko_file)[2] <- 'int'
wt_file$sample <- 'wt'
colnames(wt_file)[2] <- 'int'
merge_data <- bind(ko_file,wt_file)

ggplot(merge_data,aes(x=sample,y=int,fill=sample))+geom_boxplot(outlier.color = NA)+
  scale_y_continuous(limits = c(0,8)) +
  scale_x_discrete(limits=c('wt','ko'))+
  theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(), panel.background=element_blank()) +
  theme(
    title = element_blank() ,
    legend.key.width = unit(2,"line"),
    legend.key.height = unit(0.5,"cm"),
    legend.text = element_text(
      size = 7,
      hjust = 1.2,
      face = "italic",
      colour = "black"
    ),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    plot.title=element_text(colour="black", size=10),
    axis.text.x=element_text(#face = "bold",
      colour="black", size=12),
    axis.text.y=element_text(#face = "bold",
      colour="black", size=12),
    panel.background=element_blank(),
    plot.background=element_blank()) +
  labs(fill=NULL )
ggsave('stas_MA.pdf',width = 3.2, height = 4)
t.test(ko_file$int,wt_file$int)$p.value