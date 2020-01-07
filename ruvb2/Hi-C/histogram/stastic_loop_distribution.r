# this is to stastic the distribution of loop pet for hichip and ocean-c
#---ocean-c
library(data.table)
library(ggplot2)
setwd("/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/loop/histogram")
wt <- 'RuvbL2_wt_loop.bed'
ko <- 'RuvbL2_ko_loop.bed'
wt_oceanc <- fread(wt)
wt_ocp <- data.frame("num"=wt_oceanc$V2,
                     "class"="wt")
ko_oceanc <- fread(ko)
ko_ocp <- data.frame("num"=ko_oceanc$V2,
                     "class"="ko")
ocp <- rbind(wt_ocp,ko_ocp)
ocp$num <- round(ocp$num)
ggplot(ocp,aes(x=num,y=log2(..count..),fill=class))+ geom_bar(position = "dodge2") +
  scale_y_continuous(limits = c(0,8)) +
  scale_x_continuous(limits = c(0,50))+
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
setwd('/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/loop/histogram')
ggsave('loop_distribution.pdf',width = 5,height = 3)