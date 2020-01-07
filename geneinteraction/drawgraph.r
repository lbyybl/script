setwd("/DATA/work/lbyybl/ypjiang/geneloop/pol2_loop/old_20180817/result/graph/standard_define")
library(ggplot2)
dir()
drawheatmap <- function(data,min_va,max_va){
  ggplot(data=data, aes(x=ID1, y=ID2, fill = log(interaction))) + geom_raster() +
    #scale_fill_gradientn(colors = c("navy","royalblue4", "white", "orangered", "orangered2")) +
    scale_fill_gradientn(colors = c("royalblue4", "white", "orangered", "orangered2"),
                         values = c(0,0.1,0.5,1),limits=c(min_va,max_va)
    ) + 
    #scale_fill_brewer(palette="Pastel2", limits=c("trt1", "trt2", "ctrl"))
    # scale_fill_manual(breaks=c("0","0.5","1.5","7.5"),
    #                   values=c("dodgerblue", "white", "orange", "orangered")) +
    scale_y_discrete(expand=c(0,0),limits=c("promoter","insulator","gene_body","typical_enhancer","super_enhancer")) +
    scale_x_discrete(expand=c(0,0),position = "top",limits=c("super_enhancer","typical_enhancer","gene_body","insulator","promoter")) +
    theme(plot.title=element_text(hjust=0.5),
          panel.grid=element_blank(), panel.background=element_blank()) +
    theme(legend.position="bottom",
          title = element_blank() ,
          legend.key.width = unit(3,"line"),
          legend.key.height = unit(0.1,"cm"),
          legend.text = element_text(
            size = 7,
            hjust = 1.2,
            face = "italic",
            colour = "black",
            angle = 30
          ),
          axis.line=element_blank(),
          plot.title=element_text(colour="black", size=10),
          axis.text.x=element_text(angle=60, hjust=0,vjust = 1.5, #face = "bold",
                                   colour="black", size=12),
          axis.text.y=element_text(#face = "bold",
            colour="black", size=12),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank()) + 
    labs(fill=NULL ) 
}

#data <- read.table("interaction_wt_format_final.stastic",header = F)
#View(data)
#ggplot(data = data, aes(x=data$V1,y=data$V2,fill=data$V3)) + geom_raster()
#dense <- read.table("dense2.bed",header = F)
#ggplot(data = dense, aes(x=dense$V1,y=dense$V2,fill=dense$V3)) + geom_raster()

#--- draw meta heatmap
setwd("/DATA/work/lbyybl/ypjiang/geneloop/pol2_loop/old_20180817/result/graph/profji_adviced_define")

kointa <- read.table("interaction_wt_format_final2.stastic")
wtinta <- read.table("interaction_ko_format_final2.stastic")
colnames(kointa) <- c("ID1","ID2","interaction")
colnames(wtinta) <- c("ID1","ID2","interaction")
kointa$interaction <- kointa$interaction/1.04826
wtinta$interaction <- wtinta$interaction
va_max<-max(log(kointa$interaction),log(wtinta$interaction))
va_min<-min(log(kointa$interaction),log(wtinta$interaction))
drawheatmap(kointa,va_min,va_max)
ggsave("profji_ko.jpeg", width = 5,height = 5.3,dpi = 1000)
#jpeg("wt.jpeg")
drawheatmap(wtinta,va_min,va_max)
ggsave("profji_wt.jpeg", width = 5,height = 5.3,dpi = 1000)
#--- a function to get the divition of two matrix, data1 should bigget than data2
inte_diff <- function(data1,data2){
  data_diff <- data1
  data_diff$interaction <- (data1$interaction)/(data2$interaction)
  return(data_diff)
}
dif_kw <- inte_diff(kointa,wtinta)
# inte_diff <- kointa
# inte_diff$interaction <- (kointa$interaction/18)/(wtinta$interaction/17)
#jpeg("difference.jpeg")
drawdiffheatmap(dif_kw)
ggsave("profji_difference.jpeg", width = 5,height = 5.3,dpi = 1000)
#dev.off()
#--- a matrix remove pol2 peak
#ko_rm_p2 <- kointa[-which((kointa$ID1=="pol2_peak" | kointa$ID2=="pol2_peak")),]
#wt_rm_p2 <- wtinta[-which((wtinta$ID1=="pol2_peak" | wtinta$ID2=="pol2_peak")),]
#rm_p2_diff <- inte_diff(ko_rm_p2, wt_rm_p2)
#jpeg("rm_pol2_peak_ko.jpeg")
#drawheatmap(ko_rm_p2,"ko")
#dev.off()
#jpeg("rm_pol2_peak_wt.jpeg")
#drawheatmap(wt_rm_p2,"wt")
#dev.off()
#jpeg("rm_pol2_peak_diff.jpeg")
#drawdiffheatmap(rm_p2_diff,"ko/wt difference")
#dev.off()


drawdiffheatmap <- function(data){
  ggplot(data=data, aes(x=ID1, y=ID2, fill = interaction)) + geom_raster() +
    scale_fill_gradientn(colors = c("royalblue4", "white", "orangered"),
                         values = c(0,0.5,1),limits=c(0.85,1.15)
    ) + 
    scale_y_discrete(expand=c(0,0),limits=c("promoter","insulator","gene_body","typical_enhancer","super_enhancer")) +
    scale_x_discrete(expand=c(0,0),position = "top",limits=c("super_enhancer","typical_enhancer","gene_body","insulator","promoter")) +
    theme(plot.title=element_text(hjust=0.5),
          panel.grid=element_blank(), panel.background=element_blank()) +
    theme(legend.position="bottom",
          title = element_blank() ,
          legend.key.width = unit(3,"line"),
          legend.key.height = unit(0.1,"cm"),
          legend.text = element_text(
            size = 7,
            hjust = 1.2,
            face = "italic",
            colour = "black",
            angle = 30
          ),
          axis.line=element_blank(),
          plot.title=element_text(colour="black", size=10),
          axis.text.x=element_text(angle=60, hjust=0,vjust = 1.5, #face = "bold",
                                   colour="black", size=12),
          axis.text.y=element_text(#face = "bold",
            colour="black", size=12),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank()) + 
    labs(fill=NULL ) 
  }

#--- draw normalization heatmap

drawnormheatmap <- function(data,min_va,max_va){
  ggplot(data=data, aes(x=ID1, y=ID2, fill = interaction)) + geom_raster() +
    scale_fill_gradientn(colors = c("royalblue4", "white", "orangered", "orangered2"),
                         values = c(0,0.1,0.5,1)#,limits=c(min_va,max_va)
    ) + 
    scale_y_discrete(expand=c(0,0),limits=c("Promoter","Insulator","Gene_body","Typical_enhancer","Super_enhancer")) +
    scale_x_discrete(expand=c(0,0),position = "top",limits=c("Super_enhancer","Typical_enhancer","Gene_body","Insulator","Promoter")) +
    theme(plot.title=element_text(hjust=0.5),
          panel.grid=element_blank(), panel.background=element_blank()) +
    theme(legend.position="bottom",
          title = element_blank() ,
          legend.key.width = unit(3,"line"),
          legend.key.height = unit(0.1,"cm"),
          legend.text = element_text(
            size = 7,
            hjust = 1.2,
            face = "italic",
            colour = "black",
            angle = 30
          ),
          axis.line=element_blank(),
          plot.title=element_text(colour="black", size=10),
          axis.text.x=element_text(angle=60, hjust=0,vjust = 1.5, #face = "bold",
                                   colour="black", size=12),
          axis.text.y=element_text(#face = "bold",
            colour="black", size=12),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank()) + 
    labs(fill=NULL ) 
}


normaziinta_ko <- read.table("interactionnomalize",stringsAsFactors = F)
colnames(normaziinta_ko) <- c("ID1","ID2","interaction","ID1_LEN","ID2_LEN")
normaziinta_ko<- dplyr::mutate(normaziinta_ko, interaction=interaction/(ID1_LEN*ID2_LEN))
va_mi<-min(normaziinta_ko$interaction)
va_mA<-max(normaziinta_ko$interaction)
drawnormheatmap(normaziinta_ko,va_mi,va_ma)

