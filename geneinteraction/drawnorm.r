setwd("/DATA/work/lbyybl/ypjiang/geneloop/pol2_loop/normalization_2018_9_14/result/50K/graph")
ko_reads <- 154213532 #pol2 
wt_reads <- 147119703 #pol2 
#ko_reads <- 146275644 #pol3 
#wt_reads <- 137392176 #pol3 
#ko_reads <- 132178869 # pol1
#wt_reads <- 134624874 # pol1
#ko_reads <- 172341303 # CTCF
#wt_reads <- 174128642 # CTCF
#ko_reads <- 97049732 # YY1
#wt_reads <- 95634595 #YY1
library(ggplot2)
dir()
drawdiffheatmap <- function(data,vi,va){
  ggplot(data=data, aes(x=ID1, y=ID2, fill = interaction)) + geom_raster() +
    scale_fill_gradientn(colors = c("royalblue4", "white", "orangered"),
                         values = c(0,0.5,1),limits=c(vi,va)
    ) + 
    scale_y_discrete(expand=c(0,0),limits=c("IN","TTS","GB","TSS","TE","SE")) +
    scale_x_discrete(expand=c(0,0),position = "top",limits=c("SE","TE","TSS","GB","TTS","IN")) +
    theme(plot.title=element_text(hjust=0.5),
          panel.grid=element_blank(), panel.background=element_blank()) +
    theme(legend.position="bottom",
          title = element_blank() ,
          legend.key.width = unit(2,"line"),
          legend.key.height = unit(0.5,"cm"),
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

#--- draw chazhi heatmap

drawchazhiheatmap <- function(data,min_va,max_va){
  ggplot(data=data, aes(x=ID1, y=ID2, fill = interaction)) + geom_raster() +
    scale_fill_gradientn(colors = c("DarkCyan", "PaleTurquoise","LightCyan","white", "Linen", "PeachPuff", "Red"),
                         values = c(0,0.49,0.499,0.5,0.501,0.51,1),limits=c(min_va,max_va)
    ) + 
    scale_y_discrete(expand=c(0,0),limits=c("IN","TTS","GB","TSS","TE","SE")) +
    scale_x_discrete(expand=c(0,0),position = "top",limits=c("SE","TE","TSS","GB","TTS","IN")) +
    theme(plot.title=element_text(hjust=0.5),
          panel.grid=element_blank(), panel.background=element_blank()) +
    theme(legend.position="bottom",
          title = element_blank() ,
          legend.key.width = unit(2,"line"),
          legend.key.height = unit(0.5,"cm"),
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
          #panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank()) + 
    labs(fill=NULL ) 
}

#--- draw normalization heatmap

drawnormheatmap <- function(data,min_va,max_va){
  ggplot(data=data, aes(x=ID1, y=ID2, fill = log(interaction))) + geom_raster() +
    scale_fill_gradientn(colors = c("royalblue4", "white", "orangered", "orangered2"),
                         values = c(0,0.1,0.5,1),limits=c(min_va,max_va)
    ) + 
    scale_y_discrete(expand=c(0,0),limits=c("IN","TTS","GB","TSS","TE","SE")) +
    scale_x_discrete(expand=c(0,0),position = "top",limits=c("SE","TE","TSS","GB","TTS","IN")) +
    theme(plot.title=element_text(hjust=0.5),
          panel.grid=element_blank(), panel.background=element_blank()) +
    theme(legend.position="bottom",
          title = element_blank() ,
          legend.key.width = unit(2,"line"),
          legend.key.height = unit(0.5,"cm"),
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
#--- a function to get the divition of two matrix, data1 should bigget than data2
inte_diff <- function(data1,data2){
  data_diff <- data1
  data_diff$interaction <- (data1$interaction)/(data2$interaction)
  return(data_diff)
}
chazhi_diff <- function(data1,data2){
  data_diff <- data1
  data_diff$interaction <- (data1$interaction)-(data2$interaction)
  return(data_diff)
}
Reso="50K"
filetype=".jpeg"
filetype2=".pdf"
normaziinta_ko <- read.table(paste0("interactionnomalize_ko",Reso),stringsAsFactors = F)
colnames(normaziinta_ko) <- c("ID1","ID2","interaction","ID2_LEN","ID1_LEN")
ko_generegion <- dplyr::filter(normaziinta_ko,ID1!="NG", ID1 != "IN",ID2!="NG", ID2!= "IN") #dplyr::filter(normaziinta_ko,ID1!="NG",ID2!="NG")
ko_gene_element_len <- unique(dplyr::select(ko_generegion,ID1,ID1_LEN))
ko_GR_int_len <- data.frame("ID1"=c("GR"),
                        "ID2"=c("GR"),
                        "interaction"=c(sum(ko_generegion[ko_generegion$ID1==ko_generegion$ID2,]$interaction,
                                            ko_generegion[ko_generegion$ID1!=ko_generegion$ID2,]$interaction/2)),
                        "ID1_LEN"=c(sum(ko_gene_element_len$ID1_LEN)),
                        "ID2_LEN"=c(sum(ko_gene_element_len$ID1_LEN)),stringsAsFactors = F)
ko_nongene <- dplyr::filter(normaziinta_ko,ID1=="NG",ID2=="NG")
ko_gene_nongene <- dplyr::filter(normaziinta_ko,(ID1!="NG" & ID2=="NG") | (ID1=="NG" & ID2!="NG"))
ko_GNG_int_len <-data.frame("ID1"=c("GR"),
                         "ID2"=c("NG"),
                         "interaction"=c(sum(ko_gene_nongene[ko_gene_nongene$ID1==ko_gene_nongene$ID2,]$interaction,
                                             ko_gene_nongene[ko_gene_nongene$ID1!=ko_gene_nongene$ID2,]$interaction/2)),
                        "ID1_LEN"=c(sum(ko_gene_element_len$ID1_LEN)),
                        "ID2_LEN"=c(ko_nongene$ID1_LEN),stringsAsFactors = F)
ko_NGG_int_len <- data.frame("ID1"=c("NG"),
                          "ID2"=c("GR"),
                          "interaction"=c(sum(ko_gene_nongene[ko_gene_nongene$ID1==ko_gene_nongene$ID2,]$interaction,
                                              ko_gene_nongene[ko_gene_nongene$ID1!=ko_gene_nongene$ID2,]$interaction/2)),
                          "ID1_LEN"=c(ko_nongene$ID1_LEN),
                          "ID2_LEN"=c(sum(ko_gene_element_len$ID1_LEN)),stringsAsFactors = F)
ko_gene_nongene_matrix <- rbind(ko_GR_int_len,ko_nongene,ko_GNG_int_len,ko_NGG_int_len)
normaziinta_wt <- read.table(paste0("interactionnomalize_wt",Reso),stringsAsFactors = F)
colnames(normaziinta_wt) <- c("ID1","ID2","interaction","ID2_LEN","ID1_LEN")
wt_generegion <- dplyr::filter(normaziinta_wt,ID1!="NG",ID2!="NG")
wt_gene_element_len <- unique(dplyr::select(wt_generegion,ID1,ID1_LEN))
wt_GR_int_len <-data.frame("ID1"=c("GR"),
                           "ID2"=c("GR"),
                           "interaction"=c(sum(wt_generegion[wt_generegion$ID1==wt_generegion$ID2,]$interaction,
                                               wt_generegion[wt_generegion$ID1!=wt_generegion$ID2,]$interaction/2)),
                           "ID1_LEN"=c(sum(wt_gene_element_len$ID1_LEN)),
                           "ID2_LEN"=c(sum(wt_gene_element_len$ID1_LEN)),stringsAsFactors = F)
wt_nongene <- dplyr::filter(normaziinta_wt,ID1=="NG",ID2=="NG")
wt_gene_nongene <- dplyr::filter(normaziinta_wt,(ID1!="NG" & ID2=="NG") | (ID1=="NG" & ID2!="NG"))
wt_GNG_int_len <-data.frame("ID1"=c("GR"),
                            "ID2"=c("NG"),
                            "interaction"=c(sum(wt_gene_nongene[wt_gene_nongene$ID1==wt_gene_nongene$ID2,]$interaction,
                                                wt_gene_nongene[wt_gene_nongene$ID1!=wt_gene_nongene$ID2,]$interaction/2)),
                            "ID1_LEN"=c(sum(wt_gene_element_len$ID1_LEN)),
                            "ID2_LEN"=c(wt_nongene$ID1_LEN),stringsAsFactors = F)
wt_NGG_int_len <- data.frame("ID1"=c("NG"),
                             "ID2"=c("GR"),
                             "interaction"=c(sum(wt_gene_nongene[wt_gene_nongene$ID1==wt_gene_nongene$ID2,]$interaction,
                                                 wt_gene_nongene[wt_gene_nongene$ID1!=wt_gene_nongene$ID2,]$interaction/2)),
                             "ID1_LEN"=c(wt_nongene$ID1_LEN),
                             "ID2_LEN"=c(sum(wt_gene_element_len$ID1_LEN)),stringsAsFactors = F)
wt_gene_nongene_matrix <- rbind(wt_GR_int_len,wt_nongene,wt_GNG_int_len,wt_NGG_int_len)

normaziinta_ko<- dplyr::mutate(normaziinta_ko, interaction=interaction*10^12/(as.numeric(ID1_LEN)*as.numeric(ID2_LEN)*ko_reads))
normaziinta_wt<- dplyr::mutate(normaziinta_wt, interaction=interaction*10^12/(as.numeric(ID1_LEN)*as.numeric(ID2_LEN)*wt_reads))
va_mi<-min(log(normaziinta_ko$interaction),log(normaziinta_wt$interaction))
va_ma<-max(log(normaziinta_ko$interaction),log(normaziinta_wt$interaction))
va_mi <- -17
va_ma <- -9

drawnormheatmap(normaziinta_ko,va_mi,va_ma)
ggsave(paste0("norm_ko",Reso,filetype), width = 5,height = 5.5,dpi = 1000)
ggsave(paste0("norm_ko",Reso,filetype2), width = 5,height = 5.5,dpi = 1000)
drawnormheatmap(normaziinta_wt,va_mi,va_ma)
ggsave(paste0("norm_wt",Reso,filetype), width = 5,height = 5.5,dpi = 1000)
ggsave(paste0("norm_wt",Reso,filetype2), width = 5,height = 5.5,dpi = 1000)
ko_wt_diff <- inte_diff(normaziinta_ko,normaziinta_wt)
diff_va <- max(ko_wt_diff$interaction)
diff_mi <- min(ko_wt_diff$interaction)
diff_va <- max(diff_va,1+abs(1-diff_mi))
diff_mi <- min(diff_mi,2-diff_va)
diff_mi <- 0.6 #0.65
diff_va <- 1.4 #1.35
drawdiffheatmap(ko_wt_diff,diff_mi,diff_va)
ggsave(paste0("norm_diff",Reso,filetype), width = 5,height = 5.5,dpi = 1000)
ggsave(paste0("norm_diff",Reso,filetype2), width = 5,height = 5.5,dpi = 1000)
ko_wt_chazhi <- chazhi_diff(normaziinta_ko,normaziinta_wt)
chazhi_va <- max(ko_wt_chazhi$interaction)
chazhi_mi <- min(ko_wt_chazhi$interaction)
chazhi_va <- max(chazhi_va,abs(chazhi_mi))
chazhi_mi <- -3.9e-05 #-1.7e-5
chazhi_va <- 3.9e-05 #1.7e-5
drawchazhiheatmap(ko_wt_chazhi,chazhi_mi,chazhi_va)
ggsave(paste0("norm_chazhi",Reso,filetype), width = 5,height = 5.5,dpi = 1000)
ggsave(paste0("norm_chazhi",Reso,filetype2), width = 5,height = 5.5,dpi = 1000)

drawdiffGNGheatmap <- function(data,vi,va){
  ggplot(data=data, aes(x=ID1, y=ID2, fill = interaction)) + geom_raster() +
    scale_fill_gradientn(colors = c("royalblue4", "white", "orangered"),
                         values = c(0,0.5,1),limits=c(vi,va)
    ) + 
    scale_y_discrete(expand=c(0,0),limits=c("GR","NG")) +
    scale_x_discrete(expand=c(0,0),position = "top",limits=c("NG","GR")) +
    theme(plot.title=element_text(hjust=0.5),
          panel.grid=element_blank(), panel.background=element_blank()) +
    theme(legend.position="bottom",
          title = element_blank() ,
          legend.key.width = unit(2,"line"),
          legend.key.height = unit(0.5,"cm"),
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

#--- draw chazhi heatmap

drawchazhiGNGheatmap <- function(data,min_va,max_va){
  ggplot(data=data, aes(x=ID1, y=ID2, fill = interaction)) + geom_raster() +
    scale_fill_gradientn(colors = c("DarkCyan", "PaleTurquoise","LightCyan","white", "Linen", "PeachPuff", "Red"),
                         values = c(0,0.49,0.499,0.5,0.501,0.51,1),limits=c(min_va,max_va)
    ) + 
    scale_y_discrete(expand=c(0,0),limits=c("GR","NG")) +
    scale_x_discrete(expand=c(0,0),position = "top",limits=c("NG","GR")) +
    theme(plot.title=element_text(hjust=0.5),
          panel.grid=element_blank(), panel.background=element_blank()) +
    theme(legend.position="bottom",
          title = element_blank() ,
          legend.key.width = unit(2,"line"),
          legend.key.height = unit(0.5,"cm"),
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
          #panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank()) + 
    labs(fill=NULL ) 
}

#--- draw normalization heatmap

drawnormGNGheatmap <- function(data,min_va,max_va){
  ggplot(data=data, aes(x=ID1, y=ID2, fill = log(interaction))) + geom_raster() +
    scale_fill_gradientn(colors = c("LightSalmon", "orangered2"),
                         values = c(0,1),limits=c(min_va,max_va)
    ) + 
    scale_y_discrete(expand=c(0,0),limits=c("GR","NG")) +
    scale_x_discrete(expand=c(0,0),position = "top",limits=c("NG","GR")) +
    theme(plot.title=element_text(hjust=0.5),
          panel.grid=element_blank(), panel.background=element_blank()) +
    theme(legend.position="bottom",
          title = element_blank() ,
          legend.key.width = unit(2,"line"),
          legend.key.height = unit(0.5,"cm"),
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
GNG_normaziinta_ko<- dplyr::mutate(ko_gene_nongene_matrix, interaction=interaction*10^12/(as.numeric(ID1_LEN)*as.numeric(ID2_LEN)*ko_reads))
GNG_normaziinta_wt<- dplyr::mutate(wt_gene_nongene_matrix, interaction=interaction*10^12/(as.numeric(ID1_LEN)*as.numeric(ID2_LEN)*wt_reads))
GNG_va_mi<-min(log(GNG_normaziinta_ko$interaction),log(GNG_normaziinta_wt$interaction))
GNG_va_ma<-max(log(GNG_normaziinta_ko$interaction),log(GNG_normaziinta_wt$interaction))
GNG_va_mi <- -17.5
GNG_va_ma <- -16.2

drawnormGNGheatmap(GNG_normaziinta_ko,GNG_va_mi,GNG_va_ma)
ggsave(paste0("GNG_norm_ko",Reso,filetype), width = 5,height = 5.5,dpi = 1000)
ggsave(paste0("GNG_norm_ko",Reso,filetype2), width = 5,height = 5.5,dpi = 1000)
drawnormGNGheatmap(GNG_normaziinta_wt,GNG_va_mi,GNG_va_ma)
ggsave(paste0("GNG_norm_wt",Reso,filetype), width = 5,height = 5.5,dpi = 1000)
ggsave(paste0("GNG_norm_wt",Reso,filetype2), width = 5,height = 5.5,dpi = 1000)
GNG_ko_wt_diff <- inte_diff(GNG_normaziinta_ko,GNG_normaziinta_wt)
GNG_diff_va <- max(GNG_ko_wt_diff$interaction)
GNG_diff_mi <- min(GNG_ko_wt_diff$interaction)
GNG_diff_va <- max(GNG_diff_va,1+abs(1-GNG_diff_mi))
GNG_diff_mi <- min(GNG_diff_mi,2-GNG_diff_va)
GNG_diff_mi <- 0.93 #0.95
GNG_diff_va <- 1.07 #1.05
drawdiffGNGheatmap(GNG_ko_wt_diff,GNG_diff_mi,GNG_diff_va)
ggsave(paste0("GNG_norm_diff",Reso,filetype), width = 5,height = 5.5,dpi = 1000)
ggsave(paste0("GNG_norm_diff",Reso,filetype2), width = 5,height = 5.5,dpi = 1000)
GNG_ko_wt_chazhi <- chazhi_diff(GNG_normaziinta_ko,GNG_normaziinta_wt)
GNG_chazhi_va <- max(GNG_ko_wt_chazhi$interaction)
GNG_chazhi_mi <- min(GNG_ko_wt_chazhi$interaction)
GNG_chazhi_va <- max(GNG_chazhi_va,abs(GNG_chazhi_mi))
GNG_chazhi_mi <- min(-GNG_chazhi_va,-abs(GNG_chazhi_mi))
GNG_chazhi_mi <- -5.5e-09
GNG_chazhi_va <- 5.5e-09
drawchazhiGNGheatmap(GNG_ko_wt_chazhi,GNG_chazhi_mi,GNG_chazhi_va)
ggsave(paste0("GNG_norm_chazhi",Reso,filetype), width = 5,height = 5.5,dpi = 1000)
ggsave(paste0("GNG_chazhi",Reso,filetype2), width = 5,height = 5.5,dpi = 1000)

rm(list=ls())











# unnormaziinta_ko <- read.table(paste0("interactionnomalize_ko",Reso),stringsAsFactors = F)
# unnormaziinta_wt <- read.table(paste0("interactionnomalize_wt",Reso),stringsAsFactors = F)
# colnames(unnormaziinta_ko) <- c("ID1","ID2","interaction","ID1_LEN","ID2_LEN")
# colnames(unnormaziinta_wt) <- c("ID1","ID2","interaction","ID1_LEN","ID2_LEN")
# unnormaziinta_ko<- dplyr::mutate(unnormaziinta_ko, interaction=interaction*10^6/(ko_reads))
# unnormaziinta_wt<- dplyr::mutate(unnormaziinta_wt, interaction=interaction*10^6/(wt_reads))
# va_mi<-min(log(unnormaziinta_ko$interaction),log(unnormaziinta_wt$interaction))
# va_ma<-max(log(unnormaziinta_ko$interaction),log(unnormaziinta_wt$interaction))
# 
# drawnormheatmap(unnormaziinta_ko,va_mi,va_ma)
# ggsave(paste0("unnorm_ko",Reso,filetype), width = 5,height = 5.5,dpi = 1000)
# ggsave(paste0("unnorm_ko",Reso,filetype2), width = 5,height = 5.5,dpi = 1000)
# drawnormheatmap(unnormaziinta_wt,va_mi,va_ma)
# ggsave(paste0("unnorm_wt",Reso,filetype), width = 5,height = 5.5,dpi = 1000)
# ggsave(paste0("unnorm_wt",Reso,filetype2), width = 5,height = 5.5,dpi = 1000)
# unko_wt_diff <- inte_diff(unnormaziinta_ko,unnormaziinta_wt)
# diff_va <- max(unko_wt_diff$interaction)
# diff_mi <- min(unko_wt_diff$interaction)
# diff_va <- max(diff_va,1+abs(1-diff_mi))
# diff_mi <- min(diff_mi,2-diff_va)
# drawdiffheatmap(unko_wt_diff,diff_mi,diff_va)
# ggsave(paste0("unnorm_diff",Reso,filetype), width = 5,height = 5.5,dpi = 1000)
# ggsave(paste0("unnorm_diff",Reso,filetype2), width = 5,height = 5.5,dpi = 1000)
# ko_wt_chazhi <- chazhi_diff(unnormaziinta_ko,unnormaziinta_wt)
# chazhi_va <- max(ko_wt_chazhi$interaction)
# chazhi_mi <- min(ko_wt_chazhi$interaction)
# chazhi_va <- max(chazhi_va,abs(chazhi_mi))
# drawchazhiheatmap(ko_wt_chazhi,chazhi_mi,chazhi_va)
# ggsave(paste0("unnorm_chazhi",Reso,filetype), width = 5,height = 5.5,dpi = 1000)
# ggsave(paste0("unnorm_chazhi",Reso,filetype2), width = 5,height = 5.5,dpi = 1000)
# rm(list=ls())
# 
