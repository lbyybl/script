# it's used to call contact strength of loop
# means enrichmeant score of the median of all the loop 
library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)
library(purrr)
#--- a function to calculate enrichment score for one loop
#prefix='hichip'
#paste0(prefix,c('POLY', 'Enha', 'Supe'))
#--- calculate all enrichmeant for all loop
# setwd('/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/Contact_strength')
# hichip_wt_dir <- system('parallel echo wt/hichip_wt_{1}_{2} ::: POLY Enha Supe  :::  POLY Enha Supe',intern = T)
# hichip_ko_dir <- system('parallel echo ko/hichip_ko_{1}_{2} ::: POLY Enha Supe  :::  POLY Enha Supe',intern = T)
# setwd('/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/Contact_strength')
# oceanc_wt_dir <- system('parallel echo wt/oceanc_wt_{1}_{2} ::: POLY Enha Supe  :::  POLY Enha Supe',intern = T)
# oceanc_ko_dir <- system('parallel echo ko/oceanc_ko_{1}_{2} ::: POLY Enha Supe  :::  POLY Enha Supe',intern = T)

#ele_clas <- ' ::: Enha H3k4 Pol2 POLY Supe  :::  Enha H3k4 Pol2 POLY Supe'
ele_clas <- ' ::: POLY Enha Supe  :::  POLY Enha Supe'
#hichip_dir <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/Contact_strength/old/pd_e_hub_me3'
#oceanc_dir <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/Contact_strength/pd_e_hub_me3'
hichip_dir <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/Contact_strength/hiccups'
oceanc_dir <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/Contact_strength/hiccups'
setwd(hichip_dir)
hichip_wt_dir <- system(paste0('parallel echo wt/hichip_wt_{1}_{2}',ele_clas),intern = T)
hichip_ko_dir <- system(paste0('parallel echo ko/hichip_ko_{1}_{2}',ele_clas),intern = T)
setwd(oceanc_dir)
oceanc_wt_dir <- system(paste0('parallel echo wt/oceanc_wt_{1}_{2}',ele_clas),intern = T)
oceanc_ko_dir <- system(paste0('parallel echo ko/oceanc_ko_{1}_{2}',ele_clas),intern = T)

format_data <- function(dir){
  class <- system(paste0('parallel echo {1}_{2}',ele_clas),intern = T)
  inter_file <- data.frame("class"=class,
                           "num"=0)
  for (i in dir){
    apa <- paste0(i,'/10000/gw/measures.txt')
    apa_file <- fread(apa,col.names = c('class','num'))
    inter_file[which(inter_file$class==str_replace_all(i,".*[wt|ko]_","")),]$num <- apa_file$num[5]
    
  }
  colnames(inter_file) <- c('class','interaction')
  loc <- as.data.frame(str_split(inter_file$class, "_",simplify = T))
  colnames(loc) <- c('ID1','ID2')
  inter_file <- bind_cols(inter_file,loc)
  #ifelse(inter_file$ID1=="POLY",inter_file$ID1 <- 'PD',ifelse(inter_file$ID1=="Enha",inter_file$ID1 <- 'TE',inter_file$ID1 <- 'SE'))
  inter_file$ID2 <- str_replace_all(inter_file$ID2,c('POLY'='PD','Supe'='SE','Enha'='TE'))
  inter_file$ID1 <- str_replace_all(inter_file$ID1,c('POLY'='PD','Supe'='SE','Enha'='TE'))
  return(inter_file)
}
setwd(hichip_dir)
hichip_wt_file <- format_data(hichip_wt_dir)
hichip_ko_file <- format_data(hichip_ko_dir)
setwd(oceanc_dir)
oceanc_wt_file <- format_data(oceanc_wt_dir)
oceanc_ko_file <- format_data(oceanc_ko_dir)

max_va <- max(hichip_wt_file$interaction, hichip_ko_file$interaction,
              oceanc_wt_file$interaction, oceanc_ko_file$interaction)
min_va <- min(hichip_wt_file$interaction, hichip_ko_file$interaction,
              oceanc_wt_file$interaction, oceanc_ko_file$interaction)
max_va <- max_va*(1+0.05)
min_va <- min_va*(1-0.05)
#--- diff file

diff_file <- function(wt,ko){
  wt$interaction <- ko$interaction - wt$interaction
  return(wt)
}

hichip_diff <- diff_file(hichip_wt_file,hichip_ko_file)
oceanc_diff <- diff_file(oceanc_wt_file,oceanc_ko_file)
max_diff <- max(hichip_diff$interaction,oceanc_diff$interaction)
min_diff <- min(hichip_diff$interaction,oceanc_diff$interaction)
max_diff <- max(max_diff, abs(min_diff))
min_diff <- min(min_diff,-max_diff)
#limit <- c("PD","TE","SE",'H3k4','Pol2')
limit <- c("PD","TE","SE")
drawheatmap <- function(data,vi,va){
  ggplot(data=data, aes(x=ID1, y=ID2, fill = interaction)) + geom_raster() +
    scale_fill_gradientn(colors = c("deepskyblue1", "grey96",  "salmon"),
                         values = c(0,0.5,1),limits=c(vi,va)
    ) + 
    scale_y_discrete(expand=c(0,0),limits=limit) +
    scale_x_discrete(expand=c(0,0),position = "top",limits=rev(limit)) +
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
    labs(fill=NULL ) + 
    geom_hline(yintercept = 1.5,colour='white',size=0.5) +
    geom_hline(yintercept = 2.5,colour='white',size=0.5) +
    geom_hline(yintercept = 3.5,colour='white',size=0.5) +
    geom_hline(yintercept = 4.5,colour='white',size=0.5) +
    geom_vline(xintercept = 1.5,colour='white',size=0.5) +
    geom_vline(xintercept = 2.5,colour='white',size=0.5) +
    geom_vline(xintercept = 3.5,colour='white',size=0.5) +
    geom_vline(xintercept = 4.5,colour='white',size=0.5) 
}
drawdiffheatmap <- function(data,vi,va){
  ggplot(data=data, aes(x=ID1, y=ID2, fill = interaction)) + geom_raster() +
    scale_fill_gradientn(colors = c("deepskyblue1", "grey95", "yellow"),
                         values = c(0,0.5,1),limits=c(vi,va)
    ) + 
    scale_y_discrete(expand=c(0,0),limits=limit) +
    scale_x_discrete(expand=c(0,0),position = "top",limits=rev(limit)) +
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
    labs(fill=NULL ) + 
    geom_hline(yintercept = 0,colour='white',size=0.5) +
    geom_hline(yintercept = 1.5,colour='white',size=0.5) +
    geom_hline(yintercept = 2.5,colour='white',size=0.5) +
    geom_hline(yintercept = 3.5,colour='white',size=0.5) +
    geom_hline(yintercept = 4.5,colour='white',size=0.5) +
    geom_hline(yintercept = 5.5,colour='white',size=0.5) +
    geom_vline(xintercept = 0,colour='white',size=0.5) +
    geom_vline(xintercept = 1.5,colour='white',size=0.5) +
    geom_vline(xintercept = 2.5,colour='white',size=0.5) +
    geom_vline(xintercept = 3.5,colour='white',size=0.5) +
    geom_vline(xintercept = 4.5,colour='white',size=0.5) +
    geom_vline(xintercept = 5.5,colour='white',size=0.5)
}


savefile <- function(name){
  ggsave(paste0(name,".jpeg"), width = 5,height = 5)
  ggsave(paste0(name,".pdf"), width = 5,height = 5)
}
#--- draw heatmap
drawheatmap(hichip_wt_file,vi=min_va,va=max_va) 
savefile('hichip_wt')
drawheatmap(hichip_ko_file,vi=min_va,va=max_va) 
savefile('hichip_ko')
drawheatmap(oceanc_wt_file,vi=min_va,va=max_va) 
savefile('oceanc_wt')
drawheatmap(oceanc_ko_file,vi=min_va,va=max_va)
savefile('oceanc_ko')
drawdiffheatmap(hichip_diff,vi=min_diff,va=max_diff)
savefile('hichip_diff')
drawdiffheatmap(oceanc_diff,vi=min_diff,va=max_diff)
savefile('oceanc_diff')

