#--- to draw heatmap for product by interaction_fined.sh
setwd('/DATA/work/lbyybl/ypjiang/geneloop/pol1_loop/result/graph')
#--- only GNG plot
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)

#----------------------------------------
# read data
#----------------------------------------
pol1_num <- 75000000
pol2_num <- 75000000
pol3_num <- 75000000
CTCF_num <- 111643238
YY1_num <- 49907498
pol1_wt <- '/DATA/work/lbyybl/ypjiang/geneloop/pol1_loop/result/pol1_wt'
pol1_ko <- '/DATA/work/lbyybl/ypjiang/geneloop/pol1_loop/result/pol1_ko'
pol2_wt <- '/DATA/work/lbyybl/ypjiang/geneloop/pol2_loop/result/pol2_wt'
pol2_ko <- '/DATA/work/lbyybl/ypjiang/geneloop/pol2_loop/result/pol2_ko'
pol3_wt <- '/DATA/work/lbyybl/ypjiang/geneloop/pol3_loop/result/pol3_wt'
pol3_ko <- '/DATA/work/lbyybl/ypjiang/geneloop/pol3_loop/result/pol3_ko'
CTCF_wt <- '/DATA/work/lbyybl/ypjiang/geneloop/CTCF_loop/result/CTCF_wt'
CTCF_ko <- '/DATA/work/lbyybl/ypjiang/geneloop/CTCF_loop/result/CTCF_ko'
YY1_wt <- '/DATA/work/lbyybl/ypjiang/geneloop/YY1_loop/result/YY1_wt'
YY1_ko <- '/DATA/work/lbyybl/ypjiang/geneloop/YY1_loop/result/YY1_ko'
gene_bed <- '/DATA/work/lbyybl/genomes/mm10/DNA_elements_made_by_Boyuan/final_active_gene.bed'
gene_bed_file <- fread(gene_bed,col.names = c('chr','st','en','strand','ID'))
#--- stastic gene region length using for normalization
region_length <- gene_bed_file %>% 
  mutate(ID=str_sub(ID,1,4)) %>%
  group_by(ID) %>%
  summarise(length=sum(en-st))

#--- a function to got format data
gotnormmatrix <- function(data,num){
  data_file <- fread(data,col.names = c('loc1','loc2','num','chr1','chr2'))
  data <- data_file %>%
    mutate(class=paste0(str_sub(loc1,1,4),str_sub(loc2,1,4))) %>%
    group_by(class) %>%
    summarise(total_num=sum(num))
  data[2,2] <- data[2,2]+data[3,2]
  data[3,2] <- data[2,2]
  data$len1 <- ifelse(str_sub(data$class,1,4)=='Actg',region_length$length[1],region_length$length[2])
  data$len2 <- ifelse(str_sub(data$class,5,8)=='Actg',region_length$length[1],region_length$length[2])
  data$ID1 <- str_sub(data$class,1,4)
  data$ID2 <- str_sub(data$class,5,8)
  data[which(data$ID1=='Actg'),]$ID1 <- 'GR'
  data[which(data$ID1=='NONG'),]$ID1 <- 'NG'
  data[which(data$ID2=='Actg'),]$ID2 <- 'GR'
  data[which(data$ID2=='NONG'),]$ID2 <- 'NG'
  data$interaction <- data$total_num*1000000*1000000/(num*data$len1*data$len2)  
  return(data)
}

#---- a function to draw GNG plot
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
diff <- function(data_ko,data_wt){
  data_diff <- data_wt
  data_diff$interaction <- data_ko$interaction/data_wt$interaction
  return(data_diff)
}
chazhi <- function(data_ko,data_wt){
  data_chazhi <- data_wt
  data_chazhi$interaction <- data_ko$interaction/data_wt$interaction
  return(data_chazhi)
}
savefile <- function(name){
  ggsave(paste0(name,'.pdf'), width = 5,height = 5.5)
  ggsave(paste0(name,'.jpeg'), width = 5,height = 5.5)
}
pol1_wt <- gotnormmatrix(pol1_wt,pol1_num)
pol1_ko <- gotnormmatrix(pol1_ko,pol1_num)
pol2_wt <- gotnormmatrix(pol2_wt,pol2_num)
pol2_ko <- gotnormmatrix(pol2_ko,pol2_num)
pol3_wt <- gotnormmatrix(pol3_wt,pol3_num)
pol3_ko <- gotnormmatrix(pol3_ko,pol3_num)
CTCF_wt <- gotnormmatrix(CTCF_wt,CTCF_num)
CTCF_ko <- gotnormmatrix(CTCF_ko,CTCF_num)
YY1_wt <- gotnormmatrix(YY1_wt,YY1_num)
YY1_ko <- gotnormmatrix(YY1_ko,YY1_num)

#--- diff
pol1_diff <- diff(pol1_ko,pol1_wt)
pol2_diff <- diff(pol2_ko,pol2_wt)
pol3_diff <- diff(pol3_ko,pol3_wt)
CTCF_diff <- diff(CTCF_ko,CTCF_wt)
YY1_diff <- diff(YY1_ko,YY1_wt)
pol1vspol2diff <- diff(pol1_wt,pol2_wt)
pol1vspol3diff <- diff(pol1_wt,pol3_wt)
pol2vspol3diff <- diff(pol2_wt,pol3_wt)
max <- max(pol1_diff$interaction,pol2_diff$interaction,pol3_diff$interaction,
           CTCF_diff$interaction,YY1_diff$interaction,pol1vspol2diff$interaction,
           pol1vspol3diff$interaction,pol2vspol3diff $interaction)
min <- min(pol1_diff$interaction,pol2_diff$interaction,pol3_diff$interaction,
           CTCF_diff$interaction,YY1_diff$interaction,pol1vspol2diff$interaction,
           pol1vspol3diff$interaction,pol2vspol3diff $interaction)
max <- max(max,2-min)
min <- min(min,2-max)
drawdiffGNGheatmap(pol1_diff,min,max)
savefile('pol1_diff')
drawdiffGNGheatmap(pol2_diff,min,max)
savefile('pol2_diff')
drawdiffGNGheatmap(pol3_diff,min,max)
savefile('pol3_diff')
drawdiffGNGheatmap(CTCF_diff,min,max)
savefile('CTCF_diff')
drawdiffGNGheatmap(YY1_diff,min,max)
savefile('YY1_diff')
drawdiffGNGheatmap(pol1vspol2diff,min,max)
savefile('pol1vspol2diff_diff')
drawdiffGNGheatmap(pol1vspol3diff,min,max)
savefile('pol1vspol3diff_diff')
drawdiffGNGheatmap(pol2vspol3diff,min,max)
savefile('pol2vspol3diff_diff')
