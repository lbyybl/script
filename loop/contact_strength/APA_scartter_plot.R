rm(list=ls())
graphics.off()
cat("\f")

#---------------------------------------------------
# 1. only select the pixels of the peak center, even the peak call
# extends to multipe pixels 2. loci separaed more dis > 300k 
#---------------------------------------------------
library(ggplot2)
library(data.table)
library(dplyr)
library(purrr)
#library(Sushi)
library(HiCcompare)

#--- sparse2full() 返回的是matrix可以直接索引
#--- set resolution and minum distance 300000
reso <- 10000
t <- 300000
path <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/Contact_strength/hiccups'
#--- one function to calculate apa for one chr

apa_chr <- function(matmane,chrname,se_loop){
  #--- a function to calculate the apa; input: se_loop of specific chr and the row number
  apa_calcu <- function(i){
    data_loop_one <- chr_seloop[i,]
    bin1 <- (data_loop_one$bin1-chr_min_bin)/reso + 1
    bin2 <- (data_loop_one$bin2-chr_min_bin)/reso + 1
    peak <- densechr[bin1,bin2]
    background <- densechr[(bin1+5):(bin1+10),(bin2+5):(bin2+10)]
    apa <- peak/mean(background)
    return(apa)
  }  
  
  #--- get the coorposding chr of loop file
  
  chr_seloop <- se_loop %>% 
    filter(chr1==chrname)
  if (nrow(chr_seloop)>=1){
    #--- interaction matrix
    chr <- fread(matmane,col.names = c('bin1','bin2','interaction'))
    
    #--- convet sparse to full
    densechr <- sparse2full(chr)
    
    #--- get the start bin of dense matrix
    chr_min_bin <- min(chr$bin1,chr$bin2)
    
    chr_seloop$apa <- 1:nrow(chr_seloop) %>% 
      map_dbl(apa_calcu)
    return(chr_seloop)
  }
}

genome_apa <- function(se_loop,mefix){
  band_func <- function(i){
    return(apa_chr(paste0(i,'_',mefix,'.bed'),paste0('chr',i),se_loop))
  }
  apa <- c(1:19,'X','Y') %>%
    map_dfr(band_func)
  return(apa)
}


#--- a function to read loop file and convert to bin acording to reso
read_loop <- function(name){
  se_loop <- fread(name,col.names = c('chr1','st1','en1','chr2','st2','en2'))
  se_loop <- se_loop %>%
    filter((st2+en2)/2-(st1+en1)/2 > t)
  se_loop <- se_loop %>%
    mutate(bin1=floor(((en1+st1)/2)/reso)*reso,bin2=floor(((en2+st2)/2)/reso)*reso)
  return(se_loop)
}

#--- read chr matrix and loop file
setwd('/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/Contact_strength/hiccups')
#--- loop file            
#----all_loop.bed 的格式为chr10   21082093        21084182        chr10   21099810        21103255
hiccups_se_loop <- read_loop(paste0(path,"/all_Supe_Supe.bed"))


ko_hiccups_se_apa <- genome_apa(hiccups_se_loop,'ko')
wt_hiccups_se_apa <- genome_apa(hiccups_se_loop,'wt')
all_apa <- merge(wt_hiccups_se_apa,ko_hiccups_se_apa,by=c('chr1','st1','en1','chr2','st2','en2','bin1','bin2'))
names(all_apa)[9:10] <- c('wt','ko')
all_apa <- all_apa[which(all_apa$wt < 10000 & all_apa$ko < 1000),]
all_apa$dis <- floor((all_apa$st2-all_apa$st1)/1000000)

ggplot(all_apa,aes(x=wt, y=ko, color=dis
),alpha=0.5) +
  scale_color_gradientn(colors = c("Wheat","Khaki","Gold","Orange","Coral","Salmon","OrangeRed","Red", "MediumVioletRed","DeepPink"),
                        values = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
  ) +
  geom_point(shape=15) +
  geom_abline(intercept = 0, slope = 1,linetype="dashed", size=0.2) +
  theme(plot.title=element_text(hjust=0.5)) +
  labs(x="wt APA",y="ko APA",color="") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0))+
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
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.text.x=element_text(angle=60, hjust=0,vjust = 0,
                             colour="black", size=12),
    axis.text.y=element_text(
      colour="black", size=12),
    panel.background=element_blank()
  )
ggplot(all_apa,aes(x=wt, y=ko, color=log10(bin2-bin1)
)) +
  scale_color_gradientn(colors = c("Wheat","Khaki","Gold","Orange","Coral","Salmon","OrangeRed","Red", "MediumVioletRed","DeepPink"),
                        values = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
  ) +
  geom_point(shape=15) +
  geom_abline(intercept = 0, slope = 1,linetype="dashed", size=0.2) +
  theme(plot.title=element_text(hjust=0.5)) +
  labs(x="wt APA",y="ko APA",color="") +
  scale_y_continuous(expand=c(0,0),limits = c(0,45)) +
  scale_x_continuous(expand=c(0,0),limits = c(0,45))+
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
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.text.x=element_text(angle=60, hjust=0,vjust = 0,
                             colour="black", size=12),
    axis.text.y=element_text(
      colour="black", size=12),
    panel.background=element_blank()
  )

setwd(path)
ggsave('hichp_APA_scalter_2.png',width=6.5,height = 6)
ggsave('hichp_APA_scalter_hichipper.pdf',width=6.5,height = 6)
