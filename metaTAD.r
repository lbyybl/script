# this script is used to draw the heatmap of metaTAD
# --- this is the example, used the chr6.
  #--- read the file
  setwd("/DATA/work/lbyybl/wangcl/hic/KPNA2/trim_linker/mateTAD/kpna2_ko_tad/kpna2_ko_tad")
  library(tidyverse)
  library(ggplot2)
  library(raster)
  
  plot_tad <- function(plotdense,arg){
    ggplot(data=plotdense, aes(x=bin1, y=bin2, fill=log(interaction+0.1)))+geom_raster()+
      ggtitle(paste0("Averaged Interaction within TAD\n(KPNA2 ", arg, " merge)"))+
      scale_fill_distiller(palette = "RdYlBu",direction = -1 ,values=c(0,0.45,0.5,1)#,limit=c(0,6.8),
                           )+labs(fill="log(interaction+1)")+
      scale_x_continuous(expand=c(0, 0))+
      scale_y_reverse(expand=c(0, 0))+
      theme(plot.title=element_text(size=40,hjust=0.5), panel.border = element_rect(colour="black",  fill=NA, size=1), 
            axis.title=element_text(size=30),axis.text=element_text(size=28),legend.title=element_text(size=25),legend.text=element_text(size=22))+
      xlab("Position relative to TAD center (kb)")+ylab("Position relative to TAD center (kb)")
  }

  tad_file <- read_tsv("tad6.tad",col_names = F,)
  colnames(tad_file) <- c("bin1","bin2","interaction")
  tad_file$bin1 <- tad_file$bin1/25000
  tad_file$bin2 <- tad_file$bin2/25000
  st_bin <- min(tad_file$bin1,tad_file$bin2)
  en_bin <- max(tad_file$bin1,tad_file$bin2)
  tad_file <- sparse2dense(tad_file,st_bin,en_bin)
  drawheatmap(tad_file)
  #--- for all the tad
  tad_n <- 1169
  for (i in 1:tad_n){
    file_name <- paste0("tad",i,".tad")
    assign(paste0("tad",i),read_tsv(file_name, col_names = F))
    data <- get(paste0("tad", i))
    colnames(data) <- c("bin1", "bin2", "interaction")
    assign(paste0("tad", i),data)
  }
  
  #--- change the bin number
  resolution <- 25000
  
  for (i in 1:tad_n){
    data <- get(paste0("tad",i))
    data$bin1 <- data$bin1/resolution
    data$bin2 <- data$bin2/resolution
    assign(paste0("tad", i),data)
  }
  
  #--- sparse to dense
  for (i in 1:tad_n){
    data <- get(paste0("tad",i))
    st_bin <- min(data$bin1,data$bin2)
    en_bin <- max(data$bin1,data$bin2)
    assign(paste0("tad",i,"_dense"),sparse2dense(data, st_bin, en_bin))
  }
  date()
  
  
  #--- change the resolution
  reso_res <- 160
  
  for (i in 1:tad_n){
    data <- get(paste0("tad",i,"_dense"))
    assign(paste0("tad",i,"_cR"),transform_resolution(data,reso_res))
  }
  date()
  
  #--- average
  data_oe <- data.frame(bin1=tad1_cR$bin1, bin2=tad1_cR$bin2,
                        interaction=c(rep(0,nrow(tad1_cR))))
  for (i in 1:tad_n){
    data <- get(paste0("tad",i,"_cR"))
    data_oe$interaction <- data_oe$interaction + data$interaction
  }
  data_oe$interaction <- data_oe$interaction/tad_n
  
  #--- heatmap
  drawheatmap(data_oe)
  plot_tad(data_oe,"ko")
  