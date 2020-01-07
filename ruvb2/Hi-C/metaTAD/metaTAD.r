# this script is used to draw the heatmap of metaTAD
# --- this is the example, used the chr6.
  #--- read the file
  setwd("/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/merge/meta_TAD/juicer/wt_tad")
  setwd("/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/merge/meta_TAD/juicer/ko_tad")
  setwd('/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/merge/meta_TAD/juicer/wt_tad/oe')
  setwd('/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/merge/meta_TAD/juicer/ko_tad/oe')
  library(tidyverse)
  library(ggplot2)
  library(raster)
  library(HiCcompare)
  transform_resolution <- function(dense, res){
    r <- raster(nrow=sqrt(nrow(dense)), ncol=sqrt(nrow(dense)))
    dense <- dense[order(dense$bin1,dense$bin2),,]
    r[] <- dense$interaction
    s <- raster(nrow=res, ncol=res)
    s <- resample(r, s, method="bilinear")
    newdense <- data.frame(bin1=rep(1:res, each=res), bin2=rep(1:res, res), interaction=s[])
    return(newdense)
  }
  
  plot_tad <- function(plotdense,arg){
    ggplot(data=plotdense, aes(x=bin1, y=bin2, fill=log(interaction+1)))+geom_raster()+
      ggtitle(paste0("Averaged Interaction within TAD\n(Ruvbl2 ", arg, " merge)"))+
      scale_fill_distiller(palette = "RdYlBu",direction = -1 ,values=c(0,0.45,0.5,1),limit=c(0.5,1.3)
                           )+labs(fill="log(interaction+1)")+
      scale_x_continuous(expand=c(0, 0))+
      scale_y_reverse(expand=c(0, 0))+
      #scale_fill_continuous()+
      theme(plot.title=element_text(size=40,hjust=0.5), panel.border = element_rect(colour="black",  fill=NA, size=1), 
            axis.title=element_text(size=30),axis.text=element_text(size=28),legend.title=element_text(size=25),legend.text=element_text(size=22))+
      xlab("Position relative to TAD center (kb)")+ylab("Position relative to TAD center (kb)")
  }
  drawheatmap <- function(data){
    ggplot(data=data, aes(x=bin1, y=bin2, fill = log(interaction+1))) + geom_raster() +
      #scale_fill_gradient2(low="dodgerblue2", mid="white", high="orangered2",midpoint=2) +
      scale_fill_gradientn(colors = c("navy","royalblue4", "white", "orangered", "orangered2"),
                           values = c(0,0.2,0.3,0.4,1)#, limits=c(0.25,1.75)
      ) +
      # scale_fill_manual(breaks=c("0","0.5","1.5","7.5"),
      #                   values=c("dodgerblue", "white", "orange", "orangered")) +
      scale_y_reverse(expand=c(0,0)) +
      scale_x_continuous(expand=c(0,0)) +
      theme(plot.title=element_text(hjust=0.5),
            #panel.border = element_rect(colour="black",  fill=NA, size=1),
            panel.grid=element_blank(), panel.background=element_blank()) +
      xlab("bin 25k") +
      ylab("bin 25k") + theme(axis.line=element_blank(),
                              axis.text.x=element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks=element_blank(),
                              # axis.title.x=element_blank(),
                              # axis.title.y=element_blank(),
                              #legend.position="none",
                              panel.background=element_blank(),
                              panel.border=element_blank(),
                              panel.grid.major=element_blank(),
                              panel.grid.minor=element_blank(),
                              plot.background=element_blank())
  }
  sparse2dense<-function(mat, st, en){
    colnames(mat) <- c("bin1", "bin2", "interaction")
    dense <- data.frame(bin1=rep(st:en, each=(en-st+1)),
                        bin2=rep(st:en, en-st+1))
    sparse <- mat
    sparse <- sparse[order(sparse$bin1,sparse$bin2),,]
    sparse <- unique(sparse)
    dense <- merge(sparse, dense, by=c("bin1", "bin2"), all=T)
    dense$interaction[is.na(dense$interaction)] <- 0
    lower <- dense[order(dense$bin2), ]
    dense$interaction <- apply(cbind(dense$interaction, lower$interaction), 1, max)
    dense$chr <- unique(dense$chr)[1]
    return(dense)
  }
  tad_file <- read_tsv("tad6.tad",col_names = F)
  colnames(tad_file) <- c("bin1","bin2","interaction")
  tad_file$bin1 <- tad_file$bin1/25000
  tad_file$bin2 <- tad_file$bin2/25000
  st_bin <- min(tad_file$bin1,tad_file$bin2)
  en_bin <- max(tad_file$bin1,tad_file$bin2)
  tad_file <- sparse2dense(tad_file,st_bin,en_bin)
  #tad_file <- sparse2full(tad_file)
  drawheatmap(tad_file)
  #--- for all the tad
  tad_n <- 1744
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
  
  wt_data_oe <- data_oe
  ko_data_oe <- data_oe
  diff_data_oe <- ko_data_oe
  diff_data_oe$interaction <- ko_data_oe$interaction - wt_data_oe$interaction
  plot_diff_tad <- function(plotdense,arg){
    ggplot(data=plotdense, aes(x=bin1, y=bin2, fill=interaction))+geom_raster()+
      ggtitle(paste0("Averaged Interaction within TAD\n(Ruvbl2 ", arg, " merge)"))+
      scale_fill_gradientn(colors = c("royalblue4", "white",  "orangered2"),
                            values = c(0,0.5,1),limits = c(-0.28,0.28)
      )+labs(fill="    interaction   ")+
      scale_x_continuous(expand=c(0, 0))+
      scale_y_reverse(expand=c(0, 0))+
      #scale_fill_continuous()+
      theme(plot.title=element_text(size=40,hjust=0.5), panel.border = element_rect(colour="black",  fill=NA, size=1), 
            axis.title=element_text(size=30),axis.text=element_text(size=28),legend.title=element_text(size=25),legend.text=element_text(size=22))+
      xlab("Position relative to TAD center (kb)")+ylab("Position relative to TAD center (kb)")
  }
  plot_diff_tad(diff_data_oe,'diff')
  setwd('/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/merge/meta_TAD/')
  ggsave('oe_diff_tad_ruvbl2.pdf',width = 10.5,height = 8)
  plot_tad(wt_data_oe,'wt')
  ggsave('oe_wt_tad_ruvbl2.pdf',width = 10.5,height = 8)
  plot_tad(ko_data_oe,'ko')
  ggsave('oe_ko_tad_ruvbl2.pdf',width = 10.5,height = 8)
  