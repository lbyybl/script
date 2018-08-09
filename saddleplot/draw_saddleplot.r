#/bin/R
# 2018-8-9
# this is used to draw compartment saddle plot, and this is used to draw according to the file produced by front steps.
#--- change the ob/ex to dense and change store it in the raster
library(getopt)
command=matrix(c("Path","p",1,"character",
				"Output","o",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)

args=getopt(command)

if (!is.null(args$help) || is.null(args$Path) || is.null(args$Output)) {

    cat(paste(getopt(command, usage = T), "\n"))

    q()

}

#--- function was used 
library(ggplot2)
library(raster)
sparse2dense<-function(mat, st, en){
  colnames(mat) <- c("bin1", "bin2", "interaction")
  dense <- data.frame(bin1=rep(st:en, each=(en-st+1)), 
                      bin2=rep(st:en, en-st+1))
  #loc <- which(mat$bin1>=st&mat$bin1<=en&mat$bin2>=st&mat$bin2<=en)
  sparse <- mat #[loc,]
  sparse <- sparse[order(sparse$bin1,sparse$bin2),,]
  dense <- merge(sparse, dense, by=c("bin1", "bin2"), all=T)
  dense$interaction[is.na(dense$interaction)] <- 0
  lower <- dense[order(dense$bin2), ]
  dense$interaction <- apply(cbind(dense$interaction, lower$interaction), 1, max)
  dense$chr <- unique(dense$chr)[1]
  return(dense)
}

transform_resolution <- function(dense, res){
  r <- raster(nrow=sqrt(nrow(dense)), ncol=sqrt(nrow(dense)))
  dense <- dense[order(dense$bin1,dense$bin2),,]
  r[] <- dense$interaction
  s <- raster(nrow=res, ncol=res)
  s <- resample(r, s, method="bilinear")
  newdense <- data.frame(bin1=rep(1:res, each=res), bin2=rep(1:res, res), interaction=s[])
  return(newdense)
}

drawheatmap <- function(data){
  ggplot(data=data, aes(x=bin1, y=bin2, fill = log(interaction+1))) + geom_raster() + 
    #scale_fill_gradient2(low="dodgerblue2", mid="white", high="orangered2",midpoint=2) +
    scale_fill_gradientn(colors = c("navy","royalblue4", "white", "orangered", "orangered2"),
                         values = c(0,0.3,0.43,0.5,1)
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

#--- rm the line with non-zero value less than 20
remove_sparse_line <- function(data){
  data_rm_sparse <- data %>%
    group_by(bin1) %>%
    filter(sum(interaction!=0)>20) 
  bin_sect <- intersect(unique(data_rm_sparse$bin2), unique(data_rm_sparse$bin1))
  data_rm_sparse <- data_rm_sparse %>%
    filter(bin2 %in% bin_sect)
  return(data_rm_sparse)
}

library(tidyverse)
setwd(args$Path)
#-- read the ob/ex file
chrname=c(1:19,"X","Y")
for (i in chrname){
  file_name <- paste0("chr",i,".oe.resort")
  assign(paste0("chr",i,"_oe"),read_csv(file_name, col_names = T))
  data <- get(paste0("chr", i, "_oe"))
  colnames(data) <- c("bin1", "bin2", "interaction")
  #data <- data[-which(data$bin1==data$bin2),,] # remove self interaction
  assign(paste0("chr", i, "_oe"),data)
}

date()

for (i in chrname){
  data <- get(paste0("chr",i,"_oe"))
  data_rm_sparse <- data %>%
    group_by(bin1) %>%
    filter(sum(interaction!=0)>20) 
  bin_sect <- intersect(unique(data_rm_sparse$bin2), unique(data_rm_sparse$bin1))
  data_rm_sparse <- data_rm_sparse %>%
    filter(bin2 %in% bin_sect)
  assign(paste0("chr",i,"_dense_rm_zero"),data_rm_sparse)
}

date()
#-- rm the blank line

remove_blank_line <- function(data){
  bin_min <- min(data$bin1,data$bin2)
  bin_max <- max(data$bin1,data$bin2)
  data <- data[order(data$bin1,data$bin2),,]
  bin_number <- length(unique(data$bin1))
  data_rm_blank <- data.frame(bin1=rep(1:bin_number, each=bin_number), 
                              bin2=rep(1:bin_number, bin_number))
  data <- data[order(data$bin1,data$bin2),,]
  data_rm_blank$interaction <- data$interaction
  return(data_rm_blank)
}

for (i in chrname){
  data <- get(paste0("chr", i, "_dense_rm_zero"))
  assign(paste0("chr",i,"_dense_rm_blank"), remove_blank_line(data))
}

date()
#--- change the resolution to 80 or 122

reso_res <- 300

for (i in chrname){
  data <- get(paste0("chr",i,"_dense_rm_blank"))
  assign(paste0("chr",i,"_cR_oe"),transform_resolution(data,reso_res))
}
date()

#-----get expect

data_oe <- data.frame(bin1=chr1_cR_oe$bin1, bin2=chr1_cR_oe$bin2,
                      interaction=c(rep(0,nrow(chr1_cR_oe))))
for (i in chrname){
  data <- get(paste0("chr",i,"_cR_oe"))
  data_oe$interaction <- data_oe$interaction + data$interaction
}
data_oe$interaction <- data_oe$interaction/21

pdf(args$Output)
drawheatmap(data_oe)
dev.off()

    