setwd("/DATA/work/lbyybl/wangcl/hic/KPNA2/trim_linker/DoStreAna/Compartment/Resulation500k/ploter/pol2")
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

chr1_data <- read.table("chr1_resorter.matrix", header = F)
#mat <- read.table("chr1_resorter.matrix", header = F)
changeresolution <- function(data,st,en,res){
  data_dense <- sparse2dense(data,st,en)
  data_changeresolution <- transform_resolution(data_dense, res)
  return(data_changeresolution)
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

# chr1_dense <- sparse2dense(chr1_data,1,391)
# chr1_changesolution <- transform_resolution(chr1_dense, 80)
drawheatmap(chr1_changesolution)

chr_changeresolution <- changeresolution(chr1_data,1,391,80)

drawheatmap(chr_changeresolution)

for (i in 1:19){
  file_name <- paste0("chr",i,"_resorter.matrix")
  assign(paste0("chr",i,"_data"),read.table(file_name, header = F)) 
}

# find max bin and min bin 

for (i in 1:19){
  file_name <- paste0("chr",i,".bed")
  assign(paste0("chr",i,"_bed"),read.table(file_name, header = F))
  data_bed <- get(paste0("chr",i,"_bed"))
  assign(paste0("chr",i,"_max"),max(data_bed[,4]))
  assign(paste0("chr",i,"_min"),min(data_bed[,4]))
}

for (i in 1:19){
  data <- get(paste0("chr",i,"_data"))
  st_bin <- get(paste0("chr",i,"_min"))
  en_bin <- get(paste0("chr",i,"_max"))
  assign(paste0("chr",i,"_changeresolution"),changeresolution(data,1,en_bin-st_bin+1,80))
}

sum_chngere <- data.frame(bin1=chr1_changeresolution$bin1, bin2=chr1_changeresolution$bin2,
                          interaction=c(rep(0,nrow(chr10_changeresolution))))
for (i in 1:19){
  data <- get(paste0("chr",i,"_changeresolution"))
  sum_chngere$interaction <- sum_chngere$interaction + data$interaction
}

sum_chngere$interaction <- sum_chngere$interaction/19

drawheatmap(sum_chngere)
drawheatmap(chr17_changeresolution)
mean_dignal <- mean(sum_chngere[(sum_chngere$bin1==sum_chngere$bin2),3])
sum_chngere$interaction <- sum_chngere$interaction/mean_dignal
