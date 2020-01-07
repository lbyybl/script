setwd("/DATA/work/lbyybl/ypjiang/4c/4c20181204/rawdata/Hspa8-Nla3-Pol2_TKD181102142/test/all_bam")
library(Gviz)
library(GenomicRanges)
# data(cpgIslands)
# class(cpgIslands)
# chr <- as.character(unique(seqnames(cpgIslands)))
# gen <- genome(cpgIslands)
# atrack <- AnnotationTrack(cpgIslands, name="CpG")
# itrack <- IdeogramTrack(genome=gen, chromosome=chr, showTitle=TRUE,
#                         col.title="black",fontcolor="black",col="red",
#                         alpha=1, showId = FALSE,bevel=1,lty=0,
#                         rotation.title=0,lwd=0)
# plotTracks(list(gtrack, atrack, itrack))

#---alignment track

afrom <- 23679420
ato <- 23679440
# alTrack <- AlignmentsTrack(system.file(package="Gviz", "extdata", "gapped.bam"), isPaired=TRUE,
#                            col="red", col.coverage="green",lwd.sashimiMax=300,
#                            col.axis="black",name="abc",col.title="black",fill="green")
# plotTracks(c(bmt, alTrack), from=afrom, to=ato, chromosome="chr12", min.height=0, coverageHeight=0.08,
#            minCoverageHeight=0,type = "coverage")
# plotTracks(c(alTrack, itrack), from=afrom, to=ato, chromosome="chr12",
#            type = "coverage",background.panel = "white", background.title = "white",
#            col.main = "balck")


#---mydata 
getaligntrack <- function(filename,name){
  alTrack <- AlignmentsTrack(filename, isPaired=FALSE,
                             col="white", col.coverage="#9999CC",
                             col.axis="black",name=name,col.title="black",fill="#9999CC",
                             coverageHeight = 0.01, min.height = 25,max.height=30,
                             heights = c(0, 100))
  return(alTrack)
}

getitrack <- function(chr){
  itrack <- IdeogramTrack(genome="hg19", chromosome=chr, showTitle=TRUE,
                          col.title="black",fontcolor="black",col="white",
                          alpha=1, showId = FALSE,bevel=1,lty=0,fill="white",
                          rotation.title=0,lwd=0)
  return(itrack)
}


gtrack <- GenomeAxisTrack()
alTrack <- getaligntrack("AD38withDOX_local.rmdump.bam","AD38+DOX")
itrack <- getitrack("chr21")
plotTracks(itrack,from = afrom,to=ato,chromosome = "chr21")

pdf("test3.pdf",width = 6,height = 5)
plotTracks(c(gtrack,gtrack,gtrack,gtrack,gtrack,alTrack, itrack), 
           from=afrom, to=ato, chromosome="chr21",
           background.panel = "white", background.title = "white",type = "coverage",
           col.main = "balck")
dev.off()
#---plot chrom
plotchr <- function(chr){
  itrack <- getitrack(chr)
  pdf(paste0(chr,".pdf"),width = 6,height = 5)
  plotTracks(c(gtrack,gtrack,gtrack,gtrack,gtrack, itrack), 
             from=afrom, to=ato, chromosome=chr,
             background.panel = "white", background.title = "white",type = "coverage",
             col.main = "balck")
  dev.off()
  png(paste0(chr,".png"),width = 600,height = 500)
  plotTracks(c(gtrack,gtrack,gtrack,gtrack,gtrack, itrack), 
             from=afrom, to=ato, chromosome=chr,
             background.panel = "white", background.title = "white",type = "coverage",
             col.main = "balck")
  dev.off()
}
for (i in c(1:22,"X","Y")){
  plotchr(paste0("chr",i))
}
itrack <- getitrack("chr21")
plotTracks(c(gtrack,gtrack,gtrack,gtrack,gtrack, itrack), 
           from=afrom, to=ato, chromosome="chr21",
           background.panel = "white", background.title = "white",type = "coverage",
           col.main = "balck")


# plotTracks(c(alTrack, itrack), from=1, to=133851895, chromosome="chr12",
#            type = "coverage",background.panel = "white", background.title = "white",
#            col.main = "balck")

#---chr
plotchr <- function(chr){
  itrack1 <- getitrack(chr)
  plotTracks(itrack1, from=1, to=ato, chromosome=chr,
             type = "coverage",background.panel = "white", background.title = "white",
             col.main = "balck")
}
plotchr("chr1")
plotchr("chr2")

#--- draw density manually
file1 <- "/DATA/work/lbyybl/ypjiang/4c/4c20181204/rawdata/Hspa8-Nla3-Pol2_TKD181102142/test/all_bam/bigwig/AD38withDOX_local.rmdump.bedgraph"
library(ggplot2)
library(readr)
library(dplyr)
density <- read.table(file1,header = FALSE)
colnames(density) <- c("chrom","start","end","value")
density <- read_tsv(file1,col_names=c("chrom","st","en","num"))
head(density)
start <- density[,c(1,2,4)]
colnames(start) <- c("chrom","loc","num")
end <- density[,c(1,3,4)]
colnames(end) <- c("chrom","loc","num")
total<-rbind(start,end)
chr21 <- total %>%
  filter(chrom=="chr21")
ggplot(chr21,aes(x=loc,y=num))+geom_line()


#---- draw track with Sushi
library(Sushi)
library(dplyr)
options(SweaveHooks=list(fig=function() par(mgp=c(3, .4, 0))))
options(continue=" ")

AD38withDOX <- "/DATA/work/lbyybl/ypjiang/4c/4c20181204/rawdata/Hspa8-Nla3-Pol2_TKD181102142/test/all_bam/bigwig/AD38withDOX_local.rmdump.bedgraph"
AD38A <- "/DATA/work/lbyybl/ypjiang/4c/4c20181204/rawdata/Hspa8-Nla3-Pol2_TKD181102142/test/all_bam/bigwig/AD38A_local.rmdump.bedgraph"
AD38noDOX11d <- "/DATA/work/lbyybl/ypjiang/4c/4c20181204/rawdata/Hspa8-Nla3-Pol2_TKD181102142/test/all_bam/bigwig/AD38noDOX11d_local.rmdump.bedgraph"
AD38noDOX6d <- "/DATA/work/lbyybl/ypjiang/4c/4c20181204/rawdata/Hspa8-Nla3-Pol2_TKD181102142/test/all_bam/bigwig/AD38noDOX6d_local.rmdump.bedgraph"
AC12withHBV <- "/DATA/work/lbyybl/ypjiang/4c/4c20181204/rawdata/Hspa8-Nla3-Pol2_TKD181102142/test/all_bam/bigwig/AC12withHBV_local.rmdump.bedgraph"
DE19withDOX <- "/DATA/work/lbyybl/ypjiang/4c/4c20181204/rawdata/Hspa8-Nla3-Pol2_TKD181102142/test/all_bam/bigwig/DE19withDOX_local.rmdump.bedgraph"

readfile <- function(file){
  density <- read.table(file,header = FALSE)
  colnames(density) <- c("chrom","start","end","value")
  return(density)
}

AD38withDOX_file <- readfile(AD38withDOX)
AD38A_file <- readfile(AD38A)
AD38noDOX11d_file <- readfile(AD38noDOX11d)
AD38noDOX6d_file <- readfile(AD38noDOX6d)
AC12withHBV_file <- readfile(AC12withHBV)
DE19withDOX_file <- readfile(DE19withDOX)
findsubmax <- function(data){
  sub_data <- data %>%
    filter(chrom=="chr2" | chrom=="chr8" | chrom=="chr21" | chrom=="chr22")
  return(max(sub_data$value))
}


getOption("SweaveHooks")[["fig"]]()
# ma_va <- max(AD38withDOX_file$value,AD38A_file$value,AD38noDOX11d_file$value,
#              AD38noDOX6d_file$value,AC12withHBV_file$value,DE19withDOX_file$value)
# ma_round <- round(ma_va/100)*100
# lables <- seq(0,ma_round,ma_round/5)
drawwhochr <- function(file,chr,name,ma_round,lables){
  if (chr=="chr1"){
    chrom            = "chr1"
    chromstart       = 1
    chromend         = 249250621
  }else if (chr=="chr2"){
    chrom            = "chr2"
    chromstart       = 1
    chromend         = 243199373
  }else if (chr=="chr3"){
    chrom            = "chr3"
    chromstart       = 1
    chromend         = 198022430
  }else if (chr=="chr4"){
    chrom            = "chr4"
    chromstart       = 1
    chromend         = 191154276
  }else if (chr=="chr5"){
    chrom            = "chr5"
    chromstart       = 1
    chromend         = 180915260
  }else if (chr=="chr6"){
    chrom            = "chr6"
    chromstart       = 1
    chromend         = 171115067
  }else if (chr=="chr7"){
    chrom            = "chr7"
    chromstart       = 1
    chromend         = 159138663
  }else if (chr=="chr8"){
    chrom            = "chr8"
    chromstart       = 1
    chromend         = 146364022
  }else if (chr=="chr9"){
    chrom            = "chr9"
    chromstart       = 1
    chromend         = 141213431
  }else if (chr=="chr10"){
    chrom            = "chr10"
    chromstart       = 1
    chromend         = 135534747
  }else if (chr=="chr11"){
    chrom            = "chr11"
    chromstart       = 1
    chromend         = 135006516
  }else if (chr=="chr12"){
    chrom            = "chr12"
    chromstart       = 1
    chromend         = 133851895
  }else if (chr=="chr13"){
    chrom            = "chr13"
    chromstart       = 1
    chromend         = 115169878
  }else if (chr=="chr14"){
    chrom            = "chr14"
    chromstart       = 1
    chromend         = 107349540
  }else if (chr=="chr15"){
    chrom            = "chr15"
    chromstart       = 1
    chromend         = 102531392
  }else if (chr=="chr16"){
    chrom            = "chr16"
    chromstart       = 1
    chromend         = 90354753
  }else if (chr=="chr17"){
    chrom            = "chr17"
    chromstart       = 1
    chromend         = 81195210
  }else if (chr=="chr18"){
    chrom            = "chr18"
    chromstart       = 1
    chromend         = 78077248
  }else if (chr=="chr19"){
    chrom            = "chr19"
    chromstart       = 1
    chromend         = 59128983
  }else if (chr=="chr20"){
    chrom            = "chr20"
    chromstart       = 1
    chromend         = 63025520
  }else if (chr=="chr21"){
    chrom            = "chr21"
    chromstart       = 1
    chromend         = 48129895
  }else if (chr=="chr22"){
    chrom            = "chr22"
    chromstart       = 1
    chromend         = 51304566
  }else if (chr=="chrX"){
    chrom            = "chrX"
    chromstart       = 1
    chromend         = 155270560
  }else if (chr=="chrY"){
    chrom            = "chrY"
    chromstart       = 1
    chromend         = 59373566
  }
  #pdf(paste0(name,chr,".pdf"),width = 6,height = 3)
  png(paste0(name,chr,".png"),width = 1000,height = 300)
  plotBedgraph(file,chrom,chromstart,
               chromend,color= "orange3",range = c(0,ma_round))
  axis(side = 2,lables,lables,cex.axis=1,las=1)
  dev.off()
  pdf(paste0(name,chr,".pdf"),width = 6,height = 3)
  plotBedgraph(file,chrom,chromstart,
               chromend,color= "orange3",range = c(0,ma_round))
  axis(side = 2,lables,lables,cex.axis=1,las=1)
  dev.off()
}
drawlocalchr <- function(file,chr,name,ma_round,lables){
  if (chr=="chr2"){
    chrom            = "chr2"
    chromstart       = 27232670
    chromend         = 28071452
  }else if (chr=="chr8"){
    chrom            = "chr8"
    chromstart       = 27232670
    chromend         = 28071452
  }else if (chr=="chr21"){
    chrom            = "chr21"
    chromstart       = 23679420
    chromend         = 24160718
  }else if (chr=="chr22"){
    chrom            = "chr22"
    chromstart       = 29174883
    chromend         = 29554916
  }
  #pdf(paste0(name,chr,".pdf"),width = 6,height = 3)
  png(paste0(name,chr,".png"),width = 1000,height = 300)
  plotBedgraph(file,chrom,chromstart,
               chromend,color= "orange3",range = c(0,ma_round))
  axis(side = 2,lables,lables,cex.axis=1,las=1)
  dev.off()
  pdf(paste0(name,chr,".pdf"),width = 6,height = 3)
  plotBedgraph(file,chrom,chromstart,
               chromend,color= "orange3",range = c(0,ma_round))
  axis(side = 2,lables,lables,cex.axis=1,las=1)
  dev.off()
}
setwd("/DATA/work/lbyybl/ypjiang/4c/4c20181204/rawdata/Hspa8-Nla3-Pol2_TKD181102142/test/all_bam/grahp")
draw4 <- function(file,name){
  ma_va <- findsubmax(file) # ma_va <- findsubmax(AD38withDOX_file)
  ma_round <- round(ma_va/100)*100
  lables <- seq(0,ma_round,ma_round/5)
  drawwhochr(file,"chr1",name,ma_round,lables)
  drawwhochr(file,"chr2",name,ma_round,lables)
  drawwhochr(file,"chr3",name,ma_round,lables)
  drawwhochr(file,"chr4",name,ma_round,lables)
  drawwhochr(file,"chr5",name,ma_round,lables)
  drawwhochr(file,"chr6",name,ma_round,lables)
  drawwhochr(file,"chr7",name,ma_round,lables)
  drawwhochr(file,"chr8",name,ma_round,lables)
  drawwhochr(file,"chr9",name,ma_round,lables)
  drawwhochr(file,"chr10",name,ma_round,lables)
  drawwhochr(file,"chr11",name,ma_round,lables)
  drawwhochr(file,"chr12",name,ma_round,lables)
  drawwhochr(file,"chr13",name,ma_round,lables)
  drawwhochr(file,"chr14",name,ma_round,lables)
  drawwhochr(file,"chr15",name,ma_round,lables)
  drawwhochr(file,"chr16",name,ma_round,lables)
  drawwhochr(file,"chr17",name,ma_round,lables)
  drawwhochr(file,"chr18",name,ma_round,lables)
  drawwhochr(file,"chr19",name,ma_round,lables)
  drawwhochr(file,"chr20",name,ma_round,lables)
  drawwhochr(file,"chr21",name,ma_round,lables)
  drawwhochr(file,"chr22",name,ma_round,lables)
  drawwhochr(file,"chrX",name,ma_round,lables)
  drawwhochr(file,"chrY",name,ma_round,lables)
}
drawlocal4 <- function(file,name){
  ma_va <- findsubmax(file)
  ma_round <- round(ma_va/100)*100
  lables <- seq(0,ma_round,ma_round/5)
  drawlocalchr(file,"chr2",name,ma_round,lables)
  drawlocalchr(file,"chr8",name,ma_round,lables)
  drawlocalchr(file,"chr21",name,ma_round,lables)
  drawlocalchr(file,"chr22",name,ma_round,lables)
}
#---local
drawlocal4(AD38withDOX_file,"AD38withDOX_loc")
drawlocal4(AD38A_file,"AD38A_loc")
drawlocal4(AD38noDOX11d_file,"AD38noDOX11d_loc")
drawlocal4(AD38noDOX6d_file,"AD38noDOX6d_loc")
drawlocal4(AC12withHBV_file,"AC12withHBV_loc")
drawlocal4(DE19withDOX_file,"DE19withDOX_loc")


#----whole
draw4(AD38withDOX_file,"AD38withDOX")
draw4(AD38A_file,"AD38A")
draw4(AD38noDOX11d_file,"AD38noDOX11d")
draw4(AD38noDOX6d_file,"AD38noDOX6d")
draw4(AC12withHBV_file,"AC12withHBV")
draw4(DE19withDOX_file,"DE19withDOX")
#---histon marker
# downsample <- function(data){
#   # downsample for plotting
#   signaltrack <- read.table(data,header = FALSE)
#   signaltrack <- signaltrack[order(signaltrack[,1],signaltrack[,2],signaltrack[,3]),]
#   while (nrow(signaltrack) > 256000)
#   {
#     # downsample for plotting if neccesary
#     if (nrow(signaltrack) %% 2 != 0)
#     {
#       signaltrack = signaltrack[1:(nrow(signaltrack)-1),]
#     }
#     chr=signaltrack[seq(1, nrow(signaltrack), 2),1]
#     starts = signaltrack[seq(1, nrow(signaltrack), 2),2]
#     stops  = signaltrack[seq(2, nrow(signaltrack), 2),3]
#     meanval = apply(cbind(signaltrack[seq(1, nrow(signaltrack), 2),4],signaltrack[seq(2, nrow(signaltrack), 2),4]), 1, mean)
#     signaltrack = data.frame("chr"=chr,"starts"=starts,
#                              "stops"=stops,"meanval"=meanval)
#   }
#   return(signaltrack)
# }
setwd("/DATA/work/lbyybl/ypjiang/4c/4c20181204/rawdata/Hspa8-Nla3-Pol2_TKD181102142/test/all_bam/correlation/encode_hepg2/bigwig_hg19/")
H3K27ac <- readfile("H3K27ac.merge.bedgraph")
colnames(H3K27ac) <- c("chrom","start","end","value")
H3K36me3 <- readfile("H3K36me3.merge.bedgraph")
colnames(H3K36me3) <- c("chrom","start","end","value")
H3K4me1 <- readfile("H3K4me1.merge.bedgraph")
colnames(H3K4me1) <- c("chrom","start","end","value")
H3K4me2 <- readfile("H3K4me2.merge.bedgraph")
colnames(H3K4me2) <- c("chrom","start","end","value")
H3K4me3 <- readfile("H3K4me3.merge.bedgraph")
colnames(H3K4me3) <- c("chrom","start","end","value")
H3K9ac <- readfile("H3K9ac.merge.bedgraph")
colnames(H3K9ac) <- c("chrom","start","end","value")
H3K9me3 <- readfile("H3K9me3.merge.bedgraph")
colnames(H3K9me3) <- c("chrom","start","end","value")
H4K20me1 <- readfile("H4K20me1.merge.bedgraph")
colnames(H4K20me1) <- c("chrom","start","end","value")
SALP_HepG2 <- readfile("SALP_HepG2.bedgraph")
colnames(SALP_HepG2) <- c("chrom","start","end","value")

#---local
drawlocal4(H3K27ac,"H3K27ac_loc")
drawlocal4(H3K36me3,"H3K36me3_loc")
drawlocal4(H3K4me1,"H3K4me1_loc")
drawlocal4(H3K4me2,"H3K4me2_loc")
drawlocal4(H3K4me3,"H3K4me3_loc")
drawlocal4(H3K9ac,"H3K9ac_loc")
drawlocal4(H3K9me3,"H3K9me3_loc")
drawlocal4(H4K20me1,"H4K20me1_loc")
drawlocal4(SALP_HepG2,"SALP_HepG2_loc")
#---whole

draw4(H3K27ac,"H3K27ac")
draw4(H3K36me3,"H3K36me3")
draw4(H3K4me1,"H3K4me1")
draw4(H3K4me2,"H3K4me2")
draw4(H3K4me3,"H3K4me3")
draw4(H3K9ac,"H3K9ac")
draw4(H3K9me3,"H3K9me3")
draw4(H4K20me1,"H4K20me1")
draw4(SALP_HepG2,"SALP_HepG2")


#----- test
drawwhochr(AD38withDOX_file,"chr21","AD38withDOX")

chrom            = "chr2"
chromstart       = 1
chromend         = 243199373

png(paste0("AD38withDOX","chr2",".png"),width = 1000,height = 300)
plotBedgraph(AD38withDOX_file,chrom,chromstart,
             chromend,colorbycol= SushiColors(5),range = c(0,round((max(AD38withDOX_file$value))/100)*100))
#round((max(density$value))/10)
seq<-seq(0,round((max(AD38withDOX_file$value))/100)*100,round((max(AD38withDOX_file$value))/100)*20)
axis(side = 2,seq,seq,cex.axis=1,las=1)
dev.off()


#--------------------------------------------------------------
# stastic the distribution of reads in diffrent chrom
#--------------------------------------------------------------
library(ggplot2)
library(easyGgplot2)
library(reshape2)
setwd("/DATA/work/lbyybl/ypjiang/4c/4c20181204/rawdata/Hspa8-Nla3-Pol2_TKD181102142/test/all_bam/chr_stas")
#---- calculate expected ration according to chr length
  chr_lenth <- "/home/boyuanli/tools/hic-pro/HiC-Pro_2.9.0/annotation/chrom_hg19.sizes"
  chr_lenth_file <- read.table(chr_lenth,header = FALSE)
  colnames(chr_lenth_file) <- c("chrom","len")
  chr_lenth_file$len_raio <- chr_lenth_file$len/sum(chr_lenth_file$len) 
#---- calculate expcted reads
  gotexpected <- function(file,chr_len_file){
    #--get expect
    data <- read.table(file,header = FALSE)
    colnames(data) <- c("chrom","observed")
    chr_len_file$expected <- data[1,2]*chr_len_file$len_raio
    #--merge
    merge_reads <- merge(data,chr_len_file, by="chrom")
    rage_ma <- max(merge_reads$observed,merge_reads$expected)
    rage_mi <- min(merge_reads$observed,merge_reads$expected)
    #---rm chrM
    merge_reads$ove <- merge_reads$observed/merge_reads$expected
    merge_reads <- merge_reads[-which(merge_reads$chrom=="chrM"),]
    return(merge_reads)
  }
  ratio_bar <- function(file){
    ggplot(file,aes(x=chrom,y=ove)) + geom_col(fill="deepskyblue4") +
      labs(x="",y="observed/expected",color="") +
      scale_x_discrete(limits=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))+
      geom_hline(yintercept = 1,linetype="dashed", size=0.2)+
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
        axis.text.x=element_text(angle=0, hjust=0.5,vjust = 0,
                                 colour="black", size=12),
        axis.text.y=element_text(#face = "bold",
          colour="black", size=12),
        panel.background=element_blank()
      )
  }
  
  dodge_bar <- function(file){
    file <- melt(file,id.var = c("chrom","len","len_raio","ove"),
                        variable.name = "class",value.name = "num")
    ggplot(file,aes(x=chrom,y=num,fill=class))+geom_col(position = "dodge2")+
      labs(x="",y="reads count",color="") +
      scale_y_continuous(expand=c(0,0)) +
      scale_fill_discrete(name="")+
      scale_x_discrete(limits=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))+
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
        axis.text.x=element_text(angle=0, hjust=0.5,vjust = 0,
                                 colour="black", size=12),
        axis.text.y=element_text(#face = "bold",
          colour="black", size=12),
        panel.background=element_blank()
      )
  }
  fill_bar <- function(file){
    file <- melt(file,id.var = c("chrom","len","len_raio","ove"),
                 variable.name = "class",value.name = "num")
    ggplot(file,aes(x=chrom,y=num,fill=class))+geom_col()+
      labs(x="",y="reads count",color="") +
      scale_y_continuous(expand=c(0,0)) +
      scale_fill_discrete(name="")+
      scale_x_discrete(limits=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))+
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
        axis.text.x=element_text(angle=0, hjust=0.5,vjust = 0,
                                 colour="black", size=12),
        axis.text.y=element_text(#face = "bold",
          colour="black", size=12),
        panel.background=element_blank()
      )
  }
  savefile <- function(name){
    ggsave(paste0(name,".png"),width=12,height = 6)
    ggsave(paste0(name,".pdf"),width=12,height = 6)
  }
  drawsample <- function(file,name){
    merge_file <- gotexpected(file,chr_lenth_file)
    write.csv(merge_file,paste0(name,"merge.csv"),quote=FALSE,row.names = FALSE)
    ratio_bar(merge_file)
    savefile(paste0(name,"ratio"))
    dodge_bar(merge_file)
    savefile(paste0(name,"dodge"))
    fill_bar(merge_file)
    savefile(paste0(name,"fill"))
  }
  drawsample("AC12withHBV.txt","AC12withHBV")
  drawsample("AD38A.txt","AD38A")
  drawsample("AD38noDOX11d.txt","AD38noDOX11d")
  drawsample("AD38noDOX6d.txt","AD38noDOX6d")
  drawsample("AD38withDOX.txt","AD38withDOX")
  drawsample("DE19withDOX.txt","DE19withDOX")
  

  
#-----------------------------------------------------------
# test
#----------------------------------------------------------- 
  AC12withHBV <- read.table("AC12withHBV.txt",header = FALSE)
  colnames(AC12withHBV) <- c("chrom","observed")
  chr_lenth_file$expected <- AC12withHBV[1,2]*chr_lenth_file$len_raio  
#---merge expeted and observed
  merge_reads <- merge(AC12withHBV,chr_lenth_file, by="chrom")
  rage_ma <- max(merge_reads$observed,merge_reads$expected)
  rage_mi <- min(merge_reads$observed,merge_reads$expected)
#--- draw observed/expected
  merge_reads$ove <- merge_reads$observed/merge_reads$expected
  merge_reads <- merge_reads[-which(merge_reads$chrom=="chrM"),]
  #--ratio bar
  ggplot(merge_reads,aes(x=chrom,y=ove)) + geom_col(fill="deepskyblue4") +
    labs(x="",y="observed/expected",color="") +
    scale_x_discrete(limits=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))+
    geom_hline(yintercept = 1,linetype="dashed", size=0.2)+
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
      axis.text.y=element_text(#face = "bold",
        colour="black", size=12),
      panel.background=element_blank()
    )
  #--- draw observed and expected
  merge_reads <- melt(merge_reads,id.var = c("chrom","len","len_raio","ove"),
                      variable.name = "class",value.name = "num")
  ggplot(merge_reads,aes(x=chrom,y=num,fill=class))+geom_col(position = "dodge2")+
    labs(x="expect reads number",y="observed reads number",color="") +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_discrete(name="")+
    scale_x_discrete(limits=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))+
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
      axis.text.y=element_text(#face = "bold",
        colour="black", size=12),
      panel.background=element_blank()
    )
  
  ggplot(merge_reads,aes(x=chrom,y=num,fill=class))+geom_col()+
    labs(x="expect reads number",y="observed reads number",color="") +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_discrete(name="")+
    scale_x_discrete(limits=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))+
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
      axis.text.y=element_text(#face = "bold",
        colour="black", size=12),
      panel.background=element_blank()
    )
  
  
  AC12withHBV_file <- gotexpected("AC12withHBV.txt",chr_lenth_file)
  ratio_bar(AC12withHBV_file)
  dodge_bar(AC12withHBV_file)
  fill_bar(AC12withHBV_file)
  ggsave("AC12withHBV.png",width=12,height = 6)
  ggsave("AC12withHBV.pdf",width=12,height = 6)
  
  
  
  # ggplot2.barplot(data=merge_reads, xName='chrom', yName="num",
  #                 groupName='class',# groupColors=c('#999999','#E69F00'),
  #                 position=position_dodge(),
  #                 #background and line colors
  #                 backgroundColor="white", color="black", 
  #                 ytitle="reads count", 
  #                 xtitle="",
  #                 #mainTitle="Total bill\n per time of day",
  #                 removePanelGrid=TRUE,removePanelBorder=TRUE,
  #                 axisLine=c(0.5, "solid", "black"),
  #                 showLegend=FALSE,
  #                 xlim=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
  # )
  # #--- ratio bar
  # ggplot2.barplot(data=merge_reads, xName="chrom", yName='ove',groupName="chrom",
  #                 backgroundColor="white", fill='dogerblue4', color="black",
  #                 removePanelGrid=TRUE,removePanelBorder=TRUE,
  #                 axisLine=c(0.5, "solid", "black"))
  # 
  # df<-restaurant[c(1,4), c("time", "total_bill")]
  # ggplot2.barplot(data=df, xName="time", yName='total_bill', 
  #                 groupName="time", 
  #                 #background and line colors
  #                 backgroundColor="white", color="black",
  #                 brewerPalette="Paired")
#--- fiail to draw point plot
  ggplot(merge_reads,aes(x=expected,y=observed)) + geom_point(shape=19,color="red") +
    geom_label(aes(label=chrom),hjust = 0.1, nudge_x = 0.1,size=0.5) +
    geom_abline(intercept = 0, slope = 1,linetype="dashed", size=0.2) +
    theme(plot.title=element_text(hjust=0.5)) +
    labs(x="expect reads number",y="observed reads number",color="") +
    scale_y_continuous(expand=c(0,0),limits = c(rage_mi-300,rage_ma+300)) +
    scale_x_continuous(expand=c(0,0),limits = c(rage_mi-300,rage_ma+300))+
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
      axis.text.y=element_text(#face = "bold",
        colour="black", size=12),
      panel.background=element_blank()
    )
  
  
  #----------------------------------------------------------------------
  # track circos plot
  #----------------------------------------------------------------------
  options(stringsAsFactors = FALSE);
  library(OmicCircos);
  setwd("/DATA/work/lbyybl/ypjiang/4c/4c20181204/rawdata/Hspa8-Nla3-Pol2_TKD181102142/test/all_bam/bigwig")
  AC12withHBV <- read.table("AC12withHBV_local.rmdump.bedgraph",header = FALSE)
  AD38A <- read.table("AD38A_local.rmdump.bedgraph",header = FALSE)
  AD38noDOX11d <- read.table("AD38noDOX11d_local.rmdump.bedgraph",header = FALSE)
  AD38noDOX6d <- read.table("AD38noDOX6d_local.rmdump.bedgraph",header = FALSE)
  AD38withDOX <- read.table("AD38withDOX_local.rmdump.bedgraph",header = FALSE)
  DE19withDOX <- read.table("DE19withDOX_local.rmdump.bedgraph",header = FALSE)
  
  downsample <- function(data){
    # downsample for plotting
    signaltrack <- read.table(data,header = FALSE)
    signaltrack <- signaltrack[order(signaltrack[,1],signaltrack[,2],signaltrack[,3]),]
    while (nrow(signaltrack) > 256000)
    {
      # downsample for plotting if neccesary
      if (nrow(signaltrack) %% 2 != 0)
      {
        signaltrack = signaltrack[1:(nrow(signaltrack)-1),]
      }
      chr=signaltrack[seq(1, nrow(signaltrack), 2),1]
      starts = signaltrack[seq(1, nrow(signaltrack), 2),2]
      stops  = signaltrack[seq(2, nrow(signaltrack), 2),3]
      meanval = apply(cbind(signaltrack[seq(1, nrow(signaltrack), 2),4],signaltrack[seq(2, nrow(signaltrack), 2),4]), 1, mean)
      signaltrack = data.frame("chr"=chr,"starts"=starts,
                               "stops"=stops,"meanval"=meanval)
    }
    return(signaltrack)
  }
  setwd("/DATA/work/lbyybl/ypjiang/4c/4c20181204/rawdata/Hspa8-Nla3-Pol2_TKD181102142/test/all_bam/correlation/encode_hepg2/bigwig_hg19/")
  H3K27ac <- downsample("H3K27ac.merge.bedgraph")
  H3K36me3 <- downsample("H3K36me3.merge.bedgraph")
  H3K4me1 <- downsample("H3K4me1.merge.bedgraph")
  H3K4me2 <- downsample("H3K4me2.merge.bedgraph")
  H3K4me3 <- downsample("H3K4me3.merge.bedgraph")
  H3K9ac <- downsample("H3K9ac.merge.bedgraph")
  H3K9me3 <- downsample("H3K9me3.merge.bedgraph")
  H4K20me1 <- downsample("H4K20me1.merge.bedgraph")
  SALP_HepG2 <- downsample("SALP_HepG2.bedgraph")
  draw_circle <- function(data,name,i){
    colors   <- rainbow(15, alpha=0.9)
    pdf(paste0(name,".pdf"),width = 12, height = 12)
    par(mar=c(0, 0, 0, 0));
    
    plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
    
    circos(R=400, type="chr", cir="hg19", col=colors, print.chr.lab=TRUE, W=4, scale=TRUE);
    circos(R=360, cir="hg19", W=40, mapping=data, col.v=4, type="l",   B=TRUE, col=colors[i], lwd=0.1, scale=TRUE);
    dev.off()
  }
  setwd("/DATA/work/lbyybl/ypjiang/4c/4c20181204/rawdata/Hspa8-Nla3-Pol2_TKD181102142/test/all_bam/correlation/encode_hepg2/bigwig_hg19/")
  draw_circle(H3K27ac,"H3K27ac",1)
  draw_circle(H3K36me3,"H3K36me3",2)
  draw_circle(H3K4me1,"H3K4me1",3)
  draw_circle(H3K4me2,"H3K4me2",4)
  draw_circle(H3K4me3,"H3K4me3",5)
  draw_circle(H3K9ac,"H3K9ac",6)
  draw_circle(H3K9me3,"H3K9me3",7)
  draw_circle(H4K20me1,"H4K20me1",8)
  draw_circle(SALP_HepG2,"SALP_HepG2",9)
  setwd("/DATA/work/lbyybl/ypjiang/4c/4c20181204/rawdata/Hspa8-Nla3-Pol2_TKD181102142/test/all_bam/bigwig")
  draw_circle(AC12withHBV,"AC12withHBV_o",10)
  draw_circle(AD38A,"AD38A",11)
  draw_circle(AD38noDOX11d,"AD38noDOX11d_o",12)
  draw_circle(AD38noDOX6d,"AD38noDOX6d_o",13)
  draw_circle(AD38withDOX,"AD38withDOX_o",14)
  draw_circle(DE19withDOX,"DE19withDOX_o",15)
 
  circos_plot <- function(data){
    colors   <- rainbow(10, alpha=0.7)
    par(mar=c(0, 0, 0, 0));
    
    plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
    
    circos(R=400, type="chr", cir="hg19", col=colors, print.chr.lab=TRUE, W=4, scale=TRUE);
    circos(R=360, cir="hg19", W=40, mapping=H3K27ac, col.v=4, type="l",   B=TRUE, col=colors[1], lwd=0.1, scale=TRUE);
    circos(R=320, cir="hg19", W=40, mapping=H3K36me3, col.v=4, type="l",   B=TRUE, col=colors[2], lwd=0.1, scale=TRUE);
    circos(R=280, cir="hg19", W=40, mapping=H3K4me1, col.v=4, type="l",   B=TRUE, col=colors[3], lwd=0.1, scale=TRUE);
    circos(R=240, cir="hg19", W=40, mapping=H3K4me2, col.v=4, type="l",   B=TRUE, col=colors[4], lwd=0.1, scale=TRUE);
    circos(R=240, cir="hg19", W=40, mapping=H3K4me3, col.v=4, type="l",   B=TRUE, col=colors[5], lwd=0.1, scale=TRUE);
    circos(R=200, cir="hg19", W=40, mapping=H3K9ac, col.v=4, type="l",   B=TRUE, col=colors[6], lwd=0.1, scale=TRUE);
    circos(R=160, cir="hg19", W=40, mapping=H3K9me3, col.v=4, type="l",   B=TRUE, col=colors[7], lwd=0.1, scale=TRUE);
    circos(R=120, cir="hg19", W=40, mapping=H4K20me1, col.v=4, type="l",   B=TRUE, col=colors[8], lwd=0.1, scale=TRUE);
    circos(R=80, cir="hg19", W=40, mapping=SALP_HepG2, col.v=4, type="l",   B=TRUE, col=colors[9], lwd=0.1, scale=TRUE);
    circos(R=40, cir="hg19", W=40, mapping=data, col.v=4, type="l",   B=TRUE, col=colors[9], lwd=0.1, scale=TRUE);  
}
  
  savefile <- function(data,name){
    pdf(paste0(name,".pdf"),width = 12, height = 12)
    circos_plot(data)
    dev.off()
  }
  savefile(AC12withHBV,"AC12withHBV")
  savefile(AD38A,"AD38A")
  savefile(AD38noDOX11d,"AD38noDOX11d")
  savefile(AD38noDOX6d,"AD38noDOX6d")
  savefile(AD38withDOX,"AD38withDOX")
  savefile(DE19withDOX,"DE19withDOX")
  
  #------ test
  
  options(stringsAsFactors = FALSE);
  library(OmicCircos);
  set.seed(1234);
  
  ## initial values for simulation data 
  seg.num     <- 10;
  ind.num     <- 20;
  seg.po      <- c(20:50);
  link.num    <- 10;
  link.pg.num <- 4;
  ## output simulation data
  sim.out <- sim.circos(seg=seg.num, po=seg.po, ind=ind.num, link=link.num, 
                        link.pg=link.pg.num);
  
  seg.f     <- sim.out$seg.frame;
  seg.v     <- sim.out$seg.mapping;
  link.v    <- sim.out$seg.link
  link.pg.v <- sim.out$seg.link.pg
  seg.num   <- length(unique(seg.f[,1]));
  
  ## select segments
  seg.name <- paste("chr", 1:seg.num, sep="");
  db       <- segAnglePo(seg.f, seg=seg.name);
  
  colors   <- rainbow(seg.num, alpha=0.5);
  
  
  ###################################################
  ### code chunk number 13: OmicCircos4vignette2
  ###################################################
  pdf("test.pdf",width = 12, height = 12)
  par(mar=c(0, 0, 0, 0));
  
  plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
  
  circos(R=400, type="chr", cir=db, col=colors, print.chr.lab=TRUE, W=4, scale=TRUE);
  circos(R=360, cir=db, W=40, mapping=seg.v, col.v=8, type="box",   B=TRUE, col=colors[1], lwd=0.1, scale=TRUE);
  circos(R=320, cir=db, W=40, mapping=seg.v, col.v=8, type="hist",  B=TRUE, col=colors[3], lwd=0.1, scale=TRUE);
  circos(R=280, cir=db, W=40, mapping=seg.v, col.v=8, type="ms",  B=TRUE, col=colors[7], lwd=0.1, scale=TRUE);
  circos(R=240, cir=db, W=40, mapping=seg.v, col.v=3, type="h",  B=FALSE,  col=colors[2], lwd=0.1);
  circos(R=200, cir=db, W=40, mapping=seg.v, col.v=3, type="s", B=TRUE, col=colors, lwd=0.1);
  circos(R=160, cir=db, W=40, mapping=seg.v, col.v=3, type="b", B=FALSE, col=colors, lwd=0.1);
  circos(R=120, cir=db, W=40, mapping=link.v, type="link", lwd=2, col=colors[c(1,7)]);
  circos(R=80, cir=db, W=40, mapping=link.pg.v, type="link.pg", lwd=2, col=sample(colors,link.pg.num));
  
  dev.off()
  
    