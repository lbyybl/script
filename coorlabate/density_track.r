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
  if (chr=="chr2"){
    chrom            = "chr2"
    chromstart       = 1
    chromend         = 243199373
  }else if (chr=="chr8"){
    chrom            = "chr8"
    chromstart       = 1
    chromend         = 146364022
  }else if (chr=="chr21"){
    chrom            = "chr21"
    chromstart       = 1
    chromend         = 48129895
  }else if (chr=="chr22"){
    chrom            = "chr22"
    chromstart       = 1
    chromend         = 51304566
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
  ma_va <- 5000 #findsubmax(file) # ma_va <- findsubmax(AD38withDOX_file)
  ma_round <- round(ma_va/100)*100
  lables <- seq(0,ma_round,ma_round/5)
  drawwhochr(file,"chr2",name,ma_round,lables)
  drawwhochr(file,"chr8",name,ma_round,lables)
  drawwhochr(file,"chr21",name,ma_round,lables)
  drawwhochr(file,"chr22",name,ma_round,lables)
}
drawlocal4 <- function(file,name){
  ma_va <- 5000 #findsubmax(file)
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

