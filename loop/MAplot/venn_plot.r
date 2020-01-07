#--- this is plot venn to see the overlap between ocean hihcip
#--- it to the the loop identified by MAplot.r
setwd('/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/MAplot')
library(data.table)
library(VennDiagram)
library(dplyr)
hichip <- "/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/MAplot/pol2_hichip.bed"
oceanc <- "/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/MAplot/pol2_oceanc.bed"

#hichip_file <- fread(hichip)
#oceanc_file <- fread(oceanc)
got_allupdown<- function(file){
  all <- fread(file)
  all_up <- all %>%
    filter(change=="up")
  all_down <- all %>%
    filter(change=="down")
  return(list(all,all_up,all_down))
}
hichip_li <- got_allupdown(hichip)
oceanc_li <- got_allupdown(oceanc)

venn.diagram(list(HiChip=hichip_li[[1]]$V1,Ocean_C=oceanc_li[[1]]$V1), fill=c("red","green"), alpha=c(0.5,0.5), cex=2, cat.fontface=4, fontfamily=2,
             filename="Vennall.tiff")
venn.diagram(list(HiChip=hichip_li[[2]]$V1,Ocean_C=oceanc_li[[2]]$V1), fill=c("red","green"), alpha=c(0.5,0.5), cex=2, cat.fontface=4, fontfamily=2,
             filename="Vennup.tiff")
venn.diagram(list(HiChip=hichip_li[[3]]$V1,Ocean_C=oceanc_li[[3]]$V1), fill=c("red","green"), alpha=c(0.5,0.5), cex=2, cat.fontface=4, fontfamily=2,
             filename="Venndown.tiff")

#--- 做venn图，使用减法判断上升，下降，并可以做显著性检验；
#--- 先使用相减相加之后的平均值作为cut off; 然后找出下降，和上升的区域
#--- 做venn图；如果检验可以检验对应loop ko和wt的均值是否相等；
#--- 我觉得不该做Z-score normalize
setwd("/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/MAplot/")
#setwd("/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/MAplot/")
# 读取表达量的表格
options(stringsAsFactors = FALSE)
#library(edgeR)
library(ggplot2)
library(data.table)
library(VennDiagram)
library(dplyr)
#--- read 4 file and merge to one file
hichip_wt1 <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/MAplot/pol2_h3k27_hichip_WT_1_loop.bed'
hichip_wt2 <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/MAplot/pol2_h3k27_hichip_WT_2_loop.bed'
hichip_ko1 <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/MAplot/pol2_h3k27_hichip_KO_1_loop.bed'
hichip_ko2 <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/MAplot/pol2_h3k27_hichip_KO_2_loop.bed'
oceanc_wt1 <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/MAplot/pol2_oceanc_WT_1_loop.bed'
oceanc_wt2 <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/MAplot/pol2_oceanc_WT_2_loop.bed'
oceanc_ko1 <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/MAplot/pol2_oceanc_KO_1_loop.bed'
oceanc_ko2 <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/MAplot/pol2_oceanc_KO_2_loop.bed'
readfile <- function(file,name){
  data <- read.table(file,header = FALSE)
  colnames(data) <- c("loc",name)
  return(data)
}
up_merge_file <- function(wt1,wt2,ko1,ko2){
  wt1_file <- readfile(wt1,"wt1")
  wt2_file <- readfile(wt2,"wt2")
  ko1_file <- readfile(ko1,"ko1")
  ko2_file <- readfile(ko2,"ko2")
  
  wt_merge <- merge(wt1_file,wt2_file,by="loc", sort=TRUE,all=TRUE)
  ko_merge <- merge(ko1_file,ko2_file,by="loc", sort=TRUE,all=TRUE)
  all_merge <- unique(merge(wt_merge,ko_merge,by="loc",all=TRUE))
  all_merge[is.na(all_merge)] <- 0
  all_merge <- all_merge %>%
    mutate(change=ko1+ko2-wt1-wt2)
  up_mean <- all_merge %>% 
    filter(change >0) %>%
    summarise(mean(change))
  up_data <- all_merge %>% 
    filter(change > up_mean[1,1]) 
  return(up_data)
}

down_merge_file <- function(wt1,wt2,ko1,ko2){
  wt1_file <- readfile(wt1,"wt1")
  wt2_file <- readfile(wt2,"wt2")
  ko1_file <- readfile(ko1,"ko1")
  ko2_file <- readfile(ko2,"ko2")
  
  wt_merge <- merge(wt1_file,wt2_file,by="loc", sort=TRUE,all=TRUE)
  ko_merge <- merge(ko1_file,ko2_file,by="loc", sort=TRUE,all=TRUE)
  all_merge <- unique(merge(wt_merge,ko_merge,by="loc",all=TRUE))
  all_merge[is.na(all_merge)] <- 0
  all_merge <- all_merge %>%
    mutate(change=ko1+ko2-wt1-wt2)
  down_mean <- all_merge %>% 
    filter(change <0) %>%
    summarise(mean(change))
  down_data <- all_merge %>% 
    filter(change < down_mean[1,1]) 
  return(down_data)
}
all_merge_file <- function(wt1,wt2,ko1,ko2){
  wt1_file <- readfile(wt1,"wt1")
  wt2_file <- readfile(wt2,"wt2")
  ko1_file <- readfile(ko1,"ko1")
  ko2_file <- readfile(ko2,"ko2")
  
  wt_merge <- merge(wt1_file,wt2_file,by="loc", sort=TRUE,all=TRUE)
  ko_merge <- merge(ko1_file,ko2_file,by="loc", sort=TRUE,all=TRUE)
  all_merge <- unique(merge(wt_merge,ko_merge,by="loc",all=TRUE))
  all_merge[is.na(all_merge)] <- 0
  all_merge <- all_merge %>%
    mutate(change=ko1+ko2-wt1-wt2)
  return(all_merge)
}
up_zscore_file <- function(wt1,wt2,ko1,ko2){
  wt1_file <- readfile(wt1,"wt1")
  wt2_file <- readfile(wt2,"wt2")
  ko1_file <- readfile(ko1,"ko1")
  ko2_file <- readfile(ko2,"ko2")
  
  wt_merge <- merge(wt1_file,wt2_file,by="loc", sort=TRUE,all=TRUE)
  ko_merge <- merge(ko1_file,ko2_file,by="loc", sort=TRUE,all=TRUE)
  all_merge <- unique(merge(wt_merge,ko_merge,by="loc",all=TRUE)) 
  all_merge[is.na(all_merge)] <- 0
  all_merge <- all_merge %>%
    mutate(wt1=(wt1-mean(wt1))/var(wt1),wt2=(wt2-mean(wt2))/var(wt2),ko1=(ko1-mean(ko1))/var(ko1),ko2=(ko2-mean(ko2))/var(ko2))
  all_merge <- all_merge %>%
    mutate(change=ko1+ko2-wt1-wt2)
  up_mean <- all_merge %>% 
    filter(change >0) %>%
    summarise(mean(change))
  up_data <- all_merge %>% 
    filter(change > up_mean[1,1]) 
  return(up_data)
}
down_zscore_file <- function(wt1,wt2,ko1,ko2){
  wt1_file <- readfile(wt1,"wt1")
  wt2_file <- readfile(wt2,"wt2")
  ko1_file <- readfile(ko1,"ko1")
  ko2_file <- readfile(ko2,"ko2")
  
  wt_merge <- merge(wt1_file,wt2_file,by="loc", sort=TRUE,all=TRUE)
  ko_merge <- merge(ko1_file,ko2_file,by="loc", sort=TRUE,all=TRUE)
  all_merge <- unique(merge(wt_merge,ko_merge,by="loc",all=TRUE)) 
  all_merge[is.na(all_merge)] <- 0
  all_merge <- all_merge %>%
    mutate(wt1=(wt1-mean(wt1))/var(wt1),wt2=(wt2-mean(wt2))/var(wt2),ko1=(ko1-mean(ko1))/var(ko1),ko2=(ko2-mean(ko2))/var(ko2))
  all_merge <- all_merge %>%
    mutate(change=ko1+ko2-wt1-wt2)
  down_mean <- all_merge %>% 
    filter(change <0) %>%
    summarise(mean(change))
  down_data <- all_merge %>% 
    filter(change < down_mean[1,1]) 
  return(down_data)
}
hichip_all <- all_merge_file(wt1 = hichip_wt1, wt2 = hichip_wt2, ko1 = hichip_ko1, ko2 = hichip_ko2)
hichip_up <- up_merge_file(wt1 = hichip_wt1, wt2 = hichip_wt2, ko1 = hichip_ko1, ko2 = hichip_ko2)
hichip_down <- down_merge_file(wt1 = hichip_wt1, wt2 = hichip_wt2, ko1 = hichip_ko1, ko2 = hichip_ko2)
hichip_zscore_up <- up_zscore_file(wt1 = hichip_wt1, wt2 = hichip_wt2, ko1 = hichip_ko1, ko2 = hichip_ko2)
hichip_zscore_down <- down_zscore_file(wt1 = hichip_wt1, wt2 = hichip_wt2, ko1 = hichip_ko1, ko2 = hichip_ko2)

oceanc_all <- all_merge_file(wt1 = oceanc_wt1, wt2 = oceanc_wt2, ko1 = oceanc_ko1, ko2 = oceanc_ko2)
oceanc_up <- up_merge_file(wt1 = oceanc_wt1, wt2 = oceanc_wt2, ko1 = oceanc_ko1, ko2 = oceanc_ko2)
oceanc_down <- down_merge_file(wt1 = oceanc_wt1, wt2 = oceanc_wt2, ko1 = oceanc_ko1, ko2 = oceanc_ko2)
oceanc_zscore_up <- up_zscore_file(wt1 = oceanc_wt1, wt2 = oceanc_wt2, ko1 = oceanc_ko1, ko2 = oceanc_ko2)
oceanc_zscore_down <- down_zscore_file(wt1 = oceanc_wt1, wt2 = oceanc_wt2, ko1 = oceanc_ko1, ko2 = oceanc_ko2)

#---venn plot
venn.diagram(list(HiChip=hichip_all$loc,Ocean_C=oceanc_all$loc), fill=c("red","green"), alpha=c(0.5,0.5), cex=2, cat.fontface=4, fontfamily=2,
             filename="Vennall.tiff")
venn.diagram(list(HiChip=hichip_up$loc,Ocean_C=oceanc_up$loc), fill=c("red","green"), alpha=c(0.5,0.5), cex=2, cat.fontface=4, fontfamily=2,
             filename="Vennup.tiff")
venn.diagram(list(HiChip=hichip_down$loc,Ocean_C=oceanc_down$loc), fill=c("red","green"), alpha=c(0.5,0.5), cex=2, cat.fontface=4, fontfamily=2,
             filename="Venndown.tiff")

venn.diagram(list(HiChip=hichip_zscore_up$loc,Ocean_C=oceanc_zscore_up$loc), fill=c("red","green"), alpha=c(0.5,0.5), cex=2, cat.fontface=4, fontfamily=2,
             filename="Venn_zscore_up.tiff")
venn.diagram(list(HiChip=hichip_zscore_down$loc,Ocean_C=oceanc_zscore_down$loc), fill=c("red","green"), alpha=c(0.5,0.5), cex=2, cat.fontface=4, fontfamily=2,
             filename="Venn_zscore_down.tiff")

#---- Venn plot used FDR
#---- 1. 取交集；2. 用交集看上调显著的，和下调显著的；
  #--- 1. 取交集；
  setwd("/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/MAplot/")
  #setwd("/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/MAplot/")
  # 读取表达量的表格
  options(stringsAsFactors = FALSE)
  #library(edgeR)
  library(ggplot2)
  library(data.table)
  library(VennDiagram)
  library(dplyr)
  library(purrr)
  #--- read 4 file and merge to one file
  hichip_wt1 <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/MAplot/pol2_h3k27_hichip_WT_1_loop.bed'
  hichip_wt2 <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/MAplot/pol2_h3k27_hichip_WT_2_loop.bed'
  hichip_ko1 <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/MAplot/pol2_h3k27_hichip_KO_1_loop.bed'
  hichip_ko2 <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/MAplot/pol2_h3k27_hichip_KO_2_loop.bed'
  oceanc_wt1 <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/MAplot/pol2_oceanc_WT_1_loop.bed'
  oceanc_wt2 <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/MAplot/pol2_oceanc_WT_2_loop.bed'
  oceanc_ko1 <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/MAplot/pol2_oceanc_KO_1_loop.bed'
  oceanc_ko2 <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/MAplot/pol2_oceanc_KO_2_loop.bed'
  readfile <- function(file,name){
    data <- read.table(file,header = FALSE)
    colnames(data) <- c("loc",name)
    return(data)
  }
  
  all_merge_file <- function(wt1,wt2,ko1,ko2){
    wt1_file <- readfile(wt1,"wt1")
    wt2_file <- readfile(wt2,"wt2")
    ko1_file <- readfile(ko1,"ko1")
    ko2_file <- readfile(ko2,"ko2")
    
    wt_merge <- merge(wt1_file,wt2_file,by="loc", sort=TRUE,all=TRUE)
    ko_merge <- merge(ko1_file,ko2_file,by="loc", sort=TRUE,all=TRUE)
    all_merge <- unique(merge(wt_merge,ko_merge,by="loc",all=TRUE))
    all_merge[is.na(all_merge)] <- 0
    all_merge <- all_merge %>%
      mutate(change=ko1+ko2-wt1-wt2)
    ttest <- function(i){
      a <- as.vector(t(all_merge[i,2:3]))
      b <- as.vector(t(all_merge[i,4:5]))
      c <- t.test(a,b,paired = T)
      return(c$p.value)
    }
    all_merge$pvalue <- 1:nrow(all_merge) %>%
      map_dbl(ttest)
    all_merge$FDR <- p.adjust(all_merge$pvalue,method="fdr",n=length(all_merge$pvalue))
    return(all_merge)
  }
  up_sig_data <- function(all_merge){
    all_merge_sig_up <- all_merge %>%
      filter(change > 0) %>%
      filter(FDR < 0.05)
  } 
  down_sig_data <- function(all_merge){
    all_merge_sig_down <- all_merge %>%
      filter(change < 0) %>%
      filter(FDR < 0.05)
  } 
  hichip_all <- all_merge_file(wt1 = hichip_wt1, wt2 = hichip_wt2, ko1 = hichip_ko1, ko2 = hichip_ko2)
  
  oceanc_all <- all_merge_file(wt1 = oceanc_wt1, wt2 = oceanc_wt2, ko1 = oceanc_ko1, ko2 = oceanc_ko2)
  
  #---venn plot
  venn.diagram(list(HiChip=hichip_all$loc,Ocean_C=oceanc_all$loc), fill=c("red","green"), alpha=c(0.5,0.5), cex=2, cat.fontface=4, fontfamily=2,
               filename="Vennall.tiff")
  venn.diagram(list(HiChip=hichip_up$loc,Ocean_C=oceanc_up$loc), fill=c("red","green"), alpha=c(0.5,0.5), cex=2, cat.fontface=4, fontfamily=2,
               filename="Vennup.tiff")
  venn.diagram(list(HiChip=hichip_down$loc,Ocean_C=oceanc_down$loc), fill=c("red","green"), alpha=c(0.5,0.5), cex=2, cat.fontface=4, fontfamily=2,
               filename="Venndown.tiff")
  
  venn.diagram(list(HiChip=hichip_zscore_up$loc,Ocean_C=oceanc_zscore_up$loc), fill=c("red","green"), alpha=c(0.5,0.5), cex=2, cat.fontface=4, fontfamily=2,
               filename="Venn_zscore_up.tiff")
  venn.diagram(list(HiChip=hichip_zscore_down$loc,Ocean_C=oceanc_zscore_down$loc), fill=c("red","green"), alpha=c(0.5,0.5), cex=2, cat.fontface=4, fontfamily=2,
               filename="Venn_zscore_down.tiff")
