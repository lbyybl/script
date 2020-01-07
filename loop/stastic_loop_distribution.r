# this is to stastic the distribution of loop pet for hichip and ocean-c
#---ocean-c
library(data.table)
library(ggplot2)
setwd("/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/pol2_oceanc_hichipper")
wt <- '0h.filt.intra.loop_counts.bedpe'
ko <- '6h.filt.intra.loop_counts.bedpe'
wt_oceanc <- fread(wt)
wt_ocp <- data.frame("num"=wt_oceanc$V8,
                     "class"="wt")
ko_oceanc <- fread(ko)
ko_ocp <- data.frame("num"=ko_oceanc$V8,
                     "class"="ko")
ocp <- rbind(wt_ocp,ko_ocp)
ggplot(ocp,aes(x=num,y=log2(..count..),fill=class))+ geom_bar(position = "dodge2") 

setwd("/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/pol2_hicip_hichipper")
wt <- 'pol2_wt.filt.intra.loop_counts.bedpe'
ko <- 'pol2_ko.filt.intra.loop_counts.bedpe'
wt_oceanc <- fread(wt)
wt_ocp <- data.frame("num"=wt_oceanc$V8,
                     "class"="wt")
ko_oceanc <- fread(ko)
ko_ocp <- data.frame("num"=ko_oceanc$V8,
                     "class"="ko")
ocp <- rbind(wt_ocp,ko_ocp)
ggplot(ocp,aes(x=num,y=log2(..count..),fill=class))+ geom_bar(position = "dodge2") 
# ggplot(wt_ocp,aes(x=wt_n,y=log2(..count..))) + geom_bar(color="green") + 
#   geom_bar(data=ko_ocp,aes(x=ko_n,y=log2(..count..)),color="red")