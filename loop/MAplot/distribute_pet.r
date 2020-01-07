#--- it's used to draw reads distribution graph
options(stringsAsFactors = FALSE)
# oceanc_wt <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/pol2_oceanc_hichipper/0h.filt.intra.loop_counts.bedpe'
# oceanc_ko <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/pol2_oceanc_hichipper/6h.filt.intra.loop_counts.bedpe'
# 
# hichip_wt <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/pol2_hicip_hichipper/pol2_wt.filt.intra.loop_counts.bedpe'
# hichip_ko <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/pol2_hicip_hichipper/pol2_ko.filt.intra.loop_counts.bedpe'
oceanc_wt <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/MAplot/Hiccups/oceanc_wt.bedpe'
oceanc_ko <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/MAplot/Hiccups/oceanc_ko.bedpe'

hichip_wt <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/MAplot/Hiccups/hichip_wt.bedpe'
hichip_ko <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/MAplot/Hiccups/hichip_ko.bedpe'

pol2_hic_wt <- '/DATA/work/lbyybl/ypjiang/geneloop/pol2_loop/allvalidpairs/PET_distribution/wt_loop.bed'
pol2_hic_ko <- '/DATA/work/lbyybl/ypjiang/geneloop/pol2_loop/allvalidpairs/PET_distribution/ko_loop.bed'


library(dplyr)
library(data.table)
readdata <- function(wt,ko){
  wt_file <- fread(wt, col.names = c('chr1','st1','en1','chr2','st2','en2','score','pets'))
  wt_file$class <- 'wt'
  ko_file <- fread(ko, col.names = c('chr1','st1','en1','chr2','st2','en2','score','pets'))
  ko_file$class <- 'ko'
  data <- bind_rows(wt_file, ko_file)
  return(data)
}
oceanc_data <- readdata(wt=oceanc_wt,ko=oceanc_ko)
hichip_data <- readdata(wt=hichip_wt,ko=hichip_ko)
pol2_hic_data <- readdata(wt=pol2_hic_wt,ko=pol2_hic_ko)
# oceanc_wt_file <- fread(oceanc_wt,col.names = c('chr1','st1','en1','chr2','st2','en2','score','pets'))
# oceanc_wt_file$class <- 'wt'
# oceanc_ko_file <- fread(oceanc_ko,col.names = c('chr1','st1','en1','chr2','st2','en2','score','pets'))
# oceanc_ko_file$class <- 'ko'
# data <- bind_rows(oceanc_wt_file,oceanc_ko_file)
#hichip_wt_file <- fread(hichip_wt,col.names = c('chr1','st1','en1','chr2','st2','en2','score','pets'))
#hichip_ko_file <- fread(hichip_ko,col.names = c('chr1','st1','en1','chr2','st2','en2','score','pets'))

ggplot(oceanc_data,aes(x=pets,y=log(..count..),fill=class)) + geom_bar(position = 'dodge2') +
  labs(y="log(counts of PETs)",x="PETs")+
  scale_x_continuous(limits=c(0,160))+
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
    axis.line.x = element_blank(),#element_line(color="black", size = 0.5),
    axis.line.y = element_blank(), #element_line(color="black", size = 0.5),
    axis.text.x=element_text(angle=0, hjust=0.5,vjust = 0,
                             colour="black", size=12),
    axis.text.y=element_text(#face = "bold",
      colour="black", size=12),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    panel.background=element_blank()
  )
filename="pol2_oceanc_distribution"
#fwrite(ma_input,paste0(filename,".bed"),sep = "\t",row.names = TRUE)
ggsave(paste0(filename,".jpeg"),width =8, height = 4 )
ggsave(paste0(filename,".pdf"),width =8, height = 4 )
ggplot(hichip_data,aes(x=pets,y=log(..count..),fill=class)) + geom_bar(position = 'dodge2') +
  labs(y="log(counts of PETs)",x="PETs")+
  scale_x_continuous(limits=c(0,200))+
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
    axis.line.x = element_blank(),#element_line(color="black", size = 0.5),
    axis.line.y = element_blank(), #element_line(color="black", size = 0.5),
    axis.text.x=element_text(angle=0, hjust=0.5,vjust = 0,
                             colour="black", size=12),
    axis.text.y=element_text(#face = "bold",
      colour="black", size=12),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    panel.background=element_blank()
  )
filename="pol2_hichip_distribution"
#fwrite(ma_input,paste0(filename,".bed"),sep = "\t",row.names = TRUE)
ggsave(paste0(filename,".jpeg"),width =8, height = 4 )
ggsave(paste0(filename,".pdf"),width =8, height = 4 )

ggplot(pol2_hic_data,aes(x=pets,y=log(..count..),fill=class)) + geom_bar(position = 'dodge2') +
  labs(y="log(counts of PETs)",x="PETs")+
  scale_x_continuous(limits=c(0,200))+
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
    axis.line.x = element_blank(),#element_line(color="black", size = 0.5),
    axis.line.y = element_blank(), #element_line(color="black", size = 0.5),
    axis.text.x=element_text(angle=0, hjust=0.5,vjust = 0,
                             colour="black", size=12),
    axis.text.y=element_text(#face = "bold",
      colour="black", size=12),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    panel.background=element_blank()
  )
filename="pol2_hic_distribution"
#fwrite(ma_input,paste0(filename,".bed"),sep = "\t",row.names = TRUE)
ggsave(paste0(filename,".jpeg"),width =8, height = 4 )
ggsave(paste0(filename,".pdf"),width =8, height = 4 )
