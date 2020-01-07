# draw heatmap like plot
library(dplyr)
library(pheatmap)
setwd('/DATA/work/lbyybl/wh/ruvb2/chip-ms')
complex <- read.table('complex_list.txt',header = F, stringsAsFactors = F)
names(complex) <- c('comp','protein_id','gene_id')
protein_ms <- fread('sub_merge.txt',header = T, stringsAsFactors = F)
protein_ms <- protein_ms %>%
  dplyr::select(Accession,unique_peptides,rep1,rep2,rep3,gfp,input,ipvsgfp,ipvsinput)
ms_gfpthan2 <- protein_ms %>%
  filter(unique_peptides >=2 ) %>%
  filter(ipvsgfp >= 2)

merge <- left_join(complex,ms_gfpthan2,by=c('protein_id'='Accession'))
merge[is.na(merge)] <- 0

drawcomplex <- function(name){
  data <- merge %>%
    filter(comp==name) %>%
    arrange(-ipvsgfp)
  # data_matrix <- matrix(data$ipvsgfp,nrow=1,ncol=nrow(data))
  # colnames(data_matrix) <- data$gene_id
  # rownames(data_matrix) <- 'TIP60'
  # pheatmap(data_matrix,cluster_rows = FALSE, cluster_cols = FALSE,legend_breaks=seq(-80,80,2))
  ggplot(data=data, aes(x=gene_id, y=comp, fill=ipvsgfp))+geom_raster()+
    #ggtitle(paste0("Averaged Interaction within TAD\n(KPNA2 ", arg, " merge)"))+
    scale_fill_gradientn(colors=c('white','red1'),values=c(0,0.45,0.5,1),limit=c(0,30),
    )+scale_x_discrete(expand=c(0, 0),limits = data$gene_id)+
    scale_y_discrete(expand=c(0, 0))+
    theme(panel.border = element_rect(colour="black",  fill=NA, size=1),
          axis.title.y = element_blank(),
          axis.title.x = element_blank())                                       
  ggsave(paste0(name,'.pdf'),width=6,height = 1)
  
}

for (i in unique(merge$comp)){
  drawcomplex(i)
}
