# it's used to draw boxplot and network for protein complex
library(dplyr)
library(data.table)
library(readr)
library(ggplot2)
setwd('/DATA/work/lbyybl/wh/ruvb2/chip-ms/complex')
comp <- fread('complex_all.csv',header = F,stringsAsFactors = F)
comp <- comp[,1:4]
colnames(comp) <- c('compx','name','exit','id')
info <- fread('all_info.csv',header = T,stringsAsFactors = F)
info <- info %>%
  select(Accession,ipvsinput)
#---for boxplot
comp_det <- comp %>%
  filter(exit == 1)
comp_m <- merge(comp_det,info,by.x='id',by.y='Accession')
ggplot(comp_m,aes(x=compx,y=ipvsinput)) + geom_boxplot()+
  geom_jitter(aes(color=compx),width = 0.3)+
  scale_x_discrete(limits=c('SWR','TIP','ASTRA','pre-snoRNA','prefoldin',
                            'INO','R2TP','TFTC','LUBAC','POL_I','POL_III','replication_fork',
                            'GISN','ribose','endosome','POL_II'))+
  labs(x = NULL, y = "Protein enrichment")+
  theme(
    #title = element_blank() ,
    legend.key.width = unit(2,"line"),
    legend.key.height = unit(0.5,"cm"),
    legend.text = element_text(
      size = 7,
      hjust = 1.2,
      face = "italic",
      colour = "black"
    ),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    plot.title=element_text(colour="black", size=10),
    axis.text.x=element_text(angle = -30,hjust = -0.05,
      colour="black", size=7),
    axis.text.y=element_text(#face = "bold",
      colour="black", size=12),
    panel.background=element_blank())
  
comp_m %>%
  group_by(compx) %>%
  summarise(x=mean(ipvsinput)) %>%
  arrange(x)

ggsave('complex_enrichment.pdf',width = 6,height = 2)
#-----------------------------------------------------------------
# draw network graph
comp_all <- merge(comp,info,by.x='id',by.y='Accession',all.x=T)
comp_all[is.na(comp_all)] <- 0
fwrite(comp_all,'complex_all.csv')

library(igraph)
test <- fread('complex_all.csv')

matrix <- matrix(0,nrow = length(unique(c(test$compx,test$name))),
                 ncol = length(unique(c(test$compx,test$name))))
colnames(matrix) <- unique(c(test$compx,test$name))
rownames(matrix) <- unique(c(test$compx,test$name))
for (i in 1:nrow(test)){
  matrix[test$name[i],test$compx[i]] <- 1
  matrix[test$compx[i],test$name[i]] <- 1
}
gr <- graph.adjacency(matrix,mode='undirected')
palf <- colorRampPalette(c("white", 'palegreen','chartreuse', "yellow",'gold',"red"))
#V(g)$color <- palf[V(test)$ipvsi]
#range1.100 <- function(x){1 + 99*(x-min(x))/(max(x)-min(x))}
n <- 270
colr <- palf(n);
V(gr)$ipvsi <- as.numeric(test$ipvsinput[match(V(gr)$name,test$name)])
V(gr)$ipvsi[is.na(V(gr)$ipvsi)] <- 0
V(gr)$color <- colr[log(as.numeric(V(gr)$ipvsi)+1)*100+1]
pdf('all_complex.pdf')
#par(mfrow=c(1,2))
plot(gr,
     #layout=locs, 
     vertex.size=5, 
     vertex.label.cex=0.2,
     main="Original",
     vertex.color=V(gr)$color)
dev.off()
c_clust <- cluster_fast_greedy(gr)
pdf('cluster.pdf')
plot(c_clust,gr,
     vertex.size=5, 
     vertex.label.cex=0.2,
     main="Original",
     vertex.color=V(gr)$color)
dev.off()
pdf('legent.pdf')
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
plot(1:n, 1:n, pch = 19, cex=2, col = colr)

legend_image <- as.raster(matrix(colr, ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, 0, 0, 1,1)
dev.off()

#-------------------------
# 各个复合体分开画图
draw_compx <- function(copx){
  #copx='TIP'
  output_name=paste0(copx,'.pdf')
  data <- test %>%
    filter(compx==copx)
  
  data_matrix <- matrix(0,nrow = length(unique(c(data$compx,data$name))),
                        ncol = length(unique(c(data$compx,data$name))))
  colnames(data_matrix) <- unique(c(data$compx,data$name))
  rownames(data_matrix) <- unique(c(data$compx,data$name))
  for (i in 1:nrow(data)){
    data_matrix[data$name[i],data$compx[i]] <- 1
    data_matrix[data$compx[i],data$name[i]] <- 1
  }
  data_gr <- graph.adjacency(data_matrix,mode='undirected')
  #range1.100 <- function(x){1 + 99*(x-round(min(test$ipvsinput)))/(round(max(test$ipvsinput))-round(min(test$ipvsinput)))}
  colr <- palf(n);
  V(data_gr)$ipvsi <- as.numeric(test$ipvsinput[match(V(data_gr)$name,test$name)])
  V(data_gr)$ipvsi[is.na(V(data_gr)$ipvsi)] <- 0
  V(data_gr)$color <- colr[log(as.numeric(V(data_gr)$ipvsi)+1)*100+1]
  pdf(output_name)
  plot(data_gr,
       #layout=locs, 
       vertex.size=30, 
       vertex.label.cex=0.5,
       main="Original",
       vertex.color=V(data_gr)$color)
  dev.off()
}

for (i in unique(test$compx)){
  draw_compx(i)
}
