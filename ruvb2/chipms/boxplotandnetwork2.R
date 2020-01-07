setwd('/DATA/work/lbyybl/wh/ruvb2/chip-ms/')
library(igraph)
test <- fread('complex_list.txt',col.names = c('compx','id','name'))
test <- test %>%
  filter(compx== 'SRCAP' | compx== 'TIP60')

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
plot(gr,
     #layout=locs, 
     vertex.size=5, 
     vertex.label.cex=0.2,
     main="Original",
     vertex.color=V(gr)$color)
dev.off()
