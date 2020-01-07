# Get the library.
library(plotrix)

# Create data for the graph.
x <-  c(21, 62, 10,53)
lbl <-  c("London","New York","Singapore","Mumbai")

# Give the chart file a name.
png(file = "3d_pie_chart.jpg")

# Plot the chart.
pie3D(x,labels = lbl,explode = 0.1, main = "Pie Chart of Countries ")

# Save the file.
dev.off()

library(data.table)
library(R.utils)
library(dplyr)
setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/Ruvbl1/sample20190509/graph/subtract_heatmap/cluster')
matrix <- fread('distribution_complex.gz',skip = 1,header = F,
                stringsAsFactors = F)
# colname <- c('chr','st','en','loc','p1','p2',
#              paste0('Ino',1:1000),paste0('R1',1:1000),paste0('R2',1:1000),
#              paste0('Tip',1:1000),paste0('Znh',1:1000))
# colnames(matrix) <- colname
matrix_x <- data.frame('ino'=rowSums(matrix[,7:1006]),
                       'R1'=rowSums(matrix[,1007:2006]),
                       'R2'=rowSums(matrix[,2007:3006]),
                       'Tip'=rowSums(matrix[,3007:4006]),
                       'znh'=rowSums(matrix[,4007:5006]))
rownames(matrix_x) <- matrix$V4

matrix2<-as.matrix(matrix_x)
library(pheatmap)
library(stats)
#matrix3 <- matrix2 
matrix3 <- matrix2[1:300,]
a <- pheatmap(matrix2,cluster_rows = T,cluster_cols = F,scale='row',clustering_method = 'ward')

pheatmap(matrix2)

#---------------------------------------
setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/Ruvbl1/sample20190509/graph/subtract_heatmap/Ruvb1_2_overlap_motif')
known_motif <- fread('knownResults.txt')
known_motif$tar <- as.numeric(gsub('%',"",known_motif$`% of Target Sequences with Motif`))
known_motif$bag <- as.numeric(gsub('%','',known_motif$`% of Background Sequences with Motif`))
known_motif$enrich <- known_motif$tar/known_motif$bag
#---------------------------------------
setwd('/DATA/work/lbyybl/wh/ruvb2/chip-ms')
chipms <- fread('WH_20190620_R1_R2_ChIP-MS.csv',header = T,stringsAsFactors = F)
unipro <- fread('uniprot/uniprot_tab.tsv',col.names = c('id','gene_name','loc','go'),stringsAsFactors = F)
merge <- left_join(chipms,unipro,by=c("Accession"="id"))
fwrite(merge,'chipms_loc_go.txt',sep = '\t')

#names(merge)[21:25] <- c('rep1','rep2','rep3','gfp','input')
merge[is.na(merge)] <- 0
merge <- merge %>%
  mutate(ipvsgfp=(rep1+rep2+rep3+3)/((gfp+1)*3),ipvsinput=(rep1+rep2+rep3+3)/((input+1)*3)) %>%
  mutate(enrich=ipvsgfp+ipvsinput)
# merge <- merge %>%
#   filter(rep1>0,rep2>0,rep3>0)
# merge$enrich <- if_else(rep1>0 & rep2>0 & rep3>0,enrich,0)
# for (i in 1:nrow(merge)){
#   if (merge$rep1[i]>0 & merge$rep2[i]>0 & merge$rep2[i]>0){
#     merge$enrich[i] <- merge$enrich[i]
#   }else{
#     merge$enrich[i] <- 0
#     merge$ipvsgfp[i] <- 0
#     merge$ipvsinput[i] <- 0
#   }
# }
# merge <- merge %>% 
#   filter(rep1>0,rep2>0,rep3>0)
# merge <- merge %>%
#   arrange(-enrich,-(rep1+rep2+rep3))
merge <- merge %>%
  mutate(INO=' ',TIP=' ',SRCAP=' ',R2TP=' ','haha'=' ')
a<-b<-c<-d<-1
srcap_l<-tip_l<-ino_l<-r2tp_l<-list()
# srcap<-c("SRCAP","P60122","Q9WTM5","Q62481","Q9Z2N8","Q9CR11","P60710","Q9JI44","Q9D864")
# tip <- c("Q8CHK4","P60122","Q9WTM5","Q62481","Q9Z2N8","Q9CR11","P60710","Q9JI44","Q80YV3","P60762","Q2VPQ9","Q9DAT2","Q8C9X6","Q8VEK6","Q8R3B7","Q8CHI8")
# ino<-c("Q6ZPV2","P60122","Q9WTM5","Q80US4","Q8R2S9","Q9Z2N8","Q99PT3","Q8BHA0","Q3U1J1","Q6PIJ4","Q99L90")
# r2tp<-c("Q9CQJ2","P60122","Q9WTM5","Q9D706")
srcap<-c("A0A087WQ44","P60122","Q9WTM5","Q62481","Q9Z2N8","Q9CR11","P60710","Q9JI44","Q9D864","Q8R331")
tip <- c("Q8CHK4","P60122","Q9WTM5","Q62481","Q9Z2N8","Q9CR11","P60710","Q9JI44","Q80YV3","P60762","Q2VPQ9","Q9DAT2","Q8C9X6","Q8VEK6","Q8R3B7","Q8CHI8")
ino<-c("Q6ZPV2","P60122","Q9WTM5","Q9Z2N8","Q80US4","Q8R2S9","Q8BHA0","Q99PT3","Q66JY2","A0A0U1RP99","Q99L90","Q6PIJ4","Q3U1J1","Q9WUP7","Q00899")
r2tp<-c("Q9WTM5","P60122","Q9CQJ2","Q9D706")
for (i in 1:nrow(merge)){
  if (merge$Accession[i] %in% srcap){
    merge$SRCAP[i] <- 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
    srcap_l[a]<-i
    a<-a+1
  }
  if (merge$Accession[i] %in% tip){
    merge$TIP[i] <- 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
    tip_l[b]<-i
    b<-b+1
  }
  if (merge$Accession[i] %in% ino){
    merge$INO[i] <- 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
    ino_l[c]<-i
    c<-c+1
  }
  if (merge$Accession[i] %in% r2tp){
    merge$R2TP[i] <- 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
    r2tp_l[d]<-i
    d<-d+1
  }
}

fwrite(merge,'chipms_loc_go.txt',sep = '\t')
names(merge)[11] <- 'unique_peptides'
merge_sub <- merge %>%
  filter(unique_peptides >= 2)
merge_sub <- merge_sub %>%
  arrange(-ipvsgfp,-ipvsinput)
fwrite(merge_sub,'sub_merge.txt',sep='\t')

#--------------------------------------------------------------------
# find cutoff
data <- fread('sub_merge.txt',stringsAsFactors = F,header = T)
data <- data %>%
  filter((rep1>0)+(rep2>0)+(rep3>0)>=2)
names(data)
getz <- function(data,x,y){
  data <- data %>%
    filter(data$ipvsgfp>x,data$ipvsinput>y)
  z <- nrow(data)
  return(z)
}
getcmp <- function(data,x,y){
  data <- data %>%
    filter(data$ipvsgfp>x,data$ipvsinput>y)
  data <- data %>%
    filter(data$INO!='' | data$TIP !='' | data$SRCAP!='' | data$R2TP!='')
  z<-length(unique(data$Accession))
  return(z)
}
da_se <- data.frame('x_s'=seq(0,1,0.05),
                    'y_s'=seq(-1,1,0.1))
df <- data.frame('x'=rep(0,nrow(da_se)*nrow(da_se)),
                 'y'=rep(0,nrow(da_se)*nrow(da_se)),
                 'z'=rep(0,nrow(da_se)*nrow(da_se)))
for (i in 1:nrow(da_se)){
  for (j in 1:nrow(da_se)){
    df$x[(i-1)*nrow(da_se)+j] <- (da_se$x_s[i])
    df$y[(i-1)*nrow(da_se)+j] <- (da_se$y_s[j])
    df$z[(i-1)*nrow(da_se)+j] <- getz(data,exp(da_se$x_s[i]),exp(da_se$y_s[j]))
  }
}
df2=df
df2$z <- as.character(df2$z)
#image(df$x,df$y,df$z)
df$x <- as.character(df$x)
df$y <- as.character(df$y)
scaleFUN <- function(x) sprintf("%.2f", exp(as.numeric(x)))
ggplot(df2,aes(x,y,fill=z))+geom_raster()+geom_text(label=df$z)+
  scale_x_discrete(labels = scaleFUN,limits=as.character(sort(as.numeric(df$x))))+
  scale_y_discrete(labels = scaleFUN,limits=as.character(sort(as.numeric(df$y))))
data2 <- data %>%
  filter(ipvsgfp>=1.49,ipvsinput>=0.5)
fwrite(data2,'filter_gfp_input.txt',sep = '\t')

#---------------------------------------------------------------------
# select the protein exit in JX2015 and find the percentage of Ruvbl2 related complex protein in these protein with this chip-ms
setwd('/DATA/work/lbyybl/wh/ruvb2/chip-ms')
select <- fread('sub_merge.txt',header = T,stringsAsFactors = F)
names(select)[14] <- 'MW'
select <- select %>%
  dplyr::select(gene_name,Accession,rep1,rep2,rep3,MW)
JXprotein <- fread('JX2015.txt',header = F,stringsAsFactors = F)
names(JXprotein) <- 'gene_name'
merge <- merge(select,JXprotein,by='gene_name')
nrow(merge)
merge[is.na(merge)] <- 0
merge <- merge %>%
  dplyr::mutate(qantity=(rep1+rep2+rep3)/MW)
merge$stat <- 'notin'
head(merge)
srcap<-c("A0A087WQ44","P60122","Q9WTM5","Q62481","Q9Z2N8","Q9CR11","P60710","Q9JI44","Q9D864","Q8R331")
tip <- c("Q8CHK4","P60122","Q9WTM5","Q62481","Q9Z2N8","Q9CR11","P60710","Q9JI44","Q80YV3","P60762","Q2VPQ9","Q9DAT2","Q8C9X6","Q8VEK6","Q8R3B7","Q8CHI8")
ino<-c("Q6ZPV2","P60122","Q9WTM5","Q9Z2N8","Q80US4","Q8R2S9","Q8BHA0","Q99PT3","Q66JY2","A0A0U1RP99","Q99L90","Q6PIJ4","Q3U1J1","Q9WUP7","Q00899")
Ruvbl_chrom <- c(srcap,tip,ino)
for (i in 1:nrow(merge)){
  if (merge$Accession[i] %in% Ruvbl_chrom){
    merge$stat[i] <- 'in'
  }
}
merge %>%
  group_by(stat) %>%
  summarise(n=n())
explain_per<-sum(merge[merge$stat=="in",]$qantity)/sum(merge$qantity)
sum(merge[merge$stat=="notin",]$qantity)
fwrite(merge,'Ruvbl_chrom_detect.csv')
#---------------------------------------------------------------------
setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/Ruvbl1/sample20190509/graph/subtract_heatmap/distribution_heatmap')
tad <- fread('pol2_wt_25kb_domain.merged.bed')
mean(tad$V3-tad$V2)
loop <- fread('merged_loops.bedpe',skip = 2)
mean(loop$V5-loop$V2)

#---------------------------------------------------------------------
setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/Ruvbl1/sample20190509/graph/subtract_heatmap/distribution_heatmap')
library(R.utils)
matrix <- fread('distribution_heatmap.gz',skip=1,header = F)
matrix2 <- matrix[,7:4006]
plot(1:4,matrix_x[2,])
matrix_x <- data.frame('CTCF'=rowSums(matrix[,7:1006]),
                       'R1'=rowSums(matrix[,1007:2006]),
                       'R2'=rowSums(matrix[,2007:3006]),
                       'Cohesin'=rowSums(matrix[,3007:4006]))
matrix_fiter <- matrix_x %>%
  filter(CTCF+R1+R2+Cohesin > 0)
matrix_x$class <- 'c'
matrix_x %>%
  #group_by(class)
  summarise_all(mean)
#---------------------------------------------------------------------
library(networkD3)
library(data.table)
setwd('/DATA/work/lbyybl/wh/ruvb2/chip-ms/go/')
net <- fread('complex_det.csv',header = F,stringsAsFactors = F)
simpleNetwork(net,height = 1,width = 1,zoom = T)

# Load data
data(MisLinks)
data(MisNodes)

# Plot
forceNetwork(Links = MisLinks, Nodes = MisNodes,
             Source = "source", Target = "target",
             Value = "value", NodeID = "name",
             Group = "group", opacity = 0.8,fontSize = 25)
#--------------------------------------------------------------
# igraph test
library(igraph)
# Erdos
par(mfrow=c(1,2))
g <- erdos.renyi.game(100, 1/100)
V(g)$size<-seq(0.05,5,0.05)
betweenness(g)

# Draw nodes and save positions
d=degree(g) 
cols=setNames(colorRampPalette(c("blue", "red"))(length(unique(d))), unique(d))
locs <- layout.fruchterman.reingold(g)
plot(g, 
     layout=locs, 
     vertex.label=NA, 
     main="Original",
     vertex.color=cols[degree(g)])
g


vertex.color=degree(g)

test <- graph(c('ino','rvb','tip','rvb','swr','rvb'),directed = F)
degree(test)
ivi <- c(1.31,2.46,1.01,5.02)
ivg <- c(1,2,3,4)
cols=setNames(colorRampPalette(c("blue", "red"))((5.03-1.01)/0.01), seq(1.01,5.02,0.01))

a <- data.frame('id'=c('ino','rvb','tip','swr'),
                'ipvsi'=c(1.31,2.46,1.01,5.02))
V(test)$ipvsi <- as.numeric(a$ipvsi[match(V(test)$name,a$id)])

palf <- colorRampPalette(c("gray80", "dark red"))
#V(g)$color <- palf[V(test)$ipvsi]
range1.100 <- function(x){1 + 99*(x-min(x))/(max(x)-min(x))}
colr <- palf(10000);
V(test)$color <- colr[round(range1.100(as.numeric(V(test)$ipvsi)))]

plot(test,vertex.color=V(test)$color)
c_clust <- cluster_fast_greedy(test)
plot(c_clust,test)
V(test)

rm(list=ls())
test <- fread('complex_all.csv')
gr <- graph.adjacency(matrix)

a <- test[,c(2,1)]
a$val <- 1

matrix <- matrix(0,nrow = length(unique(c(test$compx,test$name))),
                 ncol = length(unique(c(test$compx,test$name))))
colnames(matrix) <- unique(c(test$compx,test$name))
rownames(matrix) <- unique(c(test$compx,test$name))
for (i in 1:nrow(test)){
  matrix[test$name[i],test$compx[i]] <- 1
  matrix[test$compx[i],test$name[i]] <- 1
}
gr <- graph.adjacency(matrix,mode='undirected')
palf <- colorRampPalette(c("white", "blue", "yellow","red"))
#V(g)$color <- palf[V(test)$ipvsi]
range1.100 <- function(x){1 + 99*(x-min(x))/(max(x)-min(x))}
colr <- palf(100);
V(gr)$ipvsi <- as.numeric(test$ipvsinput[match(V(gr)$name,test$name)])
V(gr)$ipvsi[is.na(V(gr)$ipvsi)] <- 0
V(gr)$color <- colr[range1.100(round(as.numeric(V(gr)$ipvsi)))]
pdf('test.pdf')
plot(gr,
     #layout=locs, 
     vertex.size=5, 
     vertex.label.cex=0.2,
     main="Original",
     vertex.color=V(gr)$color)
dev.off()
legend("bottomleft", inset=.02, title="Number of Cylinders",
      'ys', fill=colr, horiz=TRUE, cex=0.1)