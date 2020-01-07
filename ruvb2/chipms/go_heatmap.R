#--- draw the Gene ontology enrichment plot
setwd('/DATA/work/lbyybl/wh/ruvb2/chip-ms')
library(data.table)
library(dplyr)
library(pheatmap)

go_data <- fread('GO_enrich.csv',header = T,stringsAsFactors = F)
names(go_data) <- c('go_id','class','ref_n','enrich_n','expect','envsex','enrich_fold','pvalue','fdr','type')
go_data <- go_data %>%
  arrange(class,type,pvalue)
go_data %>%
  group_by(class)%>%
  summarise(n=n())
#go_data[43,]$class<-'chromatin remolder and histon modification'
which(go_data$class=='histon modification')
annotation_col2 <- go_data %>%
  select(class,type)
rownames(annotation_col2) <- go_data$go_id

enrich_data <- matrix(nrow = 1,ncol = nrow(go_data))

enrich_data[1,] <- go_data$pvalue
colnames(enrich_data) <- go_data$go_id
rownames(enrich_data) <- 'go'
#enrich_data['go2',] <- enrich_data['go',]
pheatmap(enrich_data,cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(enrich_data, annotation_col = annotation_col2, cluster_rows = FALSE, cluster_cols = FALSE)

ann_colors2 = list(
  type = c(cmp="powderblue", pathway="midnightblue",func="thistle")
  #CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
  #GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
  #class=c("white", "firebrick","#1B9E77","#7570B3")
)
enrich_data['go',] <- log(enrich_data['go',])
pheatmap(enrich_data, annotation_col = annotation_col2, cluster_rows = FALSE, cluster_cols = FALSE, 
         gaps_col = c(8, 14,38,43,62),annotation_colors = ann_colors2,
         color = colorRampPalette(c("violetred1","yellow","limegreen"))(50),
         filename = 'go_enrich.pdf')
