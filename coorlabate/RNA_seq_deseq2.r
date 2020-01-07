#--- this script is used Deseq2 to analysis the differential of YB's RNA-seq
#----deseq2
library(DESeq2)
library(dplyr)
library(gplots)
library(clusterProfiler)
library(org.Hs.eg.db)
rm(list=ls())
setwd("/DATA2/work/lbyybl/coorlaborate/YB/YB_RNA_seq/total/output/graph")
mycounts<-read.csv("/DATA2/work/lbyybl/coorlaborate/YB/YB_RNA_seq/total/output/38-doxvs38joadox_rna.csv",row.names = 1)
#mycounts<-read.csv("/DATA2/work/lbyybl/coorlaborate/YB/YB_RNA_seq/total/output/38-doxvsAc12_rna.csv",row.names = 1)

head(mycounts)

head(mycounts)
condition <- factor(c(rep("38nodox",3),rep("38jiadox",3)), levels = c("38nodox","38jiadox"))
#condition <- factor(c(rep("38nodox",3),rep("AC12withHBV",3)), levels = c("38nodox","AC12withHBV"))
condition
colData <- data.frame(row.names=colnames(mycounts), condition)

colData
#----构建dds对象，开始DEseq流程；
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
dds
#---总体结果查看

res= results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res) 
#所有结果先进行输出
write.csv(res,file="All_results.csv")
table(res$padj<0.05)
summary(res)
diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 2)
up_gene_deseq2 <-subset(res, padj < 0.05 & log2FoldChange > 2)
down_gene_deseq2 <-subset(res, padj < 0.05 & log2FoldChange < -2)
data.table::fwrite(as.data.frame(down_gene_deseq2),'down.tsv',quote = F,sep = '\t',col.names = F,row.names = T)
data.table::fwrite(as.data.frame(up_gene_deseq2),'up.tsv',quote = F,sep = '\t',col.names = F,row.names = T)


diff_name <- c(rownames(up_gene_deseq2),rownames(down_gene_deseq2))

calcu_cpm <- function(data){
  total_data <- data %>% 
    summarise_all(sum)
  for (i in 1:ncol(data)){
    data[,i] <- (data[,i]+1)/total_data[,i]*1e6
  }
  data <- log(data)
  return(data)
}

mycounts <- calcu_cpm(mycounts)
dim(diff_gene_deseq2)
dim(up_gene_deseq2)
dim(down_gene_deseq2)
i <- which(rownames(mycounts) %in% diff_name)
mycounts[i,]
mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(as.matrix(mycounts[i,]), scale="row",
#           hclustfun = function(x) hclust(x,method = 'average'),
#           #labCol=samplenames[c(4:9)], 
#           col=mycol, trace="none", density.info="none", 
#           margin=c(8,6), lhei=c(2,10), dendrogram="column")
# library(ComplexHeatmap)
# pheatmap(as.matrix(mycounts[i,]), scale="row", cluster_rows=T, cluster_cols=T, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete", color=colors)

distance.row = dist(as.matrix(mycounts[i,]), method = "manhattan")
cluster.row = hclust(distance.row, method = "median")
distance.col = dist(t(as.matrix(mycounts[i,])), method = "manhattan")
cluster.col = hclust(distance.col, method = "median")
# pdf('AD38noDOXAD38withDOX_heatmap.pdf',width = 6,height = 12)
# heatmap.2(as.matrix(mycounts[i,]), scale="row", trace="none", 
#           dendrogram="both", Rowv=as.dendrogram(cluster.row), 
#           Colv=as.dendrogram(cluster.col), col=mycol,
#           keysize = 1
#           #colsep=c(key.title=NA, keysize=0.06,key.xlab="relative abundance")
#           )
# 
# dev.off()
read_data <- function(file,name){
  data <- fread(file,col.names = c('t_id','chr','strand','start','end','t_name','num_exons','length','gene_id','gene_name','cov','FPKM')) %>%
    dplyr::select(t_name,FPKM)
  names(data)[2] <- name
  return(unique(data))
}
merge_data <- function(file1, file2,file3){
  merge <- merge(file1,file2,by='t_name')
  merge <- unique(merge)
  merge <- merge(merge,file3,by='t_name')
  merge <- unique(merge)
  return(merge)
}
nodox38r1 <- read_data('/DATA2/work/lbyybl/coorlaborate/YB/YB_RNA_seq/total/output/ballgown/38-dox-1/t_data.ctab','nodoxr38r1')
nodox38r2 <- read_data('/DATA2/work/lbyybl/coorlaborate/YB/YB_RNA_seq/total/output/ballgown/38-dox-2/t_data.ctab','nodoxr38r2')
nodox38r3 <- read_data('/DATA2/work/lbyybl/coorlaborate/YB/YB_RNA_seq/total/output/ballgown/38-dox-3/t_data.ctab','nodoxr38r3')
nodox38 <- merge_data(nodox38r1,nodox38r2,nodox38r3)
jiadox38r1 <- read_data('/DATA2/work/lbyybl/coorlaborate/YB/YB_RNA_seq/total/output/ballgown/38jiadox-1/t_data.ctab','jiadoxr38r1')
jiadox38r2 <- read_data('/DATA2/work/lbyybl/coorlaborate/YB/YB_RNA_seq/total/output/ballgown/38jiadox-2/t_data.ctab','jiadoxr38r2')
jiadox38r3 <- read_data('/DATA2/work/lbyybl/coorlaborate/YB/YB_RNA_seq/total/output/ballgown/38jiadox-3/t_data.ctab','jiadoxr38r3')
jiadox38 <- merge_data(jiadox38r1,jiadox38r2,jiadox38r3)

dox38 <- merge(nodox38,jiadox38,by='t_name')
up_gene <- as.data.frame(up_gene_deseq2)
up_gene$t_name <- rownames(up_gene)
up_gene <- up_gene %>%
  arrange(lfcSE)
up_gene <- up_gene[1:100,]
down_gene <- as.data.frame(down_gene_deseq2)
down_gene$t_name <- rownames(down_gene)
down_gene <- down_gene %>%
  arrange(lfcSE)
down_gene <- down_gene[1:100,]
diff_name2 <- c(up_gene$t_name,down_gene$t_name)
dox38 <- unique(dox38) %>%
  filter(t_name %in% diff_name2) 
NM_001164267 <- dox38 %>%
  filter(t_name=='NM_001164267')
NM_001363780 <- dox38 %>%
  filter(t_name=='NM_001363780')
NM_006625 <- dox38 %>%
  filter(t_name=='NM_006625')
NM_009588 <- dox38 %>%
  filter(t_name=='NM_009588')
dox38 <- dox38 %>%
  filter(t_name!='NM_006625')
rownames(dox38) <- dox38$t_name
pheatmap(as.matrix(dox38[2:7]),color = mycol,clustering_method = 'mcquitty',scale = 'row',
         show_rownames = F,filename = 'AD38noDOXAD38withDOX_heatmap.pdf')
library(pheatmap)
breaksList = seq(-3, 3, by = 1)
#pdf('AD38noDOXAD38withDOX_heatmap.pdf',width = 6,height = 12)
data <- pheatmap(as.matrix(mycounts[i,]),color = mycol,clustering_method = 'ward',
         scale = 'row',filename = 'AD38noDOXAD38withDOX_heatmap.pdf')
# pheatmap(as.matrix(mycounts[i,]),color = mycol,breaks = c(1,2),clustering_method = 'ward')
# pheatmap(as.matrix(mycounts[i,]),color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
#          breaks = breaksList,clustering_method = 'ward')#,filename = 'AD38noDOXAD38withDOX_heatmap.tiff')

#dev.off()
data <- as.matrix(mycounts[i,])[order(mycounts[i,1]+mycounts[i,2]+mycounts[i,3]-(mycounts[i,4]+mycounts[i,5]+mycounts[i,6])),]

heatmap.2(data,col=mycol,trace="none",scale="row",
          Rowv=FALSE,
          Colv=FALSE)

ego2_CC <- enrichGO(gene         = rownames(up_gene_deseq2),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'REFSEQ',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 1)
ego2_BP <- enrichGO(gene         = rownames(up_gene_deseq2),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'REFSEQ',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 1)
ego2_MF <- enrichGO(gene         = rownames(up_gene_deseq2),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'REFSEQ',
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 1)

pdf('up_CC_enrichment.pdf')
dotplot(ego2_CC, showCategory=30)
dev.off()
pdf('up_BP_enrichment.pdf')
dotplot(ego2_BP, showCategory=30)
dev.off()
pdf('up_MF_enrichment.pdf')
dotplot(ego2_MF, showCategory=30)
dev.off()

ego2_CC <- enrichGO(gene         = rownames(down_gene_deseq2),
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'REFSEQ',
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1)
ego2_BP <- enrichGO(gene         = rownames(down_gene_deseq2),
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'REFSEQ',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1)
ego2_MF <- enrichGO(gene         = rownames(down_gene_deseq2),
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'REFSEQ',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1)

pdf('down_CC_enrichment.pdf')
dotplot(ego2_CC, showCategory=30)
dev.off()
pdf('down_BP_enrichment.pdf')
dotplot(ego2_BP, showCategory=30)
dev.off()
pdf('down_MF_enrichment.pdf')
dotplot(ego2_MF, showCategory=30)
dev.off()