## ----setup, message=FALSE, echo = FALSE-------------------------------------------------
library(BiocStyle)
library(knitr)
library(clusterProfiler)
options(digits=3)
options(width=90)
setwd('/DATA2/work/lbyybl/coorlaborate/YB/NC_hbv_RNA/output/graph')
## ----setup2, message=FALSE, eval=TRUE---------------------------------------------------
library(limma)
library(Glimma)
library(edgeR)
#library(Mus.musculus)
library(org.Hs.eg.db)
## ----import1----------------------------------------------------------------------------
YB_file <- '/DATA2/work/lbyybl/coorlaborate/YB/NC_hbv_RNA/output/transcript_count_matrix.csv'
YB_data <- read.csv(YB_file, stringsAsFactors = F)
YB_data2 <- as.matrix(YB_data[,2:7])
rownames(YB_data2) <- as.vector(YB_data[,1])
groups <- as.factor(c('infect','infect','infect','noinfect','noinfect','noinfect'))
YB_file <- DGEList(counts=YB_data2,group=groups,remove.zeros=T)
samplenames <- c('infect1','infect2','infect3','noinfect1','noinfect2','noinfect3')
## ----import2----------------------------------------------------------------------------
x <- YB_file
class(x)
dim(x)
#rownames <- gsub('MSTRG.','',rownames(x))

#rownames(x) <- rownames
## ----annotatesamples--------------------------------------------------------------------

lane <- as.factor(rep(c("time1",'time2','time3'), each=2))
x$samples$lane <- lane
x$samples

## ----annotategenes, message=FALSE-------------------------------------------------------
geneid <- rownames(x)
# genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"),
#                 keytype="REFSEQ")
genes <- bitr(geneid, fromType = "REFSEQ", 
              toType = c("ENSEMBL", "SYMBOL",'ENTREZID'),
              OrgDb = org.Hs.eg.db)

head(genes)

# ## ----removedups-------------------------------------------------------------------------
genes <- genes[!duplicated(genes$REFSEQ),]

# ## ----assigngeneanno---------------------------------------------------------------------
x$genes <- genes
x

## ----cpm--------------------------------------------------------------------------------
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE, prior.count=2)

## ----lcpm-------------------------------------------------------------------------------
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm)

## ----filter-----------------------------------------------------------------------------
keep.exprs <- rowSums(cpm>0.5)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

## ----filterplot1, fig.height=4, fig.width=8, fig.cap="每个样本过滤前的原始数据（A）和过滤后（B）的数据的log-CPM值密度。竖直虚线标出了过滤步骤中所用阈值（相当于CPM值为约0.2）。"----
#pdf('gene_filter.pdf',width = 12,height = 6)
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
#dev.off()
## ----normalize--------------------------------------------------------------------------
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

## ----MDS1, fig.height=4, fig.width=8, fig.cap="log-CPM值在维度1和2的MDS图，以样品分组上色并标记（A）和维度3和4的MDS图，以测序道上色并标记（B）。图中的距离对应于最主要的倍数变化（fold change），默认情况下也就是前500个在每对样品之间差异最大的基因的平均（均方根）log2倍数变化。"----
#pdf('MDS_plot.pdf',width = 12,height = 6)
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- groups
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=groups, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")
#dev.off()
## ----GlimmaMDSplot----------------------------------------------------------------------
glMDSPlot(lcpm, labels=paste(groups, lane, sep="_"), 
          groups=x$samples[,c(1,4)], launch=FALSE)

## ----design-----------------------------------------------------------------------------
design <- model.matrix(~0+groups)
colnames(design) <- gsub("group", "", colnames(design))
design

## ----contrasts--------------------------------------------------------------------------
contr.matrix <- makeContrasts(
  InfectvsNoinfect = sinfect-snoinfect,
  levels = colnames(design))
contr.matrix

## ----voom, fig.height=4, fig.width=8, fig.cap="图中绘制了每个基因的均值（x轴）和方差（y轴），显示了在该数据上使用`voom`前它们之间的相关性（左），以及当运用`voom`的精确权重后这种趋势是如何消除的（右）。左侧的图是使用`voom`函数绘制的，它为进行log-CPM转换后的数据拟合线性模型从而提取残差方差。然后，对方差取平方根（或对标准差取平方根），并相对每个基因的平均表达作图。均值通过平均计数加上2再进行log2转换计算得到。右侧的图使用`plotSA`绘制了log2残差标准差与log-CPM均值的关系。平均log2残差标准差由水平蓝线标出。在这两幅图中，每个黑点表示一个基因，红线为对这些点的拟合。"----
#pdf('voom_plot.pdf',width = 12,height = 6)
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
#dev.off()
## ----decidetests------------------------------------------------------------------------
summary(decideTests(efit))

## ----treat------------------------------------------------------------------------------
tfit <- treat(vfit, lfc=0)
dt <- decideTests(tfit)
summary(dt)

## ----venn, fig.height=6, fig.width=6, fig.cap="韦恩图展示了仅basal和LP（左）、仅basal和ML（右）的对比的DE基因数量，还有两种对比中共同的DE基因数量（中）。在任何对比中均不差异表达的基因数量标于右下。"----
de.common <- which(dt[,1]!=0)
length(de.common)
head(tfit$genes$SYMBOL[de.common], n=20)
#pdf('venn_plot.pdf',width = 16,height = 6)
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))
#dev.off()
write.fit(tfit, dt, file="results.txt")

## ----toptables--------------------------------------------------------------------------
Infect.vs.Noinfect <- topTreat(tfit, coef=1, n=Inf)
head(Infect.vs.Noinfect)

## ----MDplot, fig.keep='none'------------------------------------------------------------
#pdf('AC12vsAD38nodoxMA_plot.pdf',width = 10,height = 6)
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-5,15))
#dev.off()

## ----GlimmaMDplot-----------------------------------------------------------------------
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],anno = x,
         side.main="REFSEQ", counts=lcpm, groups=groups, launch=FALSE)

## ----heatmap, fig.height=8, fig.width=5, fig.cap="在basal和LP的对比中前100个DE基因log-CPM值的热图。经过缩放调整后，每个基因（每行）的表达均值为0，并且标准差为1。给定基因相对高表达的样本被标记为红色，相对低表达的样本被标记为蓝色。浅色和白色代表中等表达水平的基因。样本和基因已通过分层聚类的方法重新排序。图中显示有样本聚类的树状图。", message=FALSE----
library(gplots)
Infect.vs.Noinfect.topgenes <- Infect.vs.Noinfect$REFSEQ[1:200]
up <- Infect.vs.Noinfect[Infect.vs.Noinfect$logFC >3 & Infect.vs.Noinfect$P.Value < 0.1,]
down <- Infect.vs.Noinfect[Infect.vs.Noinfect$logFC < -3 & Infect.vs.Noinfect$P.Value < 0.1,]
Infect.vs.Noinfect.topgenes <- c(up$REFSEQ,down$REFSEQ)
i <- which(rownames(lcpm) %in% Infect.vs.Noinfect.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(lcpm[i,], scale="row",
#           labRow=rownames(v$E)[i], labCol=group, 
#           col=mycol, trace="none", density.info="none", 
#           margin=c(8,6), lhei=c(2,10), dendrogram="column")
#pdf('38nodoxvsac12heatmap.pdf',width = 6,height = 16)
heatmap.2(lcpm[i,], scale="row",
          distfun = dist,
          hclustfun = hclust,
          labCol=samplenames, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")
library(ggplot2)
newdata <- data.frame('gene' = c(rep(rownames(lcpm[i,]),6)),
                      'sample' = c(rep(colnames(lcpm[i,]),each=nrow(lcpm[i,]))),
                      'value' = c(rep(0,6*nrow(lcpm[i,]))))
for (g in (1:nrow(newdata))){
  newdata$value[g] <- lcpm[i,][newdata$gene[g],newdata$sample[g]]
}
ggplot(newdata,aes(y=gene,x=sample,fill=value)) + geom_raster()
a <- as.data.frame(lcpm[i,])
heatmap.2(as.matrix(a[order(-(a$infect1+a$infect2+a$infect3)),]),col=mycol,trace="none", density.info="none",
          Rowv=FALSE,
          Colv=FALSE)
heatmap.2(lcpm[i,],col=mycol,trace="none", density.info="none",
          margin=c(8,6), lhei=c(2,10), dendrogram="column")
Rowv=FALSE,
Colv=FALSE)

#----富集分析
down <- tfit$genes$REFSEQ[dt[,1]==-1]
up <- tfit$genes$REFSEQ[dt[,1]==1]
ego2 <- enrichGO(gene         = down,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'REFSEQ',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 1)
#pdf('AC12vsAD38noDox_downenrich.pdf',width = 9,height = 6)
dotplot(ego2, showCategory=30)
#dev.off()


#----deseq2
library(DESeq2)
library(dplyr)
setwd("/DATA2/work/lbyybl/coorlaborate/YB/NC_hbv_RNA/output/graph")
mycounts<-read.csv("/DATA2/work/lbyybl/coorlaborate/YB/NC_hbv_RNA/output/transcript_count_matrix.csv")
head(mycounts)
mycounts <- mycounts %>% 
  filter(sum(infect1+infect2+infect3+noinfect1+noinfect2+noinfect3)>10) %>%
  filter(infect1>1,infect2>1,infect3>1,noinfect1>1,noinfect2>1,noinfect3>1)
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,2:7]
head(mycounts)
condition <- factor(c(rep("infect",3),rep("noinfect",3)), levels = c("infect","noinfect"))
condition
colData <- data.frame(row.names=colnames(mycounts), condition)
#colData <- data.frame('rawname'=colnames(mycounts),'condition'=c(rep("infect",3),rep("noinfect",3)))
#colData <- as.matrix(colData)
#rownames(colData) <- colnames(mycounts)
colData
#----构建dds对象，开始DEseq流程；
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)
dds
#---总体结果查看
#res = results(dds, contrast=c("condition", "control", "treat")) #或下面命令
res= results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res) 
#所有结果先进行输出
write.csv(res,file="All_results.csv")
table(res$padj<0.05)
summary(res)
diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 3)
up_gene_deseq2 <-subset(res, padj < 0.05 & log2FoldChange > 2)
down_gene_deseq2 <-subset(res, padj < 0.05 & log2FoldChange < -2)
data.table::fwrite(as.data.frame(down_gene_deseq2),'down.tsv',quote = F,sep = '\t',col.names = F,row.names = T)
data.table::fwrite(as.data.frame(up_gene_deseq2),'up.tsv',quote = F,sep = '\t',col.names = F,row.names = T)

#diff_gene_deseq2 <- cbind(up_gene_deseq2,down_gene_deseq2)
diff_name <- c(rownames(up_gene_deseq2),rownames(down_gene_deseq2))
#down_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) < 0)
total_reads <- mycounts %>% 
  summarise_all(sum)
cpm_reads <- mycounts %>%
  mutate(infect1=infect1/total_reads$infect1,infect2=infect2/total_reads$infect2,infect3=infect3/total_reads$infect3,
         noinfect1=noinfect1/total_reads$noinfect1,noinfect2=noinfect2/total_reads$noinfect2,noinfect3=noinfect3/total_reads$noinfect3)
cpm_reads <- log(cpm_reads)*1e6
rownames(cpm_reads) <- rownames(mycounts)
mycounts <- cpm_reads
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
pdf('heatmap.pdf',width = 6,height = 12)
heatmap.2(as.matrix(mycounts[i,]), scale="row", trace="none", dendrogram="both", Rowv=as.dendrogram(cluster.row), Colv=as.dendrogram(cluster.col), col=mycol)

# heatmap.2(as.matrix(mycounts[i,]),col=mycol,trace="none",scale="row",
#          # Rowv=FALSE,
#           Colv="NA")
dev.off()
data <- as.matrix(mycounts[i,])[order(mycounts[i,1]+mycounts[i,2]+mycounts[i,3]-(mycounts[i,4]+mycounts[i,5]+mycounts[i,6])),]

heatmap.2(data,col=mycol,trace="none",scale="row",
          Rowv=FALSE,
          Colv=FALSE)

ego2 <- enrichGO(gene         = rownames(up_gene_deseq2),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'REFSEQ',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 1)
# ego <- enrichKEGG(
#   gene          = rownames(up_gene_deseq2),
#   keyType     = "kegg",
#   organism   = 'hsa',
#   pvalueCutoff      = 0.05,
#   pAdjustMethod     = "BH",
#   qvalueCutoff  = 0.05
# )
pdf('enrichment.pdf')
dotplot(ego2, showCategory=30)
dev.off()
# 
# geneid <- rownames(x)
# # genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"),
# #                 keytype="REFSEQ")
# genes <- bitr(geneid, fromType = "REFSEQ", 
#               toType = c("ENSEMBL", "SYMBOL",'ENTREZID'),
#               OrgDb = org.Hs.eg.db)
# keytypes(org.Hs.eg.db)