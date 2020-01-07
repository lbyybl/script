## ----setup, message=FALSE, echo = FALSE-------------------------------------------------
library(BiocStyle)
library(knitr)
options(digits=3)
options(width=90)

## ----setup2, message=FALSE, eval=TRUE---------------------------------------------------
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
## ----import1----------------------------------------------------------------------------
YB_file <- '/DATA/work/lbyybl/YB/YB_gene_count_matrix.csv'
YB_data <- read.csv(YB_file, stringsAsFactors = F)
YB_data2 <- as.matrix(YB_data[,2:10])
rownames(YB_data2) <- as.vector(YB_data[,1])
groups <- as.factor(c('AD38_dox','AD38_dox','AD38_dox','AD38jiadox','AD38jiadox','AD38jiadox','Ac12','Ac12','Ac12'))
YB_file <- DGEList(counts=YB_data2,group=groups,remove.zeros=T)
samplenames <- c('38_dox_1','38_dox_2','38_dox_3','38jiadox_1','38jiadox_2','38jiadox_3','Ac12_1','Ac12_2','Ac12_3')
## ----import2----------------------------------------------------------------------------
x <- YB_file
class(x)
dim(x)
#rownames <- gsub('MSTRG.','',rownames(x))

#rownames(x) <- rownames
## ----annotatesamples--------------------------------------------------------------------

lane <- as.factor(rep(c("time1",'time2','time3'), 3))
x$samples$lane <- lane
x$samples

## ----annotategenes, message=FALSE-------------------------------------------------------
geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENTREZID")
head(genes)

## ----removedups-------------------------------------------------------------------------
genes <- genes[!duplicated(genes$ENTREZID),]

## ----assigngeneanno---------------------------------------------------------------------
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

## ----zeroes-----------------------------------------------------------------------------
table(rowSums(x$counts==0)==9)

## ----filter-----------------------------------------------------------------------------
keep.exprs <- filterByExpr(x, min.total.count=10,min.count=2)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

## ----filterplot1, fig.height=4, fig.width=8, fig.cap="每个样本过滤前的原始数据（A）和过滤后（B）的数据的log-CPM值密度。竖直虚线标出了过滤步骤中所用阈值（相当于CPM值为约0.2）。"----
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
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

## ----normalize--------------------------------------------------------------------------
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

## ----normalizemodifieddata--------------------------------------------------------------
x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

## ----plotmodifieddata, fig.height=4, fig.width=8, fig.cap="样例数据：log-CPM值的箱线图展示了未经归一化的数据（A）及归一化后的数据（B）中每个样本的表达分布。数据集经过调整，样本1和2中的表达计数分别被缩放到其原始值的5%和500%。"----
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")

## ----MDS1, fig.height=4, fig.width=8, fig.cap="log-CPM值在维度1和2的MDS图，以样品分组上色并标记（A）和维度3和4的MDS图，以测序道上色并标记（B）。图中的距离对应于最主要的倍数变化（fold change），默认情况下也就是前500个在每对样品之间差异最大的基因的平均（均方根）log2倍数变化。"----
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

## ----GlimmaMDSplot----------------------------------------------------------------------
glMDSPlot(lcpm, labels=paste(groups, lane, sep="_"), 
          groups=x$samples[,c(1,4)], launch=FALSE)

## ----design-----------------------------------------------------------------------------
design <- model.matrix(~0+groups+lane)
colnames(design) <- gsub("group", "", colnames(design))
design

## ----contrasts--------------------------------------------------------------------------
contr.matrix <- makeContrasts(
  AC12vsAD38noDox = sAc12-sAD38_dox, 
  AC12vsAD38jiaDox = sAc12-sAD38jiadox, 
  AD38noDoxvsAD38jiaDox = sAD38_dox-sAD38jiadox, 
  levels = colnames(design))
contr.matrix

## ----voom, fig.height=4, fig.width=8, fig.cap="图中绘制了每个基因的均值（x轴）和方差（y轴），显示了在该数据上使用`voom`前它们之间的相关性（左），以及当运用`voom`的精确权重后这种趋势是如何消除的（右）。左侧的图是使用`voom`函数绘制的，它为进行log-CPM转换后的数据拟合线性模型从而提取残差方差。然后，对方差取平方根（或对标准差取平方根），并相对每个基因的平均表达作图。均值通过平均计数加上2再进行log2转换计算得到。右侧的图使用`plotSA`绘制了log2残差标准差与log-CPM均值的关系。平均log2残差标准差由水平蓝线标出。在这两幅图中，每个黑点表示一个基因，红线为对这些点的拟合。"----
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

## ----decidetests------------------------------------------------------------------------
summary(decideTests(efit))

## ----treat------------------------------------------------------------------------------
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

## ----venn, fig.height=6, fig.width=6, fig.cap="韦恩图展示了仅basal和LP（左）、仅basal和ML（右）的对比的DE基因数量，还有两种对比中共同的DE基因数量（中）。在任何对比中均不差异表达的基因数量标于右下。"----
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
head(tfit$genes$SYMBOL[de.common], n=20)
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))
write.fit(tfit, dt, file="results.txt")

## ----toptables--------------------------------------------------------------------------
basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)
head(basal.vs.lp)
head(basal.vs.ml)

## ----MDplot, fig.keep='none'------------------------------------------------------------
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))

## ----GlimmaMDplot-----------------------------------------------------------------------
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],anno = x,
         side.main="ENTREZID", counts=lcpm, groups=groups, launch=FALSE)

## ----heatmap, fig.height=8, fig.width=5, fig.cap="在basal和LP的对比中前100个DE基因log-CPM值的热图。经过缩放调整后，每个基因（每行）的表达均值为0，并且标准差为1。给定基因相对高表达的样本被标记为红色，相对低表达的样本被标记为蓝色。浅色和白色代表中等表达水平的基因。样本和基因已通过分层聚类的方法重新排序。图中显示有样本聚类的树状图。", message=FALSE----
library(gplots)
basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(lcpm[i,], scale="row",
#           labRow=rownames(v$E)[i], labCol=group, 
#           col=mycol, trace="none", density.info="none", 
#           margin=c(8,6), lhei=c(2,10), dendrogram="column")
heatmap.2(lcpm[i,], scale="row",
          labRow=rownames(v$E)[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")


## ----camera-----------------------------------------------------------------------------
load(system.file("extdata", "mouse_c2_v5p1.rda", package = "RNAseq123"))
idx <- ids2indices(Mm.c2,id=rownames(v))
cam.BasalvsLP <- camera(v,idx,design,contrast=contr.matrix[,1])
head(cam.BasalvsLP,5)
cam.BasalvsML <- camera(v,idx,design,contrast=contr.matrix[,2])
head(cam.BasalvsML,5)
cam.LPvsML <- camera(v,idx,design,contrast=contr.matrix[,3])
head(cam.LPvsML,5)

## ----barcodeplot, fig.height=6, fig.width=6, fig.cap="`LIM_MAMMARY_LUMINAL_MATURE_UP` （红色条形，图表上方）和`LIM_MAMMARY_LUMINAL_MATURE_DN`（蓝色条形，图表下方）基因集在LP和ML的对比中的条码图，每个基因集都有一条富集线展示了竖直条形在图表每部分的相对富集程度。Lim等人的实验[@Lim:BreastCancerRes:2010]非常类似于我们的，用了相同的分选方式来获取不同的细胞群，只是他们使用的是微阵列而不是RNA-seq来测定基因表达。需要注意的是，上调基因集发生下调而下调基因集发生上调的逆相关性来自于对比的设定方式（LP相比于ML），如果将其对调，方向性将会吻合。"----
barcodeplot(efit$t[,3], index=idx$LIM_MAMMARY_LUMINAL_MATURE_UP, 
            index2=idx$LIM_MAMMARY_LUMINAL_MATURE_DN, main="LPvsML")

## ----softwareinfo-----------------------------------------------------------------------
sessionInfo()

