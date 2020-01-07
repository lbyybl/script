#--QC for chip-seq
#source("http://bioconductor.org/biocLite.R")
#biocLite("ChIPQC")
require(DiffBind)
library(ChIPQC)
library(BiocParallel)
library(DiffBind)
register(SerialParam())
data_blacklist <- read.table('/DATA/work/lbyybl/wh/ruvb2/chip-seq/sample190312/local_result/R2-IAA-Pol2-C1_FKDL190724841-1a/bigwig/R2-IAA-Pol2-C1_blacklist')
data_blacklist <- GRanges(
  seqnames = Rle(data_blacklist$V1),
  ranges = IRanges(data_blacklist$V2,data_blacklist$V3)#,
#  strand = Rel(strand("."))
)
setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/bw_bam/unique/bam')
samples <- read.table('sample.tsv',header = T)
## Create ChIPQC object
chipObj <- ChIPQC(samples, annotation="mm10",blacklist = data_blacklist) 
## create ChIPQC report
ChIPQCreport(chipObj, reportName="ChIP_QC_report_RUVB2_POL2", reportFolder="ChIPQCreport")

#--- diffbind
#---一旦读入了peaksets，合并函数就找到所有重叠的peaks，并导出一致性的peaksets。
setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/bw_bam/Pol2/bam')
dbObj <- dba(sampleSheet="sample.csv")
#---计算每个peaks/regions的count信息。先对一致性的peaks数据集进行标准化，然后根据他们的峰值（point of greatest read overlap）
#---再次中心化并修剪一些peaks，最终得到更加标准的peak间隔。使用函数dba.count()，参数bUseSummarizeOverlaps可以得到更加标准的计算流程。
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)
plot(dbObj)
olaps <- dba.overlap(dbObj,mode=DBA_OLAP_RATE)

dba.plotVenn(dbObj,c(1:4))
#---差异分析
# Establishing a contrast 
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
#  summary of results
dba.show(dbObj, bContrasts=T)
#  overlapping peaks identified by the two different tools (DESeq2 and edgeR)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)

comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)

# EdgeR
out <- as.data.frame(comp1.edgeR)
write.table(out, file="diff_peak/RUVB2_POL2_edgeR.txt", sep="\t", quote=F, col.names = NA)
# DESeq2
out <- as.data.frame(comp1.deseq)
write.table(out, file="diff_peak/RUVB2_POL2_deseq2.txt", sep="\t", quote=F, col.names = NA)

# Create bed files for each keeping only significant peaks (p < 0.05)
# EdgeR
out <- as.data.frame(comp1.edgeR)
edge.bed <- out[ which(out$FDR < 0.05), 
                 c("seqnames", "start", "end", "strand", "Fold")]
write.table(edge.bed, file="diff_peak/RUVB2_POL2_edgeR_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)

# DESeq2
out <- as.data.frame(comp1.deseq)
deseq.bed <- out[ which(out$FDR < 0.05), 
                  c("seqnames", "start", "end", "strand", "Fold")]
write.table(deseq.bed, file="diff_peak/RUVB2_POL2_deseq2_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)

setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/Ruvbl1/sample20190509/graph/compare/diff_peak')
pdf('MA_plot.pdf',width = 5,height = 5)
dba.plotMA(dbObj,bXY=F,xrange = c(3,11),yrange = c(-3,3),bSmooth = T,fold = log(2))
dev.off()
pdf('scatter_plot.pdf',width = 5,height = 5.5)
dba.plotMA(dbObj,bXY=T,xrange = c(3,11),yrange = c(3,11),bSmooth = T,fold = log(2))
dev.off()
pdf('volcano_plot.pdf',width = 7,height = 5.5)
dba.plotVolcano(dbObj,bLabels = F)
dev.off()
pdf('heatmap.pdf')
dba.plotHeatmap(dbObj,contrast=1,correlations=FALSE)
dev.off()