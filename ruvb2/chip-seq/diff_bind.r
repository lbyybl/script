#--QC for chip-seq
#source("http://bioconductor.org/biocLite.R")
#biocLite("ChIPQC")
require(DiffBind)
#library(ChIPQC)
library(BiocParallel)
library(DiffBind)
register(SerialParam())
data_blacklist <- read.table('/DATA/work/lbyybl/wh/ruvb2/chip-seq/sample190312/local_result/R2-IAA-Pol2-C1_FKDL190724841-1a/bigwig/R2-IAA-Pol2-C1_blacklist')
data_blacklist <- GRanges(
  seqnames = Rle(data_blacklist$V1),
  ranges = IRanges(data_blacklist$V2,data_blacklist$V3)#,
  #  strand = Rel(strand("."))
)
#--- diffbind
#---一旦读入了peaksets，合并函数就找到所有重叠的peaks，并导出一致性的peaksets。
#file <- 'r1doxsample.csv'
#name <- 'R21dox'
diffbind <- function(file,name){
  dbObj <- dba(sampleSheet=file)
  #---计算每个peaks/regions的count信息。先对一致性的peaks数据集进行标准化，然后根据他们的峰值（point of greatest read overlap）
  #---再次中心化并修剪一些peaks，最终得到更加标准的peak间隔。使用函数dba.count()，参数bUseSummarizeOverlaps可以得到更加标准的计算流程。
  dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE,filter = 10)
  # dba.plotPCA(dbObj,  attributes=DBA_CONDITION, label=DBA_ID)
  # plot(dbObj)
  # olaps <- dba.overlap(dbObj,mode=DBA_OLAP_RATE)
  # 
  # dba.plotVenn(dbObj,c(1:4))
  #---差异分析
  # Establishing a contrast 
  dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION,minMembers = 2)
  dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
  #  summary of results
  # dba.show(dbObj, bContrasts=T)
  #  overlapping peaks identified by the two different tools (DESeq2 and edgeR)
  # dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
  
  comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
  comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
  
  # EdgeR
  out <- as.data.frame(comp1.edgeR)
  write.table(out, file=paste0("diff_peak/",name,"_edgeR.txt"), sep="\t", quote=F, col.names = NA)
  # DESeq2
  out <- as.data.frame(comp1.deseq)
  write.table(out, file=paste0("diff_peak/",name,"_deseq2.txt"), sep="\t", quote=F, col.names = NA)
  
  # Create bed files for each keeping only significant peaks (p < 0.05)
  # EdgeR
  out <- as.data.frame(comp1.edgeR)
  edge.bed <- out[ which(out$FDR < 0.05), 
                   c("seqnames", "start", "end", "strand", "Fold")]
  write.table(edge.bed, file=paste0("diff_peak/",name,"_edgeR_sig.bed"), sep="\t", quote=F, row.names=F, col.names=F)
  
  # DESeq2
  out <- as.data.frame(comp1.deseq)
  deseq.bed <- out[ which(out$FDR < 0.05), 
                    c("seqnames", "start", "end", "strand", "Fold")]
  write.table(deseq.bed, file=paste0("diff_peak/",name,"_deseq2_sig.bed"), sep="\t", quote=F, row.names=F, col.names=F)
  
  #setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/Ruvbl1/sample20190509/graph/compare/diff_peak')
  pdf(paste0(name,'MA_plot.pdf'),width = 5,height = 5)
  dba.plotMA(dbObj,bXY=F,xrange = c(3,9),yrange = c(-4,4),bSmooth = T,fold = log(2))
  dev.off()
  pdf(paste0(name,'scatter_plot.pdf'),width = 5,height = 5.5)
  dba.plotMA(dbObj,bXY=T,xrange = c(2,11),yrange = c(2,11),bSmooth = T,fold = log(2))
  dev.off()
  pdf(paste0(name,'volcano_plot.pdf'),width = 7,height = 5.5)
  dba.plotVolcano(dbObj,bLabels = F)
  dev.off()
  pdf(paste0(name,'heatmap.pdf'))
  dba.plotHeatmap(dbObj,contrast=1,correlations=FALSE)
  dev.off()
  return(dbObj)
}
setwd('/WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl/Pol2')
sample <- fread('sample.csv')
dox05h <- sample[1:4,]
fwrite(dox05h,'dox05hsample.csv',sep = ',')
R05h1 <- sample[3:6,]
fwrite(R05h1,'R05h1sample.csv',sep = ',')
doxh1 <- sample[c(1:2,5:6),]
fwrite(doxh1,'doxh1sample.csv',sep = ',')

dox05h<-diffbind('dox05hsample.csv','dox05h')
R05h1<-diffbind('R05h1sample.csv','R05h1')
Rdox<-diffbind('doxh1sample.csv','R1dox')

#-----------------------------------
# this script is used to clasified the protein to different class
setwd('/WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl/Pol2/diff_peak')
library(data.table)
library(dplyr)
#-------------------------------
# got different region for Ruvbl1/2 peak
#setwd('/WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl/Pol2/diff_peak')
#file <- 'R1dox_deseq2.txt'
#name <- 'R1dox'
filter_fdr <- function(file,name){
  data <- fread(file,header = T)
  #data2 <- fread('RUVB2_POL2_deseq2_sig.bed',header = F)
  names(data)
  data <- as.data.frame(data)
  unchange <- data %>%
    filter(FDR >= 0.05 | (FDR < 0.05 & abs(Fold) < log(2))) %>%
    dplyr::select(seqnames,start,end)
  #unchange[,.(seqnames,start,end)]
  down <- data %>% 
    filter(FDR <0.05 & Fold > log(2)) %>%
    dplyr::select(seqnames,start,end)
  up <- data %>% 
    filter(FDR < 0.05 & Fold < -log(2)) %>%
    dplyr::select(seqnames,start,end)
  nrow(unchange)
  nrow(down)
  nrow(up)
  fwrite(unchange,paste0(name,'unchange_deseq2.txt'),sep = '\t',col.names = F)
  fwrite(down,paste0(name,'first_high_deseq2.txt'),sep = '\t',col.names = F)
  fwrite(up,paste0(name,'second_high_deseq2.txt'),sep = '\t',col.names = F)
  
}

filter_fdr('dox05h_deseq2.txt','dox05h')
filter_fdr('R05h1_deseq2.txt','R05h1')
filter_fdr('R1dox_deseq2.txt','R1dox')

#--- 找出文件对应的基???
setwd('/WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl/Pol2/diff_peak/diff_region')
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)
library(ReactomePA)
annopeak_fuc <- function(peak_fle,name){
  peak <- readPeakFile(peak_fle)
  pdf(paste0(name,'peak_chr.pdf'),width = 6,height = 5)
  print(covplot(peak_fle,chrs = c(paste0('chr',c(1:19,'X','Y')))))
  dev.off()

  pdf(paste0(name,'promoter_hetmap.pdf'),width = 3,height = 7)
  pdf(peakHeatmap(peak_fle, TxDb=txdb, upstream=3000, downstream=3000, color="red"))
  dev.off()

  pdf(paste0(name,'promoter_profile.pdf'),width=5,height = 5)
  pdf(plotAvgProf2(peak_fle, TxDb=txdb, upstream=3000, downstream=3000,
               xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency"))
  dev.off()
  peakAnno <- annotatePeak(peak_fle, tssRegion=c(-2000, 2000),
                           TxDb=txdb, annoDb="org.Mm.eg.db")
  colors   <- rainbow(10)
  pdf(paste0(name,'pieplot_anno.pdf'),width = 6,height = 4)
  print(plotAnnoPie(peakAnno))
  dev.off()
  pdf(paste0(name,'bar_anno.pdf'),width = 6,height = 3)
  print(plotAnnoBar(peakAnno))
  dev.off()
  pdf(paste0(name,'venn_anno.pdf'),width = 5,height = 5)
  print(vennpie(peakAnno))
  dev.off()
  pdf(paste0(name,'upset_anno.pdf'),width = 1,height = 1)
  print(upsetplot(peakAnno))
  dev.off()
  pdf(paste0(name,'unset_venn_anno.pdf'),width = 12,height = 12)
  print(upsetplot(peakAnno, vennpie=TRUE))
  dev.off()
  pdf(paste0(name,'dis2tss_anno.pdf'),width = 6,height = 3)
  print(plotDistToTSS(peakAnno,
                title="Distribution of transcription factor-binding loci\nrelative to TSS"))
  dev.off()
  pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId,organism = 'mouse')
  head(pathway1, 2)
  gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
  enrich_gene <- data.frame('name'=gene)
  fwrite(enrich_gene,paste0(name,'enrich_gene.txt'),col.names = F)
  pathway2 <- enrichPathway(gene,organism = 'mouse')
  head(pathway2, 2)
  pdf(paste0(name,'seq2gene_enrich.pdf'),width = 8,height = 5)
  print(dotplot(pathway2))
  dev.off()
  pdf(paste0(name,'enrich.pdf'),width=8,height = 5)
  print(dotplot(pathway1))
  dev.off()
  
}
for (file in system('ls ../*high_deseq2.txt',intern = T)){
  name <- strsplit(file,"\\/")[[1]][2]
  
  annopeak_fuc(file,strsplit(name,'high')[[1]][1])
}
