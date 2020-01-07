## ---- echo=FALSE, results="hide", warning=FALSE, message=FALSE-------------
suppressPackageStartupMessages({
  library(ATACseqQC)
  library(ChIPpeakAnno)
  library(data.table)
  library(GenomicRanges)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(GenomicScores)
  gsco <- getGScores("phastCons60way.UCSC.mm10")
  #library(phastCons100way.UCSC.hg19)
  library(MotifDb)
  library(GenomicAlignments)
})
## ---- eval=FALSE-----------------------------------------------------------
#  library(BiocManager)
#  BiocManager::install(c("ATACseqQC", "ChIPpeakAnno", "MotifDb", "GenomicAlignments",
#             "BSgenome.Hsapiens.UCSC.hg19", "TxDb.Hsapiens.UCSC.hg19.knownGene",
#             "phastCons100way.UCSC.hg19"))

## --------------------------------------------------------------------------
## load the library
library(ATACseqQC)
setwd('/WORK/lbyybl/WH/rvb/ATAC/sample20191014/test/')
## input the bamFile from the ATACseqQC package 
bamfile <- 'test_sort_dup.bam'
bam_index <- paste0(bamfile,'.bai')
bamfile.labels <- gsub(".bam", "", basename(bamfile))

## ---- eval=FALSE-----------------------------------------------------------
#  source(system.file("extdata", "IGVSnapshot.R", package = "ATACseqQC"))

## --------------------------------------------------------------------------
QC_list <- bamQC(bamfile,bam_index,mitochondria = 'chrM', outPath=paste0(bamfile.labels,'_clean.bam'))
con <-  file(paste0(bamfile.labels,"mtcars.txt"), open = "w")
writeLines(text = paste0("totalQNAMEs ", QC_list$totalQNAMEs), con = con )
writeLines(text = paste0("duplicateRate ", QC_list$duplicateRate), con = con )
writeLines(text = paste0("mitochondriaRate ", QC_list$mitochondriaRate), con = con )
writeLines(text = paste0("properPairRate ", QC_list$properPairRate), con = con )
writeLines(text = paste0("unmappedRate ", QC_list$unmappedRate), con = con )
writeLines(text = paste0("hasUnmappedMateRate ", QC_list$hasUnmappedMateRate), con = con )
writeLines(text = paste0("notPassingQualityControlsRate ", QC_list$notPassingQualityControlsRate), con = con )
writeLines(text = paste0("nonRedundantFraction ", QC_list$nonRedundantFraction), con = con )
writeLines(text = paste0("PCRbottleneckCoefficient_1 ", QC_list$PCRbottleneckCoefficient_1), con = con )
writeLines(text = paste0("PCRbottleneckCoefficient_2 ", QC_list$PCRbottleneckCoefficient_2), con = con )
write.table(x=QC_list$MAPQ, file = con, quote = FALSE, sep = "\t")
write.table(x=QC_list$idxstats, file = con, quote = FALSE, sep = "\t")
close(con)

duplicate_freq <- readsDupFreq(bamfile,index = bam_index)
fwrite(as.data.frame(duplicate_freq),paste0(bamfile.labels,'duplicate_frep.txt'),sep = '\t')
pdf(paste0(bamfile.labels,'estimate_complexity.pdf'))
plot <- estimateLibComplexity(duplicate_freq)
print(plot)
fwrite(plot,paste0(bamfile.labels,'estimate_complexity.txt'),sep = '\t')
dev.off()

## --------------------------------------------------------------------------
## generate fragement size distribution
pdf(paste0(bamfile.labels,'_fragement_size.pdf'))
# fragSize <- fragSizeDist(bamfile, bamfile.labels)
print(fragSizeDist(bamfile,index = bam_index, bamFiles.labels = bamfile.labels))
dev.off()

## --------------------------------------------------------------------------
## bamfile tags to be read in
bamfile <- 'test_sort.bam'
bam_index <- 'test_sort.bam.bai'
outPath <- paste0(bamfile.labels,"_splited")
possibleTag <- combn(LETTERS, 2)
possibleTag <- c(paste0(possibleTag[1, ], possibleTag[2, ]),
                 paste0(possibleTag[2, ], possibleTag[1, ]))
library(Rsamtools)
bamTop100 <- scanBam(BamFile(bamfile, yieldSize = 100),
                     param = ScanBamParam(tag=possibleTag))[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100)>0]
tags
## files will be output into outPath
#setwd('/DATA/work/lbyybl/tools/ATACseqQC/splited/splited/')

dir.create(outPath)
## shift the coordinates of 5'ends of alignments in the bam file
#library(BSgenome.Hsapiens.UCSC.hg19)
seqlev <- paste0('chr',c((1:19),'X','Y')) ## subsample data for quick run
# which <- as(seqinfo(Mmusculus)[seqlev], "GRanges")
# which <- as(which,"IntegerRangesList")
# which <- as(seqinfo(Mmusculus), "GRanges")
# seqlev <- "chr1" ## subsample data for quick run
which <- as(seqinfo(Mmusculus)[seqlev], "GRanges")
gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
shiftedBamfile <- file.path(outPath, paste0(bamfile.labels,"shifted.bam"))
gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)
#export(gal1,shiftedBamfile)

## --------------------------------------------------------------------------
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# filter <- list(tx_chrom = paste0('chr',c(1:19,'X','Y')))
# txs <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene,filter=filter)
txs <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)
pt <- PTscore(gal1, txs)
pdf(paste0(bamfile.labels,'_PTscore.pdf'))
print(plot(pt$log2meanCoverage, pt$PT_score, 
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript"))
dev.off()

## --------------------------------------------------------------------------
nfr <- NFRscore(gal1, txs)
pdf(paste0(bamfile.labels,'_NFRscore.pdf'))
print(plot(nfr$log2meanCoverage, nfr$NFR_score, 
     xlab="log2 mean coverage",
     ylab="Nucleosome Free Regions score",
     main="NFRscore for 200bp flanking TSSs",
     xlim=c(-10, 0), ylim=c(-5, 5)))
dev.off()
## --------------------------------------------------------------------------
tsse <- TSSEscore(gal1, txs)
tsse_summary <- summary(tsse$TSS.enrichment.score)
tsse_summary <- as.data.frame(as.matrix(tsse_summary))
fwrite(tsse_summary,paste0(bamfile.labels,'_tsse.txt'),sep = '\t',col.names = F)
## --------------------------------------------------------------------------
#library(phastCons100way.UCSC.hg19)
## run program for chromosome 1 only
# txs <- txs[seqnames(txs) %in% paste0('chr',c(1:19,'X','Y'))]
#txs <- txs[seqnames(txs) %in% "chr1"]
genome <- Mmusculus
## split the reads into NucleosomeFree, mononucleosome, 
## dinucleosome and trinucleosome.
## and save the binned alignments into bam files.
objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = outPath,
                              conservation=gsco)
## list the files generated by splitGAlignmentsByCut.
dir(outPath)

## ----eval=FALSE------------------------------------------------------------
# objs <- splitBam(bamfile,index = bam_index, outPath=outPath,tags='YT',
#                 txs=txs, genome=genome,seqlev = paste0("chr",c(1:19,'X','Y')),
#                 conservation=gsco)

## ----fig.height=4, fig.width=4---------------------------------------------
library(ChIPpeakAnno)
bamfiles <- file.path(outPath,
                      c("NucleosomeFree.bam",
                        "mononucleosome.bam",
                        "dinucleosome.bam",
                        "trinucleosome.bam"))
## Plot the cumulative percentage of tag allocation in nucleosome-free 
## and mononucleosome bam files.
pdf(paste0(bamfile.labels,'_nucleosome_percentage.pdf'))
print(cumulativePercentage(bamfiles[1:2], as(seqinfo(Mmusculus), "GRanges")))
dev.off()
## ----fig.height=8, fig.width=4---------------------------------------------
TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)
## estimate the library size for normalization
(librarySize <- estLibSize(bamfiles))
librarySize2 <- as.data.frame(as.matrix(librarySize))
fwrite(librarySize2,paste0(bamfile.labels,'_libaray_size.txt'),sep = '\t')
## calculate the signals around TSSs.
NTILE <- 101
dws <- ups <- 1010
sigs <- enrichedFragments(gal=objs[c("NucleosomeFree", 
                                     "mononucleosome",
                                     "dinucleosome",
                                     "trinucleosome")], 
                          TSS=TSS,
                          librarySize=librarySize,
                          seqlev=seqlev,
                          TSS.filter=0.5,
                          n.tile = NTILE,
                          upstream = ups,
                          downstream = dws)
## log2 transformed signals
sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
#plot heatmap
pdf(paste0(bamfile.labels,'_nucleosome_heatmap.pdf'))
print(featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                      zeroAt=.5, n.tile=NTILE))
dev.off()
## ----fig.show="hide"-------------------------------------------------------
## get signals normalized for nucleosome-free and nucleosome-bound regions.
pdf(paste0(bamfile.labels,'_nucleosome_prifile.pdf'))
out <- featureAlignedDistribution(sigs, 
                                  reCenterPeaks(TSS, width=ups+dws),
                                  zeroAt=.5, n.tile=NTILE, type="l", 
                                  ylab="Averaged coverage")
dev.off()
## --------------------------------------------------------------------------
## rescale the nucleosome-free and nucleosome signals to 0~1
pdf(paste0(bamfile.labels,'_nucleosome_prifile_norm.pdf'))
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
out <- apply(out, 2, range01)
matplot(out, type="l", xaxt="n", 
        xlab="Position (bp)", 
        ylab="Fraction of signal")
axis(1, at=seq(0, 100, by=10)+1, 
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
dev.off()
## --------------------------------------------------------------------------
## foot prints
library(MotifDb)
CTCF <- query(MotifDb, c("CTCF"))
CTCF <- as.list(CTCF)
#print(CTCF[[1]], digits=2)
pdf(paste0(bamfile.labels,'_CTCF_footprint.pdf'))
sigs <- factorFootprints(shiftedBamfile, pfm=CTCF[[1]], 
                         genome=genome,
                         min.score="90%", seqlev='chr1',
                         upstream=100, downstream=100)
dev.off()
## ----fig.height=6, fig.width=6---------------------------------------------
pdf(paste0(bamfile.labels,'_CTCF_footprint_heatmap.pdf'))
featureAlignedHeatmap(sigs$signal, 
                      feature.gr=reCenterPeaks(sigs$bindingSites,
                                               width=200+width(sigs$bindingSites[1])), 
                      annoMcols="score",
                      sortBy="score",
                      n.tile=ncol(sigs$signal[[1]]))
dev.off()

sigs$spearman.correlation
sigs$Profile.segmentation

## --------------------------------------------------------------------------
pdf(paste0(bamfile.labels,'_CTCF_Vplot.pdf'))
vp <- vPlot(shiftedBamfile, pfm=CTCF[[1]], 
            genome=genome, min.score="90%", seqlev='chr1',
            upstream=200, downstream=200, 
            ylim=c(30, 250), bandwidth=c(2, 1))
dev.off()
pdf(paste0(bamfile.labels,'_CTCF_Vplot_dyad.pdf'))
distanceDyad(vp, pch=20, cex=.5)
dev.off()
## --------------------------------------------------------------------------
# path <- system.file("extdata", package="ATACseqQC", mustWork=TRUE)
# bamfiles <- dir(path, "*.bam$", full.name=TRUE)
# gals <- lapply(bamfiles, function(bamfile){
#   readBamFile(bamFile=bamfile, tag=character(0), 
#               which=GRanges("chr1", IRanges(1, 1e6)), 
#               asMates=FALSE)
# })
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
# library(GenomicAlignments)
# plotCorrelation(GAlignmentsList(gals), txs, seqlev="chr1")
# 
# ## ----sessionInfo-----------------------------------------------------------
# sessionInfo()

