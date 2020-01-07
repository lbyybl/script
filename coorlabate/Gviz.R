### R code from vignette source 'Gviz.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: init
###################################################
options(width=65)
library(xtable)
source(system.file("scripts/documentation.R", package="Gviz"))
xtabDetails <- details
addParTable <- function(class, skip=c("showTitle", "size", "background.title"), add=NULL)
{
    Parameters <- data.frame("Display Parameter"=names(xtabDetails[[class]]),
                             "Description"=xtabDetails[[class]], check.names=FALSE)
    align <- "lrp{5in}"
    if(!is.null(add)){
        Parameters <- cbind(Parameters, add)
        align <- c("lp{4in}", "lp{4in}", "lp{4in}", "lp{4in}")
    }
    Parameters <- Parameters[order(Parameters[,1]),]
    Parameters <- apply(Parameters, 2, function(x) gsub("_", "\\_", x, fixed=TRUE))
    rownames(Parameters) <-  gsub("_", "\\_", rownames(Parameters), fixed=TRUE)
    sel <- Parameters[,1] %in% skip
    Parameters[,2] <- gsub("\\code{\\linkS4class{", "\\Rclass{{", Parameters[,2], fixed=TRUE)
    print(xtable(Parameters[!sel,], align=align), sanitize.text.function=function(x) x, include.rownames=FALSE,
          floating=FALSE, tabular.environment="longtable")
    return(invisible())
}

hasUcscConnection <- !is(try(rtracklayer::browserSession(), silent=TRUE), "try-error")
oto <- options(timeout=5)
hasBiomartConnection <- (!is(try(download.file("http://www.biomart.org", tempfile(), quiet=TRUE)), "try-error") &&
                         !is(try(biomaRt::listMarts(), silent=TRUE), "try-error"))
options(timeout=oto)
## Uncommenting this helps when the UCSC server has a hickup but still lets you connect:
##hasUcscConnection <- !is(try(rtracklayer::browserSession(), silent=TRUE), "try-error") && !is(try(IdeogramTrack(genome="hg19", chromosome=7), silent=TRUE), "try-error")



###################################################
### code chunk number 3: loadPackage
###################################################
library(Gviz)


###################################################
### code chunk number 4: AnnotationTrack
###################################################
library(GenomicRanges)
data(cpgIslands)
class(cpgIslands)
chr <- as.character(unique(seqnames(cpgIslands)))
gen <- genome(cpgIslands)
atrack <- AnnotationTrack(cpgIslands, name="CpG")


###################################################
### code chunk number 5: plotAnnotationTrack
###################################################
plotTracks(atrack)


###################################################
### code chunk number 6: GenomeAxisTrack
###################################################
gtrack <- GenomeAxisTrack()


###################################################
### code chunk number 7: plotGenomeAxisTrack
###################################################
plotTracks(list(gtrack, atrack))


###################################################
### code chunk number 8: showIdeogramTrack (eval = FALSE)
###################################################
## itrack <- IdeogramTrack(genome=gen, chromosome=chr)


###################################################
### code chunk number 9: doIdeogramTrack
###################################################
if(hasUcscConnection){
    itrack <- IdeogramTrack(genome=gen, chromosome=chr)
}else{
    data(itrack)
}


###################################################
### code chunk number 10: plotIdeogramTrack
###################################################
plotTracks(list(itrack, gtrack, atrack))

plotTracks(list(gtrack, atrack, itrack))
###################################################
### code chunk number 11: GeneRegionTrack
###################################################
data(geneModels) ## 此处的geneModels是data.frame
grtrack <- GeneRegionTrack(geneModels, genome=gen, chromosome=chr, name="Gene Model")
plotTracks(list(itrack, gtrack, atrack, grtrack))


###################################################
### code chunk number 12: zooming
###################################################
plotTracks(list(itrack, gtrack, atrack, grtrack), from=26700000, to=26750000)


###################################################
### code chunk number 13: zooming2
###################################################
plotTracks(list(itrack, gtrack, atrack, grtrack), extend.left=0.5, extend.right=1000000)


###################################################
### code chunk number 14: zooming3
###################################################
plotTracks(list(itrack, gtrack, atrack, grtrack), extend.left=0.5, extend.right=1000000, col=NULL)


###################################################
### code chunk number 15: zooming4
###################################################
library(BSgenome.Hsapiens.UCSC.hg19)
strack <- SequenceTrack(Hsapiens, chromosome=chr)
plotTracks(list(itrack, gtrack, atrack, grtrack, strack), from=26591822, to=26591852, cex=0.8)


###################################################
### code chunk number 16: DataTrack
###################################################
set.seed(255)
lim <- c(26700000, 26750000)
coords <- sort(c(lim[1], sample(seq(from=lim[1], to=lim[2]), 99), lim[2]))
dat <- runif(100, min=-10, max=10)
dtrack <- DataTrack(data=dat, start=coords[-length(coords)], end=coords[-1], chromosome=chr,
                    genome=gen, name="Uniform")
plotTracks(list(itrack, gtrack, atrack, grtrack, dtrack), from=lim[1], to=lim[2])


###################################################
### code chunk number 17: DataTrackHist
###################################################
plotTracks(list(itrack, gtrack, atrack, grtrack, dtrack), from=lim[1], to=lim[2], type="histogram")


###################################################
### code chunk number 18: displayPars1f
###################################################
grtrack <- GeneRegionTrack(geneModels, genome=gen, chromosome=chr, name="Gene Model",
                           transcriptAnnotation="symbol", background.title="brown")
head(displayPars(grtrack))
displayPars(grtrack) <- list(background.panel="#FFFEDB", col=NULL)
head(displayPars(grtrack))

plotTracks(list(itrack, gtrack, atrack, grtrack))


###################################################
### code chunk number 19: displayPars2f
###################################################
plotTracks(list(itrack, gtrack, atrack, grtrack),
           background.panel="#FFFEDB", background.title="darkblue")


###################################################
### code chunk number 20: displayPars3
###################################################
dp <- availableDisplayPars(grtrack)
tail(dp)


###################################################
### code chunk number 21: displayPars4
###################################################
getOption("Gviz.scheme")
scheme <- getScheme()
scheme$GeneRegionTrack$fill <- "salmon"
scheme$GeneRegionTrack$col <- NULL
scheme$GeneRegionTrack$transcriptAnnotation <- "transcript"
addScheme(scheme, "myScheme")
options(Gviz.scheme="myScheme")
grtrack <- GeneRegionTrack(geneModels, genome=gen, chromosome=chr, name="Gene Model")
plotTracks(grtrack)
options(Gviz.scheme="default")
grtrack <- GeneRegionTrack(geneModels, genome=gen, chromosome=chr, name="Gene Model",
                           transcriptAnnotation="symbol")


###################################################
### code chunk number 22: schemes (eval = FALSE)
###################################################
## .GvizSchemes <- list(myScheme=list(GeneRegionTrack=list(fill="salmon", col=NULL, transcriptAnnotation="transcript")))


###################################################
### code chunk number 23: plottingdirections
###################################################
plotTracks(list(itrack, gtrack, atrack, grtrack), reverseStrand=TRUE)


###################################################
### code chunk number 24: GenomeAxisTrackClass1
###################################################
axisTrack <- GenomeAxisTrack()
plotTracks(axisTrack, from=1e6, to=9e6)


###################################################
### code chunk number 25: GenomeAxisTrackClass2
###################################################
axisTrack <- GenomeAxisTrack(range=IRanges(start=c(2e6, 4e6), end=c(3e6, 7e6),
                                           names=rep("N-stretch", 2)))
plotTracks(axisTrack, from=1e6, to=9e6)


###################################################
### code chunk number 26: GenomeAxisTrackClass2a
###################################################
plotTracks(axisTrack, from=1e6, to=9e6, showId=TRUE)


###################################################
### code chunk number 27: GenomeAxisTrackClass3
###################################################
plotTracks(axisTrack, from=1e6, to=9e6, add53=TRUE, add35=TRUE)


###################################################
### code chunk number 28: GenomeAxisTrackClass4
###################################################
plotTracks(axisTrack, from=1e6, to=9e6, add53=TRUE, add35=TRUE, littleTicks=TRUE)


###################################################
### code chunk number 29: GenomeAxisTrackClass5
###################################################
plotTracks(axisTrack, from=1e6, to=9e6, exponent=4)


###################################################
### code chunk number 30: GenomeAxisTrackClass6
###################################################
plotTracks(axisTrack, from=1e6, to=9e6, labelPos="below")


###################################################
### code chunk number 31: GenomeAxisTrackClass7
###################################################
plotTracks(axisTrack, from=1e6, to=9e6, scale=0.5)


###################################################
### code chunk number 32: GenomeAxisTrackClass8
###################################################
plotTracks(axisTrack, from=1e6, to=9e6, scale=0.5, labelPos="below")


###################################################
### code chunk number 33: GenomeAxisTrackClassTable
###################################################
addParTable("GenomeAxisTrack")


###################################################
### code chunk number 34: IdeogramTrackClass1Show (eval = FALSE)
###################################################
## ideoTrack <- IdeogramTrack(genome="hg19", chromosome="chrX")
## plotTracks(ideoTrack, from=85e6, to=129e6)


###################################################
### code chunk number 35: IdeogramTrackClass1Do
###################################################
if(hasUcscConnection){
    ideoTrack <- IdeogramTrack(genome="hg19", chromosome="chrX")
}else{
    data(itrack)
}
plotTracks(ideoTrack, from=85e6, to=129e6)


###################################################
### code chunk number 36: IdeogramTrackClass2
###################################################
plotTracks(ideoTrack, from=85e6, to=129e6, showId=FALSE)


###################################################
### code chunk number 37: IdeogramTrackClass3
###################################################
plotTracks(ideoTrack, from=85e6, to=129e6, showId=FALSE, showBandId=TRUE, cex.bands=0.5)


###################################################
### code chunk number 38: IdeogramTrackClassTable
###################################################
addParTable("IdeogramTrack")


###################################################
### code chunk number 39: DataClass1
###################################################
data(twoGroups)
dTrack <- DataTrack(twoGroups, name="uniform")
plotTracks(dTrack)


###################################################
### code chunk number 40: <types
###################################################
types <- data.frame(Value=c("p", "l", "b", "a", "s", "S", "g", "r", "h", "confint", "smooth", "histogram", "mountain", "polygon", "boxplot", "gradient", "heatmap", "horizon"),
                    Type=c("dot plot", "lines plot", "dot and lines plot", "lines plot of average (i.e., mean) values", "stair steps (horizontal first)",
                           "stair steps (vertical first)", "add grid lines", "add linear regression line", "histogram lines", "confidence intervals for average values", "add loess curve",
                           "histogram (bar width equal to range with)", "'mountain-type' plot relative to a baseline",
                           "'polygon-type' plot relative to a baseline", "box and whisker plot",
                           "false color image of the summarized values", "false color image of the individual values",
                           "Horizon plot indicating magnitude and direction of a change relative to a baseline"))
print(xtable(types, align="lrp{5in}"), sanitize.text.function=function(x) x, include.rownames=FALSE,
          floating=FALSE, tabular.environment="longtable")


###################################################
### code chunk number 41: typePlots
###################################################
pushViewport(viewport(layout=grid.layout(nrow=9, ncol=2)))
i <- 1
for(t in types$Value)
{
    pushViewport(viewport(layout.pos.col=((i-1)%%2)+1, layout.pos.row=((i-1)%/%2)+1))
    if(t != "horizon"){
        names(dTrack) <- t
        plotTracks(dTrack, type=t, add=TRUE, cex.title=0.8, margin=0.5)
    }else{
        data(dtHoriz)
        names(dtHoriz) <- "horizon *"
        plotTracks(dtHoriz[8,], type="horizon", add=TRUE, cex.title=0.8, margin=0.5, showAxis=FALSE, horizon.origin=0.7)
    }
    i <- i+1
    popViewport(1)
}
popViewport(1)
names(dTrack) <- "uniform"


###################################################
### code chunk number 42: mutitype
###################################################
plotTracks(dTrack, type=c("boxplot", "a", "g"))


###################################################
### code chunk number 43: sampNames
###################################################
colnames(mcols(twoGroups))
plotTracks(dTrack, type=c("heatmap"), showSampleNames=TRUE, cex.sampleNames=0.6)


###################################################
### code chunk number 44: grouping
###################################################
plotTracks(dTrack, groups=rep(c("control", "treated"), each=3), type=c("a", "p", "confint"))


###################################################
### code chunk number 45: typeGroupedPlots
###################################################
pushViewport(viewport(layout=grid.layout(nrow=9, ncol=1)))
i <- 1
for(t in c("a", "s", "confint", "smooth", "histogram", "boxplot", "heatmap", "horizon"))
{
    pushViewport(viewport(layout.pos.col=((i-1)%%1)+1, layout.pos.row=((i-1)%/%1)+1))
    if(t != "horizon"){
        names(dTrack) <- t
        plotTracks(dTrack, type=t, add=TRUE, cex.title=0.8, groups=rep(1:2, each=3), margin=0.5)
    }else{
        plotTracks(dtHoriz[c(1,8),], type="horizon", add=TRUE, cex.title=0.8, margin=0.5, showAxis=FALSE, horizon.origin=0.3,
                   groups=1:2)
    }
    i <- i+1
    popViewport(1)
}
pushViewport(viewport(layout.pos.col=((i-1)%%1)+1, layout.pos.row=((i-1)%/%1)+1))
names(dTrack) <- "hor. hist."
plotTracks(dTrack, type="histogram", stackedBars=FALSE, add=TRUE, cex.title=0.8, groups=rep(1:2, each=3), margin=0.5)
popViewport(2)
names(dTrack) <- "uniform"


###################################################
### code chunk number 46: groupingLegend
###################################################
plotTracks(dTrack, groups=rep(c("control", "treated"), each=3), type=c("a", "p"), legend=TRUE)


###################################################
### code chunk number 47: horizLegend
###################################################
data(dtHoriz)
dtHoriz <- dtHoriz[1:6,]
plotTracks(dtHoriz, type="horiz", groups=rownames(values(dtHoriz)),
           showSampleNames=TRUE, cex.sampleNames = 0.6, separator=1)


###################################################
### code chunk number 48: filedt1
###################################################
bgFile <- system.file("extdata/test.bedGraph", package="Gviz")
dTrack2 <- DataTrack(range=bgFile, genome="hg19", type="l", chromosome="chr19", name="bedGraph")
class(dTrack2)
plotTracks(dTrack2)


###################################################
### code chunk number 49: filedt2
###################################################
library(rtracklayer)
dTrack3 <-  DataTrack(range=bgFile, genome="hg19", type="l", chromosome="chr19", name="bedGraph",
                      importFunction=function(file) import(con=file))
identical(dTrack2, dTrack3)


###################################################
### code chunk number 50: filedt3
###################################################
bamFile <- system.file("extdata/test.bam", package="Gviz")
dTrack4 <- DataTrack(range=bamFile, genome="hg19", type="l", name="Coverage", window=-1, chromosome="chr1")
class(dTrack4)
dTrack4
plotTracks(dTrack4, from=189990000, to=190000000)


###################################################
### code chunk number 51: filedt4
###################################################
plotTracks(dTrack4, chromosome="chr1", from=189891483, to=190087517)


###################################################
### code chunk number 52: filedt5
###################################################
myImportFun <- function(file, selection){
    ## do something here
}
DataTrack(range=bamFile, genome="hg19", type="l", name="Coverage", window=-1, chromosome="chr1",
          importFunction=myImportFun, stream=TRUE)


###################################################
### code chunk number 53: biggerdata
###################################################
dat <- sin(seq(pi, 10*pi, len=500))
dTrack.big <- DataTrack(start=seq(1,100000, len=500), width=15, chromosome="chrX",
                        genome="hg19", name="sinus",
                        data=sin(seq(pi, 5*pi, len=500))*runif(500, 0.5, 1.5))
plotTracks(dTrack.big, type="hist")


###################################################
### code chunk number 54: aggregation
###################################################
plotTracks(dTrack.big, type="hist", window=50)


###################################################
### code chunk number 55: aggregation2
###################################################
plotTracks(dTrack.big, type="hist", window=-1, windowSize=2500)


###################################################
### code chunk number 56: transformation
###################################################
plotTracks(dTrack.big, type="l", transformation=function(x){x[x<0] <- 0; x})


###################################################
### code chunk number 57: groupingAv1
###################################################
plotTracks(dTrack, groups=rep(c("control", "treated"), each=3), type=c("b"), aggregateGroups=TRUE)


###################################################
### code chunk number 58: groupingAv2
###################################################
plotTracks(dTrack, groups=rep(c("control", "treated"), each=3), type=c("b"), aggregateGroups=TRUE, aggregation="max")


###################################################
### code chunk number 59: DataTrackClassTable
###################################################
addParTable("DataTrack")


###################################################
### code chunk number 60: anntrack1
###################################################
aTrack <- AnnotationTrack(start=c(10, 40, 120), width=15, chromosome="chrX",
                          strand=c("+", "*", "-"),
                          id=c("Huey", "Dewey", "Louie"), genome="hg19", name="foo")
plotTracks(aTrack)


###################################################
### code chunk number 61: anntrack2
###################################################
plotTracks(aTrack, shape="box", featureAnnotation="id")


###################################################
### code chunk number 62: anntrack3
###################################################
plotTracks(aTrack, shape="ellipse", featureAnnotation="id", fontcolor.feature="darkblue")


###################################################
### code chunk number 63: anntrack4f
###################################################
aTrack.groups <- AnnotationTrack(start=c(50, 180, 260, 460, 860, 1240), width=c(15,20,40,100,200, 20),
                                 chromosome="chrX",
                                 strand=rep(c("+", "*", "-"), c(1,3,2)),
                                 group=rep(c("Huey", "Dewey", "Louie"), c(1,3,2)),
                                 genome="hg19", name="foo")
plotTracks(aTrack.groups, groupAnnotation="group")


###################################################
### code chunk number 64: anntrack4af
###################################################
plotTracks(aTrack.groups, groupAnnotation="group", just.group="right")


###################################################
### code chunk number 65: anntrack4bf
###################################################
plotTracks(aTrack.groups, groupAnnotation="group", just.group="above")


###################################################
### code chunk number 66: stacking1
###################################################
aTrack.stacked <- AnnotationTrack(start=c(50, 180, 260, 800, 600, 1240), width=c(15,20,40,100,500, 20),
                                 chromosome="chrX",
                                 strand="*",
                                 group=rep(c("Huey", "Dewey", "Louie"), c(1,3,2)),
                                 genome="hg19", name="foo")
plotTracks(aTrack.stacked, groupAnnotation="group")


###################################################
### code chunk number 67: stacking2
###################################################
plotTracks(aTrack.stacked, stacking="dense")


###################################################
### code chunk number 68: features
###################################################
feature(aTrack.stacked)
feature(aTrack.stacked)<- c("foo", "bar", "bar", "bar", "no", "no")


###################################################
### code chunk number 69: featuresIdPlot
###################################################
plotTracks(aTrack.stacked, featureAnnotation="feature", groupAnnotation="feature", fontcolor.feature=1, cex.feature=0.7)


###################################################
### code chunk number 70: featuresPlotf
###################################################
plotTracks(aTrack.stacked, groupAnnotation="group", foo="darkred", bar="darkgreen")


###################################################
### code chunk number 71: overplotting
###################################################
data("denseAnnTrack")
plotTracks(denseAnnTrack, showOverplotting=TRUE)


###################################################
### code chunk number 72: collapse1f
###################################################
data(collapseTrack)
plotTracks(ctrack)


###################################################
### code chunk number 73: collapse2f
###################################################
plotTracks(ctrack, min.width=1)


###################################################
### code chunk number 74: collapse3f
###################################################
plotTracks(ctrack, min.width=1, collapse=TRUE)


###################################################
### code chunk number 75: collapse4f
###################################################
plotTracks(ctrack, min.width=3, min.distance=5, collapse=TRUE)


###################################################
### code chunk number 76: collapse5f
###################################################
plotTracks(ctrack, min.width=3, min.distance=5, collapse=TRUE,
           mergeGroups=TRUE, extend.left=0.1)


###################################################
### code chunk number 77: fileat1
###################################################
aTrack2 <- AnnotationTrack(range=bamFile, genome="hg19", name="Reads", chromosome="chr1")
class(aTrack2)
aTrack2
plotTracks(aTrack2, from=189995000, to=190000000)


###################################################
### code chunk number 78: fileat2
###################################################
aTrack3 <- AnnotationTrack(range=bamFile, genome="hg19", name="Reads", chromosome="chr1", group="id")
aTrack3
plotTracks(aTrack3, from=189995000, to=190000000)


###################################################
### code chunk number 79: fileat3
###################################################
availableDefaultMapping(bamFile, "AnnotationTrack")


###################################################
### code chunk number 80: fileat4
###################################################
plotTracks(list(dTrack4, aTrack2), from=189990000, to=190000000)


###################################################
### code chunk number 81: AnnotationTrackClassTable
###################################################
addParTable("AnnotationTrack")


###################################################
### code chunk number 82: generegtrackf
###################################################
data(geneModels)
grtrack <- GeneRegionTrack(geneModels, genome=gen, chromosome=chr, name="foo")
head(gene(grtrack))
head(transcript(grtrack))
head(exon(grtrack))
head(symbol(grtrack))
plotTracks(grtrack)


###################################################
### code chunk number 83: generegtrack2af
###################################################
plotTracks(grtrack, transcriptAnnotation="symbol")


###################################################
### code chunk number 84: generegtrack2bf
###################################################
plotTracks(grtrack, transcriptAnnotation="transcript")


###################################################
### code chunk number 85: generegtrack2bf
###################################################
plotTracks(grtrack, exonAnnotation="exon", extend.left=-0.8, fontcolor.exon=1)


###################################################
### code chunk number 86: generegtrack3f
###################################################
plotTracks(grtrack, collapseTranscripts=TRUE, shape="arrow", transcriptAnnotation="symbol")


###################################################
### code chunk number 87: generegtrack3g
###################################################
plotTracks(grtrack, collapseTranscripts="longest", shape="arrow", transcriptAnnotation="symbol")


###################################################
### code chunk number 88: generegtrack3h
###################################################
plotTracks(grtrack, collapseTranscripts="meta", shape="arrow", transcriptAnnotation="symbol")


###################################################
### code chunk number 89: tdb2grt1
###################################################
library(GenomicFeatures)
samplefile <- system.file("extdata", "hg19_knownGene_sample.sqlite", package="GenomicFeatures")
txdb <- loadDb(samplefile)
GeneRegionTrack(txdb)


###################################################
### code chunk number 90: tdb2grt2
###################################################
txTr <- GeneRegionTrack(txdb, chromosome="chr6", start=35000000, end=40000000)


###################################################
### code chunk number 91: generegtrack4f
###################################################
feature(txTr)
plotTracks(txTr)


###################################################
### code chunk number 92: GeneRegionTrackClassTable
###################################################
addParTable("GeneRegionTrack")


###################################################
### code chunk number 93: BiomartGeneRegionTrackShow (eval = FALSE)
###################################################
## biomTrack <- BiomartGeneRegionTrack(genome="hg19", chromosome=chr, start=20e6, end=21e6,
##                                   name="ENSEMBL")
## plotTracks(biomTrack)


###################################################
### code chunk number 94: BiomartGeneRegionTrackDo
###################################################
if(hasBiomartConnection){
    biomTrack <- BiomartGeneRegionTrack(genome="hg19", chromosome=chr, start=20e6, end=21e6,
                                        name="ENSEMBL")
}else{
    data("biomTrack")
}
plotTracks(biomTrack)


###################################################
### code chunk number 95: BiomartGeneRegionTrackCol
###################################################
plotTracks(biomTrack, col.line=NULL, col=NULL)


###################################################
### code chunk number 96: BiomartGeneRegionTrackHeight
###################################################
plotTracks(biomTrack, col.line=NULL, col=NULL, stackHeight=0.3)


###################################################
### code chunk number 97: BiomartGeneRegionTrackFilterShow (eval = FALSE)
###################################################
## biomTrack <- BiomartGeneRegionTrack(genome="hg19", chromosome=chr, start=20e6, end=21e6,
##                                   name="ENSEMBL", filter=list(with_ox_refseq_mrna=TRUE))
## plotTracks(biomTrack, col.line=NULL, col=NULL, stackHeight=0.3)


###################################################
### code chunk number 98: BiomartGeneRegionTrackFilterDo
###################################################
if(hasBiomartConnection){
    biomTrack <- BiomartGeneRegionTrack(genome="hg19", chromosome=chr, start=20e6, end=21e6,
                                  name="ENSEMBL", filter=list(with_ox_refseq_mrna=TRUE))
    plotTracks(biomTrack, col.line=NULL, col=NULL, stackHeight=0.3)
}else{
    biomTrack@filter <- list(with_ox_refseq_mrna=TRUE)
    plotTracks(biomTrack, col.line=NULL, col=NULL, stackHeight=0.3)
    biomTrack@filter <- list()
}


###################################################
### code chunk number 99: BiomartGeneRegionTrackSymbolShow (eval = FALSE)
###################################################
## biomTrack <- BiomartGeneRegionTrack(genome="hg19", name="ENSEMBL", symbol="ABCB5")
## plotTracks(biomTrack)


###################################################
### code chunk number 100: BiomartGeneRegionTrackSymbolDo
###################################################
if(hasBiomartConnection){
    biomTrack <- BiomartGeneRegionTrack(genome="hg19", name="ENSEMBL", symbol="ABCB5")
    plotTracks(biomTrack, transcriptAnnotation="symbol")
}else{
    ranges(biomTrack) <- ranges(biomTrack)[symbol(biomTrack) == "ABCB5"]
    plotTracks(biomTrack, transcriptAnnotation="symbol")
}


###################################################
### code chunk number 101: BiomartGeneRegionTrackCustom
###################################################
library(biomaRt)
bm <- useMart(host="may2012.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",
              dataset="hsapiens_gene_ensembl")
fm <- Gviz:::.getBMFeatureMap()
fm[["symbol"]] <- "external_gene_id"
biomTrack <- BiomartGeneRegionTrack(genome="hg19", chromosome="chr7", start=20e6, end=21e6,name="ENSEMBL",
                                    featureMap=fm, biomart=bm)
plotTracks(biomTrack, col.line=NULL, col=NULL, stackHeight=0.3)


###################################################
### code chunk number 102: BiomartGeneRegionTrackClassTable
###################################################
addInfo <- t(data.frame(displayPars(biomTrack, names(details[["BiomartGeneRegionTrack"]]))))
colnames(addInfo) <- "Color"
addParTable("BiomartGeneRegionTrack", add=addInfo)


###################################################
### code chunk number 103: DetailsAnnotationTrack1
###################################################
library(GenomicRanges)
probes <- GRanges(seqnames="chr7", ranges=IRanges(start=c(2000000, 2070000, 2100000, 2160000), end=c(2050000, 2130000, 2150000, 2170000)),
                  strand=c("-", "+", "-", "-"))


###################################################
### code chunk number 104: DetailsAnnotationTrack2
###################################################
methylation <- matrix(c(rgamma(400, 1)), ncol=100, dimnames=list(paste("probe", 1:4, sep=""), NULL))
methylation[,51:100] <- methylation[,51:100] + 0:3
sgroups <- rep(c("grp1","grp2"), each=50)


###################################################
### code chunk number 105: DetailsAnnotationTrack3
###################################################
library(lattice)
details <- function(identifier, ...) {
    d <- data.frame(signal=methylation[identifier,], group=sgroups)
    print(densityplot(~signal, group=group, data=d, main=list(label=identifier, cex=0.7),
                      scales=list(draw=FALSE, x=list(draw=TRUE)), ylab="", xlab="",
                      ), newpage=FALSE, prefix="plot")
}


###################################################
### code chunk number 106: DetailsAnnotationTrack4
###################################################
deTrack <- AnnotationTrack(range=probes, genome="hg19", chromosome=7, id=rownames(methylation),
                           name="probe details", stacking="squish",
                           fun=details)
plotTracks(deTrack)


###################################################
### code chunk number 107: DetailsAnnotationTrack5
###################################################
selFun <- function(identifier, start, end, track, GdObject, ...){
    gcount <- table(group(GdObject))
    ## This computes the width of 2 pixels in genomic coordinates
    pxRange <- Gviz:::.pxResolution(min.width=20, coord="x")
    return((end-start)<pxRange && gcount[identifier]==1)
}


###################################################
### code chunk number 108: DetailsAnnotationTrack6
###################################################
detFun <- function(identifier, GdObject.original, ...){
    plotTracks(list(GenomeAxisTrack(scale=0.3, size=0.2, cex=0.7), GdObject.original[group(GdObject.original)==identifier]),
               add=TRUE, showTitle=FALSE)
}


###################################################
### code chunk number 109: DetailsAnnotationTrack7f
###################################################
data(geneDetails)
deTrack2 <- AnnotationTrack(geneDetails, fun=detFun, selectFun=selFun,
                            groupDetails=TRUE, details.size=0.5, detailsConnector.cex=0.5, detailsConnector.lty="dotted",
                            shape=c("smallArrow", "arrow"), groupAnnotation="group")
plotTracks(deTrack2, extend.left=90000)


###################################################
### code chunk number 110: DetailsAnnotationTrack5
###################################################
plotTracks(deTrack, details.size=0.75, detailsConnector.pch=NA, detailsConnector.col="darkred",
           detailsBorder.fill="#FFE3BF", detailsBorder.col="darkred", shape="box", detailsConnector.lty="dotted")


###################################################
### code chunk number 111: DetailsAnnotationTrackClassTableSec
###################################################
addParTable("DetailsAnnotationTrack")


###################################################
### code chunk number 112: SequenceTrack1
###################################################
library(BSgenome.Hsapiens.UCSC.hg19)
sTrack <- SequenceTrack(Hsapiens)
sTrack


###################################################
### code chunk number 113: SequenceTrack2
###################################################
plotTracks(sTrack, chromosome=1, from=20000, to=20050)


###################################################
### code chunk number 114: SequenceTrack3
###################################################
fcol <- c(A="darkgray", C="darkgray", T="darkgray", G="darkgray")
plotTracks(sTrack, chromosome=1, from=20000, to=20050, fontcolor=fcol)


###################################################
### code chunk number 115: SequenceTrack4
###################################################
plotTracks(sTrack, chromosome=1, from=20000, to=20050, add53=TRUE)


###################################################
### code chunk number 116: SequenceTrack5
###################################################
plotTracks(sTrack, chromosome=1, from=20000, to=20050, add53=TRUE, complement=TRUE)


###################################################
### code chunk number 117: SequenceTrack6
###################################################
plotTracks(sTrack, chromosome=1, from=20000, to=20100)


###################################################
### code chunk number 118: SequenceTrack7
###################################################
plotTracks(sTrack, chromosome=1, from=20000, to=201000)


###################################################
### code chunk number 119: SequenceTrack8
###################################################
plotTracks(sTrack, chromosome=1, from=20000, to=20100, cex=0.5)


###################################################
### code chunk number 120: SequenceTrackClassTableSec
###################################################
addParTable("SequenceTrack")


###################################################
### code chunk number 121: alignmentstrack_1_do
###################################################
afrom <- 2960000
ato <- 3160000
alTrack <- AlignmentsTrack(system.file(package="Gviz", "extdata", "gapped.bam"), isPaired=TRUE)
data(alTrackGenes)


###################################################
### code chunk number 122: alignmentstrack_1_show (eval = FALSE)
###################################################
## afrom <- 2960000
## ato <- 3160000
## alTrack <- AlignmentsTrack(system.file(package="Gviz", "extdata", "gapped.bam"), isPaired=TRUE)
## bmt <- BiomartGeneRegionTrack(genome="hg19", chromosome="chr12", start=afrom, end=ato,
##                               filter=list(with_ox_refseq_mrna=TRUE), stacking="dense")


###################################################
### code chunk number 123: alignmentstrack_2
###################################################
plotTracks(c(bmt, alTrack), from=afrom, to=ato, chromosome="chr12")

plotTracks(c(alTrack, itrack), from=afrom, to=ato, chromosome="chr12")
###################################################
### code chunk number 124: alignmentstrack_3
###################################################
plotTracks(c(bmt, alTrack), from=afrom, to=ato, chromosome="chr12", min.height=0, coverageHeight=0.08,
           minCoverageHeight=0)


###################################################
### code chunk number 125: alignmentstrack_4
###################################################
plotTracks(c(alTrack, bmt), from=afrom, to=ato, chromosome="chr12", type="coverage")


###################################################
### code chunk number 126: alignmentstrack_5
###################################################
plotTracks(c(bmt, alTrack), from=afrom+12700, to=afrom+15200, chromosome="chr12")


###################################################
### code chunk number 127: alignmentstrack_5_1
###################################################
plotTracks(c(bmt, alTrack), from=afrom+12700, to=afrom+15200, chromosome="chr12",
           type=c("coverage", "sashimi"))


###################################################
### code chunk number 128: alignmentstrack_5_2
###################################################
introns <- GRanges("chr12", IRanges(start=c(2973662, 2973919), end=c(2973848, 2974520)))
plotTracks(c(bmt, alTrack), from=afrom+12700, to=afrom+15200, chromosome="chr12",
           type=c("coverage", "sashimi"), sashimiFilter=introns)


###################################################
### code chunk number 129: alignmentstrack_5_3
###################################################
plotTracks(c(bmt, alTrack), from=afrom+12700, to=afrom+15200, chromosome="chr12",
           type=c("coverage", "sashimi"), sashimiFilter=introns, sashimiFilterTolerance=5L)


###################################################
### code chunk number 130: alignmentstrack_6
###################################################
plotTracks(c(bmt, alTrack), from=afrom+12700, to=afrom+15200, chromosome="chr12", reverseStacking=TRUE,
           col.mates="purple", col.gap="orange", type="pileup")


###################################################
### code chunk number 131: alignmentstrack_6_1
###################################################
alTrack <- AlignmentsTrack(system.file(package="Gviz", "extdata", "gapped.bam"), isPaired=FALSE)
plotTracks(c(bmt, alTrack), from=afrom+12700, to=afrom+15200, chromosome="chr12")


###################################################
### code chunk number 132: alignmentstrack_7
###################################################
afrom <- 44945200
ato <- 44947200
alTrack <- AlignmentsTrack(system.file(package="Gviz", "extdata", "snps.bam"), isPaired=TRUE)
plotTracks(alTrack, chromosome="chr21", from=afrom, to=ato)


###################################################
### code chunk number 133: alignmentstrack_8
###################################################
plotTracks(c(alTrack, sTrack), chromosome="chr21", from=afrom, to=ato)


###################################################
### code chunk number 134: alignmentstrack_9
###################################################
plotTracks(c(alTrack, sTrack), chromosome="chr21", from=44946590, to=44946660)


###################################################
### code chunk number 135: alignmentstrack_10
###################################################
plotTracks(c(alTrack, sTrack), chromosome="chr21", from=44946590, to=44946660, cex=0.5, min.height=8)


###################################################
### code chunk number 136: AlignmentsTrackClassTable
###################################################
addParTable("AlignmentsTrack")


###################################################
### code chunk number 137: ucscTrack1 (eval = FALSE)
###################################################
## from <- 65921878
## to <- 65980988
## knownGenes <- UcscTrack(genome="mm9", chromosome="chrX", track="knownGene", from=from, to=to,
##                         trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
##                         symbol="name", transcript="name", strand="strand", fill="#8282d2", name="UCSC Genes")


###################################################
### code chunk number 138: ucscTrack2 (eval = FALSE)
###################################################
## refGenes <- UcscTrack(genome="mm9", chromosome="chrX", track="xenoRefGene", from=from, to=to,
##                       trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
##                       symbol="name2", transcript="name", strand="strand", fill="#8282d2",
##                       stacking="dense", name="Other RefSeq")
## 
## ensGenes <- UcscTrack(genome="mm9", chromosome="chrX", track="ensGene", from=from, to=to,
##                       trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
##                       symbol="name2", transcript="name", strand="strand", fill="#960000",
##                       name="Ensembl Genes")


###################################################
### code chunk number 139: ucscTrack3 (eval = FALSE)
###################################################
## cpgIslands <- UcscTrack(genome="mm9", chromosome="chrX", track="cpgIslandExt", from=from, to=to,
##                         trackType="AnnotationTrack", start="chromStart", end="chromEnd", id="name",
##                         shape="box", fill="#006400", name="CpG Islands")
## 
## snpLocations <-  UcscTrack(genome="mm9", chromosome="chrX", track="snp128", from=from, to=to,
##                            trackType="AnnotationTrack", start="chromStart", end="chromEnd", id="name",
##                            feature="func", strand="strand", shape="box", stacking="dense", fill="black",
##                            name="SNPs")


###################################################
### code chunk number 140: ucscTrack4 (eval = FALSE)
###################################################
## conservation <- UcscTrack(genome="mm9", chromosome="chrX", track="Conservation", table="phyloP30wayPlacental",
##                           from=from, to=to, trackType="DataTrack", start="start", end="end", data="score",
##                           type="hist", window="auto", col.histogram="darkblue", fill.histogram="darkblue",
##                           ylim=c(-3.7, 4), name="Conservation")
## 
## gcContent <- UcscTrack(genome="mm9", chromosome="chrX", track="GC Percent", table="gc5Base",
##                        from=from, to=to, trackType="DataTrack", start="start", end="end", data="score",
##                        type="hist", window=-1, windowSize=1500, fill.histogram="black", col.histogram="black",
##                        ylim=c(30, 70), name="GC Percent")


###################################################
### code chunk number 141: ucscTrack5 (eval = FALSE)
###################################################
## axTrack <- GenomeAxisTrack()
## idxTrack <- IdeogramTrack(genome="mm9", chromosome="chrX")


###################################################
### code chunk number 142: ucscTrackLoad
###################################################
data(ucscItems)


###################################################
### code chunk number 143: ucscTrack6
###################################################
plotTracks(list(idxTrack, axTrack, knownGenes, refGenes, ensGenes, cpgIslands,
                gcContent, conservation, snpLocations), from=from, to=to, showTitle=FALSE)


###################################################
### code chunk number 144: HighlightTrack
###################################################
ht <- HighlightTrack(trackList=list(atrack, grtrack, dtrack), start=c(26705000, 26720000), width=7000, chromosome=7)
plotTracks(list(itrack, gtrack, ht), from = lim[1], to = lim[2])


###################################################
### code chunk number 145: HighlightTrack2
###################################################
ht1 <- HighlightTrack(trackList=list(itrack, gtrack, atrack), start=c(26705000, 26720000), width=7000, chromosome=7)
ht2 <- HighlightTrack(trackList=dtrack, start=c(26705000, 26720000), width=7000, chromosome=7)
plotTracks(list(ht1, grtrack, ht2), from=lim[1], to=lim[2])


###################################################
### code chunk number 146: HighlightTrackClassTable
###################################################
addParTable("HighlightTrack")


###################################################
### code chunk number 147: OverlayTrack
###################################################
dat <- runif(100, min=-2, max=22)
dtrack2 <- DataTrack(data=dat, start=coords[-length(coords)], end=coords[-1], chromosome=chr,
                     genome=gen, name="Uniform2", groups=factor("sample 2", levels=c("sample 1", "sample 2")),
                     legend=TRUE)
displayPars(dtrack) <- list(groups=factor("sample 1", levels=c("sample 1", "sample 2")), legend=TRUE)
ot <- OverlayTrack(trackList=list(dtrack2, dtrack))
ylims <- extendrange(range(c(values(dtrack), values(dtrack2))))
plotTracks(list(itrack, gtrack, ot), from=lim[1], to=lim[2], ylim=ylims, type=c("smooth", "p"))


###################################################
### code chunk number 148: OverlayTrack2
###################################################
displayPars(dtrack) <- list(alpha.title=1, alpha=0.5)
displayPars(dtrack2) <- list(alpha.title=1, alpha=0.5)
ot <- OverlayTrack(trackList=list(dtrack, dtrack2))
plotTracks(list(itrack, gtrack, ot), from=lim[1], to=lim[2], ylim=ylims, type=c("hist"), window=30)


###################################################
### code chunk number 149: multPlot1
###################################################
chroms <- c("chr1", "chr2", "chr3", "chr4")
maTrack <- AnnotationTrack(range=GRanges(seqnames=chroms, ranges=IRanges(start=1, width=c(100,400,200,1000)),
                                         strand=c("+", "+", "-", "+")), genome="mm9", chromosome="chr1", name="foo")

mdTrack <- DataTrack(range=GRanges(seqnames=rep(chroms, c(10, 40, 20, 100)),
                                   ranges=IRanges(start=c(seq(1,100,len=10), seq(1,400,len=40), seq(1, 200, len=20),
                                                          seq(1,1000, len=100)), width=9), values=runif(170)),
                     data="values", chromosome="chr1", genome="mm9", name="bar")


###################################################
### code chunk number 150: multPlot2
###################################################
mgTrack <- GenomeAxisTrack(scale=50, labelPos="below", exponent=3)
chromosome(itrack) <- "chr1"


###################################################
### code chunk number 151: multPlot3
###################################################
ncols <- 2
nrows <- length(chroms)%/%ncols
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrows, ncols)))
for(i in seq_along(chroms)){
    pushViewport(viewport(layout.pos.col=((i-1)%%ncols)+1, layout.pos.row=(((i)-1)%/%ncols)+1))
    plotTracks(list(itrack, maTrack, mdTrack, mgTrack), chromosome=chroms[i], add=TRUE)
    popViewport(1)
}


###################################################
### code chunk number 152: multPlot4
###################################################
library(lattice)
chroms <- data.frame(chromosome=chroms)
xyplot(1~chromosome|chromosome, data=chroms, panel=function(x){plotTracks(list(itrack , maTrack, mdTrack, mgTrack), chromosome=x, add=TRUE, showId=FALSE)},
       scales=list(draw=FALSE), xlab=NULL, ylab=NULL)


###################################################
### code chunk number 153: session-info
###################################################
  sessionInfo()


