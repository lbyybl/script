library(encodeChIPqc)
library(BSgenome.Mmusculus.UCSC.mm10)
setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/bw_bam/Ruvbl2/bam/test')


#-------------------------------------------------function
#------ IPstrength
counts <- tagCount(
  samples=c("RuvbL2.bam", "input.bam"),
  org="Mmusculus", assembly="UCSC", version="mm10")
IPstrength(counts['RuvbL2.bam'],counts['input.bam'])

#---- calculateBindingCharacteristics
bc <- calculateBindingCharacteristics('RuvbL2.bam','input.bam')
pp <- phantomPeak(bc$chip.data, bc$input.data, bc$bc)

        # #---- calculateIDR and plotIDR
        # 
        # chipTuples <- calculateIDR(c('RuvbL2.bam','RuvbL2_2.bam'),
        #                            c('input.bam','input2.bam'),
        #                            'Mmusculus','UCSC','mm10')
        # chipTuples <- calculateIDR(c('../Ruvb2_rep1_uniqe.bam','../Ruvb2_rep2_uniqe.bam'),
        #                            c('../Input_rep1_uniqe.bam','../Input_rep2_uniqe.bam'),
        #                            'Mmusculus','UCSC','mm10')
        # 
        # for (chipTuples in chipTuples) {
        #   pdf(paste(basename(chipTuples$rep1),'_VS_',basename(chipTuples$rep2)))
        #   plotIDR(chipTuple)
        #   dev.off()
        # }

#--- callPeaks
bc <- calculateBindingCharacteristics('RuvbL2.bam','input.bam',read.len = 150)
peaks <- callPeaks(bc$chip.data, bc$input.data,bc$bc)

#--- frip
library(rtracklayer)
frip('RuvbL2.bam',import.bed('/DATA/work/lbyybl/wh/ruvb2/chip-seq/bw_bam/Ruvbl2/peak/consistInput_GFP.bed'))

#--- gAlignmentsToSppTags
library(GenomicAlignments)
chipSample <- readGAlignments('RuvbL2.bam',param = ScanBamParam(tag='NM'))
sppTags <- gAlignmentsToSppTags(chipSample)

        # #----getCrossConsistencyRatios
        # chipTuples <- calculateIDR(c('RuvbL2.bam','RuvbL2_2.bam'),
        #                            c('input.bam','input2.bam'),
        #                            'Mmusculus','UCSC','mm10')
        # getCrossConsistencyRatios(chipTuples)
        # 
        # #---- getSelfConsistencyRatios
        # chipTuples <- calculateIDR(c('RuvbL2.bam','RuvbL2_2.bam'),
        #                            c('input.bam','input2.bam'),
        #                            'Mmusculus','UCSC','hg19')
        # getSelfConsistencyRatios(chipTuples)

        # #---- mser
        # bc <- calculateBindingCharacteristics('RuvbL2.bam','input.bam',read.len = 150)
        # minser <- mser(bc$chip.data, bc$input.data,bc$bc)

#---- normalizeTagCount and samplesHeatmap	and tagCount and pcaPlot
counts <- tagCount(
  samples=c("RuvbL2.bam",'RuvbL2_2.bam', "input.bam",'input2.bam'),
  org="Mmusculus", assembly="UCSC", version="mm10")
nCounts <- normalizeTagCount(counts)
topHeatmap(nCounts)
samplesHeatmap(nCounts)
pcaPlot(nCounts,as.factor(c('r1','r2','i1','i2')))

#--- PBC
pbc <- PBC('RuvbL2.bam')


