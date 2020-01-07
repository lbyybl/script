#--QC for chip-seq
#source("http://bioconductor.org/biocLite.R")
#biocLite("ChIPQC")
require(DiffBind)
library(ChIPQC)
library(BiocParallel)
library(DiffBind)
register(SerialParam())
data_blacklist <- read.table('/DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed')
data_blacklist <- GRanges(
  seqnames = Rle(data_blacklist$V1),
  ranges = IRanges(data_blacklist$V2,data_blacklist$V3)#,
#  strand = Rel(strand("."))
)
setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/unGFP-chip')
samples <- read.csv('sample.csv',header = T)
## Create ChIPQC object
#source('/home/boyuanli/bashscript/bin/ruvb2/ATAC-seq/QC/ChIPQC.r')
chipObj <- ChIPQC(samples, annotation="mm10",blacklist = data_blacklist) 
## create ChIPQC report
ChIPQCreport(chipObj, reportName="ChIP_QC_report_RUVB2_POL2", reportFolder="ChIPQCreport")
