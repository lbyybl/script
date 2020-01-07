options(stringsAsFactors = FALSE);
rm(list=ls())
library(OmicCircos);
setwd("/DATA/work/lbyybl/ypjiang/4c/4c20181204/rawdata/Hspa8-Nla3-Pol2_TKD181102142/test/all_bam/bigwig")
AC12withHBV <- read.table("AC12withHBV_local_1000.rmdump.bedgraph",header = FALSE)
AD38A <- read.table("AD38A_local_1000.rmdump.bedgraph",header = FALSE)

colors   <- rainbow(15, alpha=0.2)
data <- AD38A
pdf('AD38A_single_cirlcleB100.pdf',width = 12, height = 12)
par(mar=c(0, 0, 0, 0));

plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");

#circos(R=40, type="chr", cir="hg19", col=colors, print.chr.lab=TRUE, W=4, scale=TRUE);
circos(R=100, cir="hg19", W=200, mapping=data, col.v=4, type="l",   B=TRUE, col='black', lwd=0.1, scale=TRUE);
dev.off()

#-------------
colors   <- rainbow(15, alpha=0.9)
data <- AC12withHBV
pdf('AC12withHBV_single_cirlcleB1000.pdf',width = 12, height = 12)
par(mar=c(0, 0, 0, 0));

plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");

#circos(R=40, type="chr", cir="hg19", col=colors, print.chr.lab=TRUE, W=4, scale=TRUE);
circos(R=100, cir="hg19", W=200, mapping=data, col.v=4, type="l",   B=TRUE, col='black', lwd=0.1, scale=TRUE);
dev.off()
