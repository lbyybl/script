library(HiCcompare)
options(stringsAsFactors = FALSE);
library(OmicCircos);
library(dplyr)
#data(UCSC.mm10.chr);
#head(UCSC.mm10.chr);
setwd('/DATA/work/lbyybl/YB/capture20180921HBV/hicpro_result/HiC_pro_output/18R153461/hic_results/matrix/18R153461/raw/500000')
dir()
mat1 <- read.table("18R153461_500000.matrix")
bed1 <- read.table("18R153461_500000_abs.bed")
dat1 <- hicpro2bedpe(mat1, bed1)
head(dat1)
head(dat1$trans)
interwithHBV <- dat1$trans %>%
  filter((chr1=="HBV" & chr2!="HBV") | (chr1!="HBV" & chr2=="HBV"))
tail(interwithHBV)
bed2 <- bed1 %>%
  mutate(the.v=NA,NO=NA) 
colnames(bed2) <- c("seg.name", "seg.Start", "seg.End", "the.v", "NO")
chrname <- c(paste0("chr",c(1:19,"X","Y")),"HBV")
hbv_mm20mix <- segAnglePo(bed2, seg=chrname)
dat1.num <- nrow(interwithHBV)
colors <- rainbow(22, alpha = 0.5)
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="");
circos(R=300, type="chr", cir=hbv_mm20mix, col=colors, print.chr.lab=TRUE, W=4);
circos(R=290, cir=hbv_mm20mix, W=20, mapping=interwithHBV, type="link.pg", lwd=2, col=sample(colors,15));

