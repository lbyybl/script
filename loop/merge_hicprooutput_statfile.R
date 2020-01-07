#--- mgerge .RSstat file for hichipper to call loop merge
#--- for RSstat file
#--h3k27ac ko
setwd("/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/pol2_ko_hicpro/hic_results/data/pol2_ko")
#--- h3k27ac wt
setwd("/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/pol2_wt_hicpro/hic_results/data/pol2_wt")
dir()
file1 <- read.table("pol2_wt_1_mm10.bwt2pairs.RSstat",col.names = c("class","num"))
file2 <- read.table("pol2_wt_2_mm10.bwt2pairs.RSstat",col.names = c("class","num"))

total_file <- file1
total_file$num <- file1$num + file2$num

data.table::fwrite(total_file,file = "pol2_wt_mm10.bwt2pairs.RSstat",sep = "\t",col.names = FALSE)

#--- for pairstat file
#---h3k27ac ko
setwd("/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/pol2_ko_hicpro/bowtie_results/bwt2/pol2_ko")
#--- h3k27ac wt
setwd("/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/pol2_wt_hicpro/bowtie_results/bwt2/pol2_wt")
pairstatfile1 <- data.table::fread("pol2_wt_1_mm10.bwt2pairs.pairstat",col.names = c("class","num","percent"))
pairstatfile2 <- data.table::fread("pol2_wt_2_mm10.bwt2pairs.pairstat",col.names = c("class","num","percent"))

pairstat_total_file <- pairstatfile1
pairstat_total_file$num <- pairstatfile1$num + pairstatfile2$num
pairstat_total_file$percent <- round(pairstat_total_file$num/pairstat_total_file$num[1]*100,3)

data.table::fwrite(pairstat_total_file,file = "pol2_wt_mm10.bwt2pairs.pairstat",sep = "\t",col.names = FALSE)


#----------------------------------------------
#--- ocean-c
#--- mgerge .RSstat file for hichipper to call loop merge
#--- for RSstat file
#--- ocean-c
setwd("/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/pol2_oceanc_0h/hic_results/data/0h")
#--- 6h
setwd("/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/pol2_oceanc_6h/hic_results/data/6h")
dir()
file1 <- read.table("pol2_oceanc_6h_1_mm10.bwt2pairs.RSstat",col.names = c("class","num"))
file2 <- read.table("pol2_oceanc_6h_2_mm10.bwt2pairs.RSstat",col.names = c("class","num"))

total_file <- file1
total_file$num <- file1$num + file2$num

data.table::fwrite(total_file,file = "pol2_oceanc_6h_mm10.bwt2pairs.RSstat",sep = "\t",col.names = FALSE)

#--- for pairstat file
#--- ocean-c 0h
setwd("/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/pol2_oceanc_0h/bowtie_results/bwt2/0h")
#--- 6h
setwd("/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/pol2_oceanc_6h/bowtie_results/bwt2/6h")
pairstatfile1 <- data.table::fread("pol2_oceanc_6h_1_mm10.bwt2pairs.pairstat",col.names = c("class","num","percent"))
pairstatfile2 <- data.table::fread("pol2_oceanc_6h_2_mm10.bwt2pairs.pairstat",col.names = c("class","num","percent"))

pairstat_total_file <- pairstatfile1
pairstat_total_file$num <- pairstatfile1$num + pairstatfile2$num
pairstat_total_file$percent <- round(pairstat_total_file$num/pairstat_total_file$num[1]*100,3)

data.table::fwrite(pairstat_total_file,file = "pol2_oceanc_6h_mm10.bwt2pairs.pairstat",sep = "\t",,col.names = FALSE)
