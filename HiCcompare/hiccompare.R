setwd("/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare")
## ----set-options, echo=FALSE, cache=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- HiCcompare不使用inter interaction的信息所以灭有办法找出染色体间的互做差异
options(width = 400)
options(stringsAsFactors = F)

library(data.table)

for (i in c(1:19,"X","Y")){
  wtdata <- fread(paste0("Wt_",i,".observed"))
  kodata <- fread(paste0("Ko_",i,".observed"))
  assign(paste0("chr",i,".table"), create.hic.table(wtdata, kodata))
}

## ---- eval = FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- 使用bam文件找CNV，bin.size要和你后面找difference使用的大小一致，CNV.level分四个水平，由四个数字代表
#  cnv <- get_CNV(path2bam = 'path/to/bamfiles', out.file = 'path/to/bamfiles/outfile',
#                 bin.size = 1000, genome = 'hg19', CNV.level = 2)

## ---- eval = FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- 基因组上有一些区域总是能在二代测序中有一些异常的信号，ENCODE中称其为blacklist
# data('hg19_blacklist')
# 
# # combine cnv excluded regions with blacklist regions
# exclude <- cbind(cnv, hg19_blacklist)

## ---- warning=FALSE, message=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(`HiCcompare`)
# load the data
data("HMEC.chr22")
data("NHEK.chr22")
head(HMEC.chr22)  # 貌似跟juicer提出来的矩阵很像

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# create the `hic.table` object
chr22.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22')
head(chr22.table)

## ---- eval = FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  chr22.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22', exclude.regions = exclude, exclude.overlap = 0.2)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# create list of `hic.table` objects
data("HMEC.chr10")
data("NHEK.chr10")

# create the `hic.table` object
chr10.table <- create.hic.table(HMEC.chr10, NHEK.chr10, chr = 'chr10')
hic.list <- list(chr10.table, chr22.table)
head(hic.list) # 将多条染色体放到了一个列表中

## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- 从bedpe文件中创建hic.talbe对象
HMEC.chr22_BEDPE <- chr22.table[, 1:7, with=FALSE]
NHEK.chr22_BEDPE <- chr22.table[, c(1:6, 8), with=FALSE]
head(HMEC.chr22_BEDPE)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bed.hic.tab <- create.hic.table(HMEC.chr22_BEDPE, NHEK.chr22_BEDPE)
head(bed.hic.tab)

## ---- eval=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- 从hicpro的结果中创建table.hic对象
#  # first dataset
#  mat1 <- read.table("hic1_1000000.matrix")
#  bed1 <- read.table("hic1_1000000_abs.bed")
#  dat1 <- hicpro2bedpe(mat, bed)
#  dat1 <- dat1$cis # extract intrachromosomal matrices
#  # second dataset
#  mat2 <- read.table("hic2_1000000.matrix")
#  bed2 <- read.table("hic2_1000000_abs.bed")
#  dat2 <- hicpro2bedpe(mat, bed)
#  dat2 <- dat2$cis # extract intrachromosomal matrices
#  
#  # for chr1
#  hic.table <- create.hic.table(dat1[[1]], dat2[[1]])
#  
#  # for all chromosomes
#  hic.list <- mapply(create.hic.table, dat1, dat2, SIMPLIFY = FALSE)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## ---- eval = FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  hic.list <- total_sum(hic.list)

## ---- fig.width=7, fig.height=4-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- normalize样本间的差异
#--- 单条染色体
# Jointly normalize data for a single chromosome
hic.table <- hic_loess(chr22.table, Plot = TRUE, Plot.smooth = FALSE)
knitr::kable(head(hic.table))

## ---- message=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- 多条染色体
# Multiple hic.tables can be processed in parallel by entering a list of hic.tables
hic.list <- hic_loess(hic.list, parallel = TRUE)

## ---- message=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- 参数筛选
filter_params(hic.table)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- 差异分析
hic.table <- hic_compare(hic.table, A.min = 15, adjust.dist = TRUE, p.method = 'fdr', Plot = TRUE)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::kable(head(hic.table))

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- 将结果转换成GRange对象以便于后续分析
IntSet <- make_InteractionSet(hic.table)

## ----fig.width=7, fig.height=4--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- sumulate difference hic数据
number_of_unitdistances <- 100 # The dimensions of the square matrix to be simualted
number_of_changes       <- 250 # How many cells in the matrix will have changes
i.range <- sample(1:number_of_unitdistances, number_of_changes, replace = TRUE) # Indexes of cells to have controlled changes
j.range <- sample(1:number_of_unitdistances, number_of_changes, replace = TRUE) # Indexes of cells to have controlled changes

#--- 返回的结果是一个列表
sim_results <- hic_simulate(nrow = number_of_unitdistances, medianIF = 50000, sdIF = 14000, powerlaw.alpha = 1.8, fold.change = 4, i.range = i.range, j.range = j.range, Plot = TRUE, alpha = 0.1)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
names(sim_results) # 查看列表的名字

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sims <- sim_matrix(nrow = number_of_unitdistances, medianIF = 50000, sdIF = 14000, powerlaw.alpha = 1.8, fold.change = 4, 
                   i.range = i.range, j.range = j.range)

## ---- fig.height=4, fig.width=7-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
MD.plot1(M = hic.table$M, D = hic.table$D, mc = hic.table$mc, smooth = TRUE)

## ---- fig.height=4, fig.width=7-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# no p-value coloring
MD.plot2(M = hic.table$adj.M, D = hic.table$D, smooth = FALSE)

# p-value coloring
MD.plot2(M = hic.table$adj.M, D = hic.table$D, hic.table$p.value, smooth = FALSE)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
full.NHEK <- sparse2full(NHEK.chr22)
full.NHEK[1:5, 1:5]

sparse.NHEK <- full2sparse(full.NHEK)
head(sparse.NHEK)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
KR.NHEK <- KRnorm(full.NHEK)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
SCN.NHEK <- SCN(full.NHEK)

## ---- message=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
result <- MA_norm(hic.table, Plot = TRUE)

## ---- eval = FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  library(HiCcompare)
#  library(BiocParallel)
#  
#  args = commandArgs(trailingOnly=TRUE)
#  
#  dat1 <- read.table(args[1], header=FALSE, col.names=c("chr1", "start1", "end1", "chr2", "start2", "end2", "IF"))
#  dat2 <- read.table(args[2], header=FALSE, col.names=c("chr1", "start1", "end1", "chr2", "start2", "end2", "IF"))
#  
#  dat1 <- dat1[dat1$chr1==dat1$chr2, ]
#  dat2 <- dat2[dat2$chr1==dat2$chr2, ]
#  
#  dat1 <- split(dat1, dat1$chr1)
#  dat2 <- split(dat2, dat2$chr1)
#  
#  hic.list <- mapply(create.hic.table, dat1, dat2, SIMPLIFY = FALSE, scale=FALSE)
#  
#  hic.list <- total_sum(hic.list)
#  
#  register(MulticoreParam(workers = 10), default = TRUE)
#  
#  hic.list <- hic_loess(hic.list, Plot=TRUE, parallel=TRUE)
#  hic.list <- hic_compare(hic.list, A.min = NA, adjust.dist = TRUE, p.method = 'fdr', Plot = TRUE, parallel=TRUE)
#  
#  hic.list <- do.call(rbind, hic.list)
#  
#  hic.list <- hic.list[hic.list$p.adj<0.05,]
#  
#  write.table(hic.list, args[3])

