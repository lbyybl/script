rm(list=ls())
setwd('/WORK/lbyybl/WH/rvb/cor_uniq/MNase/graph/diff')
R05dox<- fread('bam_R205hMNase_merge-bam_R2DoxMNase_merge.positions.integrative.xls')
names <- 'R05dox'

fold <- 1
R05dox_occupancy <- R05dox[,1:13]
R05dox_shift <- R05dox[,c(1:8,14:18)]
R05dox_fuzz <- R05dox[,c(1:8,19:23)]

R05dox_up_ocup <- R05dox_occupancy %>%
  dplyr::filter(smt_diff_FDR<=0.05 & smt_log2FC>fold)
R05dox_down_ocup <- R05dox_occupancy %>%
  dplyr::filter(smt_diff_FDR<=0.05 & smt_log2FC< -fold)
fwrite(R05dox_up_ocup,paste0(names,'up_ocup.txt'),sep = '\t')
fwrite(R05dox_down_ocup,paste0(names,'down_ocup.txt'),sep = '\t')

R05dox_up_shift <- R05dox_shift %>%
  dplyr::filter(point_diff_FDR<=0.05 & point_log2FC >fold)
R05dox_down_shift <- R05dox_shift %>%
  dplyr::filter(point_diff_FDR<=0.05 & point_log2FC < -fold)
fwrite(R05dox_up_shift,paste0(names,'up_shift.txt'),sep = '\t')
fwrite(R05dox_down_shift,paste0(names,'down_shift.txt'),sep = '\t')

fold <- 0
R05dox_up_fuzz <- R05dox_fuzz %>%
  dplyr::filter(fuzziness_diff_FDR<=0.05 & fuzziness_log2FC >fold)
R05dox_dwon_fuzz <- R05dox_fuzz %>%
  dplyr::filter(fuzziness_diff_FDR<=0.05 & fuzziness_log2FC < -fold)
fwrite(R05dox_up_fuzz,paste0(names,'up_fuzz.txt'),sep = '\t')
fwrite(R05dox_dwon_fuzz,paste0(names,'down_fuzz.txt'),sep = '\t')


#---- annation
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)
library(ReactomePA)
#setwd('/WORK/lbyybl/WH/rvb/chip-seq/PolII/unique/diff_peak/class')
setwd('/WORK/lbyybl/WH/rvb/cor_uniq/MNase/graph/diff')
annopeak_fuc <- function(peak_fle,name){
  peak <- readPeakFile(peak_fle)
  # pdf(paste0(name,'peak_chr.pdf'),width = 6,height = 5)
  # covplot(peak_fle,chrs = c(paste0('chr',c(1:19,'X','Y'))))
  # dev.off()
  # 
  # pdf(paste0(name,'promoter_hetmap.pdf'),width = 3,height = 7)
  # peakHeatmap(peak_fle, TxDb=txdb, upstream=3000, downstream=3000, color="red")
  # dev.off()
  # 
  # pdf(paste0(name,'promoter_profile.pdf'),width=5,height = 5)
  # plotAvgProf2(peak_fle, TxDb=txdb, upstream=3000, downstream=3000,
  #              xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
  # dev.off()
  
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
  pdf(paste0(name,'upset_anno.pdf'),width = 12,height = 12)
  print(upsetplot(peakAnno))
  dev.off()
  pdf(paste0(name,'unset_venn_anno.pdf'),width = 20,height = 20)
  print(upsetplot(peakAnno, vennpie=TRUE))
  dev.off()
  pdf(paste0(name,'dis2tss_anno.pdf'),width = 6,height = 3)
  print(plotDistToTSS(peakAnno,
                title="Distribution of transcription factor-binding loci\nrelative to TSS"))
  dev.off()
  MF <- enrichGO(as.data.frame(peakAnno)$geneId,OrgDb = "org.Mm.eg.db",keyType = 'ENTREZID',
           ont = 'MF',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.05)
  pdf(paste0(name,'MF_GO.pdf'),width=10,height = 5)
  print(dotplot(MF,showCategory=20))
  dev.off()
  BP <- enrichGO(as.data.frame(peakAnno)$geneId,OrgDb = "org.Mm.eg.db",keyType = 'ENTREZID',
                 ont = 'BP',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.05)
  pdf(paste0(name,'BP_GO.pdf'),width=10,height = 5)
  print(dotplot(BP,showCategory=20))
  dev.off()
  CC <- enrichGO(as.data.frame(peakAnno)$geneId,OrgDb = "org.Mm.eg.db",keyType = 'ENTREZID',
                 ont = 'CC',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.05)
  pdf(paste0(name,'CC_GO.pdf'),width=10,height = 5)
  print(dotplot(CC,showCategory=20))
  dev.off()
  
  pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId,organism = 'mouse')
  head(pathway1, 2)
  gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
  enrich_gene <- data.frame('name'=gene)
  fwrite(enrich_gene,paste0(name,'enrich_gene.txt'),col.names = F)
  pathway2 <- enrichPathway(gene,organism = 'mouse')
  head(pathway2, 2)
  pdf(paste0(name,'seq2gene_enrich.pdf'),width = 10,height = 5)
  print(dotplot(pathway2))
  #print(plot)
  dev.off()
  pdf(paste0(name,'enrich.pdf'),width=10,height = 5)
  print(dotplot(pathway1))
  #print(plot1)
  dev.off()
  
}
annopeak_fuc('R05doxdown_ocup.txt','R05doxdown_ocup')
annopeak_fuc('R05doxup_ocup.txt','R05doxup_ocup')
annopeak_fuc('R1doxdown_ocup.txt','R1doxdown_ocup')
annopeak_fuc('R1doxup_ocup.txt','R1doxup_ocup')

#-----------------------------------------------------------------------
#---火山图
setwd('/WORK/lbyybl/WH/rvb/cor_uniq/MNase/graph/diff/valno')

valno <- function(file,name){
  names <- fread('../name_class.csv')
  
  names <- c(names(names)[1:13],'chr','st','en','gene','score','strand')
  R05_dox <- fread(file,col.names = names)
  # R01_dox <- fread('R1dox.bed',col.names = names)
  R05_dox$logCPM <- log2((R05_dox$control_smt_val + R05_dox$treat_smt_val)/2)
  R05_dox$color <-R05_dox$smt_diff_FDR
  for (i in 1:nrow(R05_dox)){
    if (R05_dox$smt_diff_FDR[i] <= 0.05 & abs(R05_dox$smt_log2FC[i])>1){
      R05_dox$color[i] <- 'red'
    }else{
      R05_dox$color[i] <- 'blue'
    }
  }
  R05_dox <- R05_dox[,c(1:13,20,21)]
  R05_dox <- unique(R05_dox)
}

R05_dox <- valno('R05dox.bed','R05dox')
name <- 'R05dox'
R05_dox_up <- R05_dox %>%
  filter(smt_diff_FDR<0.05,smt_log2FC> 1)
R05_dox_down <- R05_dox %>%
  filter(smt_diff_FDR<0.05,smt_log2FC< -1)
fwrite(R05_dox_up,paste0(name,'up_ocup.txt'),sep = '\t')
fwrite(R05_dox_down,paste0(name,'down_ocup.txt'),sep = '\t')

R1_dox <- valno('R1dox.bed','R1dox')
R1_dox_up <- R1_dox %>%
  filter(smt_diff_FDR<0.05,smt_log2FC>1)
R1_dox_down <- R1_dox %>%
  filter(smt_diff_FDR<0.05,smt_log2FC< -1)

name <- 'R1dox'
fwrite(R1_dox_up,paste0(name,'up_ocup.txt'),sep = '\t')
fwrite(R1_dox_down,paste0(name,'down_ocup.txt'),sep = '\t')

ggplot(R1_dox,aes(x=logCPM,y=smt_log2FC,color=color))+ geom_point(size=0.3,alpha=0.7) +
  geom_hline(aes(yintercept = 0),colour="black",size=0.2,linetype="dashed")+
  #geom_hline(aes(yintercept = log2(2)),color="cyan3",size=0.2,linetype="dashed")+
  #geom_hline(aes(yintercept = log2(0.5)),color="cyan3",size=0.2,linetype="dashed")+
  labs(x="log2 (Degron+untreated)/2",y="log2 Degron/untreated FC",color="") +
  #scale_y_continuous(expand=c(0,0),limits = c(-4,4)) +
  #scale_color_manual(values = c("grey67","blue","blue","red","red"))+
  scale_color_manual(values = c("blue","red"))+
  #scale_color_manual(values = c("grey67","grey67","blue","red"))+
  scale_x_continuous(limits = c(7,10))+
  theme(
    legend.key.width = unit(0.1,"line"),
    legend.key.height = unit(1,"cm"),
    legend.text = element_text(
      size = 10,
      hjust = 0,
      #face = "italic",
      colour = "black",
      angle = 0
    ),
    
    plot.title=element_text(colour="black", size=10),
    axis.line.x = element_blank(),#element_line(color="black", size = 0.5),
    axis.line.y = element_blank(), #element_line(color="black", size = 0.5),
    axis.text.x=element_text(angle=0, hjust=0.5,vjust = 0,
                             colour="black", size=12),
    axis.text.y=element_text(#face = "bold",
      colour="black", size=12),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    panel.background=element_blank()
  )
ggsave(paste0(name,'tss_diff.pdf'),width = 4,height = 2.8)
