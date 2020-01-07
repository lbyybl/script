#--- select all the transcription factor, co factor, chromatin factor, and disease
# associate protein
library(ggplot2)
library(dplyr)
library(data.table)
library(Hmisc)
# select the protein exit in JX2015 and find the percentage of Ruvbl2 related complex protein in these protein with this chip-ms
setwd('/DATA/work/lbyybl/wh/ruvb2/chip-ms')
select <- fread('sub_merge.txt',header = T,stringsAsFactors = F)
names(select)[14] <- 'MW'
select <- select %>%
  dplyr::filter(rep1+rep2+rep3>=1)
JXprotein <- fread('JX2015.txt',header = F,stringsAsFactors = F)
names(JXprotein) <- 'gene_name'
merge <- merge(select,JXprotein,by='gene_name')
nrow(merge)
merge[is.na(merge)] <- ''
merge <- merge %>%
  dplyr::mutate(qantity=(rep1+rep2+rep3)/MW)
merge$stat <- 'notin'
head(merge)
srcap<-c("A0A087WQ44","P60122","Q9WTM5","Q62481","Q9Z2N8","Q9CR11","P60710","Q9JI44","Q9D864","Q8R331")
tip <- c("Q8CHK4","P60122","Q9WTM5","Q62481","Q9Z2N8","Q9CR11","P60710","Q9JI44","Q80YV3","P60762","Q2VPQ9","Q9DAT2","Q8C9X6","Q8VEK6","Q8R3B7","Q8CHI8")
ino<-c("Q6ZPV2","P60122","Q9WTM5","Q9Z2N8","Q80US4","Q8R2S9","Q8BHA0","Q99PT3","Q66JY2","A0A0U1RP99","Q99L90","Q6PIJ4","Q3U1J1","Q9WUP7","Q00899")
Ruvbl_chrom <- c(srcap,tip,ino)
for (i in 1:nrow(merge)){
  if (merge$Accession[i] %in% Ruvbl_chrom){
    merge$stat[i] <- 'in'
  }
}
merge %>%
  group_by(stat) %>%
  summarise(n=n())
explain_per<-sum(merge[merge$stat=="in",]$qantity)/sum(merge$qantity)
sum(merge[merge$stat=="notin",]$qantity)
fwrite(merge,'Ruvbl_chrom_detect.csv')
#--------------------------------------------------------
setwd('/DATA/work/lbyybl/wh/ruvb2/chip-ms/class')
TF <- fread('Mus_musculus_TF',header = T, stringsAsFactors = F)
TF <- TF %>%
  dplyr::select(Symbol)
protein_list <- fread('../sub_merge.txt',header = T,stringsAsFactors = F)
protein_list <- protein_list %>%
  dplyr::filter(rep1+rep2+rep3>=1)
protein_TF <- merge(protein_list,TF,by.x='gene_name',by.y='Symbol')
fwrite(protein_TF,'TF_detected.csv')

protein_TF <- fread('TF_detected.txt',header = F,stringsAsFactors = F)
TF_homer <- fread('../knownResults.txt',header = F,stringsAsFactors = F)
merge_TF <- merge(protein_TF,TF_homer,by='V1')
fwrite(merge_TF,'merge_TF.txt')
#----------------------------------------------------------------
co_factor <- fread('Mus_musculus_TF_cofactors',header = T, stringsAsFactors = F)
co_factor <- co_factor %>%
  dplyr::select(Symbol)
protein_cofactor <- merge(protein_list,co_factor,by.x='gene_name',by.y='Symbol')
fwrite(protein_cofactor,'cofactor_detected.csv')

# remolder <- fread('Mus_musculus_TF_cofactors',header = T, stringsAsFactors = F)
# remolder <- remolder %>%
#   dplyr::select(Symbol)
# protein_remolder <- merge(protein_list,remolder,by.x='gene_name',by.y='Symbol')
# fwrite(protein_remolder,'cofactor_detected.txt')

#--- overlap with disease
# library(devtools)
# install_bitbucket("ibi_group/disgenet2r")
library(disgenet2r)
# omim <- fread('disease/mim2gene.txt',stringsAsFactors = F,header = F)
# omim <- omim %>%
#   dplyr::select(V4)
# #names(omim) <- 'Symbol'
# #omim[which(is.na(omim))]
# omim <- omim %>%
#   arrange(V4)
# nrow(omim)
# omim <- omim %>%
#   filter(V4!="")

# dis1 <- fread('disease/all_gene_disease_pmid_associations.tsv',
#               stringsAsFactors = F,header = T)
# dis1 <- dis1 %>%
#   dplyr::select(geneSymbol)
# head(dis1)
dis2 <- fread('disease/curated_gene_disease_associations.tsv',
              stringsAsFactors = F,header = T)
# dis2 <- dis2 %>%
#   dplyr::select(geneSymbol)
dis2$geneSymbol <- tolower(dis2$geneSymbol)
dis2$geneSymbol <- capitalize(dis2$geneSymbol)
head(dis2)

# cosmic <- fread('/DATA/work/lbyybl/ZR/cosmic/CosmicMutantExport.tsv',
#                 stringsAsFactors = F,header = T)
# names(cosmic)[1] <- 'Symbol'
# cosmic <- cosmic %>%
#   dplyr::select(Symbol)
# fwrite(cosmic,'cosmic_gene.txt')
# cosmic <- fread('cosmic_gene.txt',
#                 stringsAsFactors = F,header = T)
# # for (i in 1:nrow(cosmic)){
# #   cosmic$Symbol[i] <- strsplit(cosmic$Symbol[i],"_")[[1]][1]
# # }
# # cosmic$Symbol[6721293]
# # strsplit(cosmic$Symbol[1:100],"_")[[1]][1]
# names(dis1) <- 'gene'
# names(dis2) <- 'gene'
# names(omim) <- 'gene'
# names(cosmic) <- 'gene'
# disease <- rbind(dis1,dis2)
# disease <- rbind(disease,omim)
# disease <- rbind(disease,cosmic)
# disease2 <- data.frame('gene'=unique(disease$gene),stringsAsFactors = F)
# head(disease2)
# nrow(disease2)
# protein_list$gene_name <- toupper(protein_list$gene_name)
# disease2$gene <- toupper(disease2$gene)
protein_disease <- merge(protein_list,dis2,by.x='gene_name',by.y='geneSymbol')
protein_disease <- unique(protein_disease)
fwrite(protein_disease,'disease_detected.csv')
length(unique(protein_disease$gene_name))
length(unique(protein_disease$diseaseName))
nrow(protein_disease)
#-----------------------------------------
library(Hmisc)
cell_cycle <- fread('cell_cycle.txt',stringsAsFactors = F,header = F)
cell_cycle$V1 <- tolower(cell_cycle$V1)
cell_cycle$V1 <- capitalize(cell_cycle$V1)
protein_cycle <- merge(protein_list,cell_cycle,by.x='gene_name',by.y='V1')
fwrite(protein_cycle,'cellcycle_detected.csv')
#--------------------------------------------
# uniprot cycle cell and singal transduction
cell_cycle <- fread('uniprot_class/uniprot-cell_cycle.tab',stringsAsFactors = F,header = T)
cell_cycle <- cell_cycle %>%
  dplyr::select(Entry)
protein_cycle <- merge(protein_list,cell_cycle,by.x='Accession',by.y='Entry')
fwrite(protein_cycle,'uniprot_cellcycle_detected.csv')
#----- singal transduction
signal_trans <- fread('uniprot_class/uniprot-signal_transduction.tab',stringsAsFactors = F,header = T)
signal_trans <- signal_trans %>%
  dplyr::select(Entry)
protein_signal <- merge(protein_list,signal_trans,by.x='Accession',by.y='Entry')
fwrite(protein_signal,'uniprot_signal_detected.csv')
#----- chromatin
chromatin <- fread('uniprot_class/uniprot-chromosome.tab',stringsAsFactors = F,header = T)
chromatin <- chromatin %>%
  dplyr::select(Entry)
protein_chrom <- merge(protein_list,chromatin,by.x='Accession',by.y='Entry')
fwrite(protein_chrom,'uniprot_chrom_detected.csv')
names(protein_chrom)[14] <- 'MW'

nrow(protein_chrom)
protein_chrom[is.na(protein_chrom)] <- ''
protein_chrom <- protein_chrom %>%
  dplyr::mutate(qantity=(rep1+rep2+rep3)/MW)
protein_chrom$stat <- 'notin'
head(protein_chrom)
srcap<-c("A0A087WQ44","P60122","Q9WTM5","Q62481","Q9Z2N8","Q9CR11","P60710","Q9JI44","Q9D864","Q8R331")
tip <- c("Q8CHK4","P60122","Q9WTM5","Q62481","Q9Z2N8","Q9CR11","P60710","Q9JI44","Q80YV3","P60762","Q2VPQ9","Q9DAT2","Q8C9X6","Q8VEK6","Q8R3B7","Q8CHI8")
ino<-c("Q6ZPV2","P60122","Q9WTM5","Q9Z2N8","Q80US4","Q8R2S9","Q8BHA0","Q99PT3","Q66JY2","A0A0U1RP99","Q99L90","Q6PIJ4","Q3U1J1","Q9WUP7","Q00899")
Ruvbl_chrom <- c(srcap,tip,ino)
for (i in 1:nrow(protein_chrom)){
  if (protein_chrom$Accession[i] %in% Ruvbl_chrom){
    protein_chrom$stat[i] <- 'in'
  }
}
protein_chrom %>%
  group_by(stat) %>%
  summarise(n=n())
explain_per<-sum(protein_chrom[protein_chrom$stat=="in",]$qantity)/sum(protein_chrom$qantity)
sum(protein_chrom[protein_chrom$stat=="notin",]$qantity)
#---- cofactor
cofactor <- fread('uniprot_class/uniprot-cofactor.tab',stringsAsFactors = F,header = T)
cofactor <- cofactor %>%
  dplyr::select(Entry)
protein_cofactor <- merge(protein_list,cofactor,by.x='Accession',by.y='Entry')
fwrite(protein_cofactor,'uniprot_cofactor_detected.csv')
#---- transcription factor
transfactor <- fread('uniprot_class/uniprot-Transcription_factors.tab',stringsAsFactors = F,header = T)
transfactor <- transfactor %>%
  dplyr::select(Entry)
protein_TF <- merge(protein_list,transfactor,by.x='Accession',by.y='Entry')
fwrite(protein_TF,'uniprot_TF_detected.csv')
#--------------------------------------------
# draw graph
class <- data.frame('class'=c("TF",'cofactor','cell_cycle','chromatin','signal'),
                    'value'=c(101,369,655,317,270))
enrich_class <- data.frame('class'=c("TF",'cofactor','cell_cycle','chromatin','signal'),
                    'value'=c(14,62,61,47,26))
class <- class %>%
  arrange(-value)
enrich_class <- enrich_class %>%
  arrange(-value)
ggplot(enrich_class,aes(class,value)) + geom_col(fill='lightblue',width = 0.5) +
  geom_text(label=enrich_class$value)+
  scale_x_discrete(limits=enrich_class$class)+
  theme(
    #title = element_blank() ,
    legend.key.width = unit(2,"line"),
    legend.key.height = unit(0.5,"cm"),
    legend.text = element_text(
      size = 12,
      hjust = 1.2,
      face = "italic",
      colour = "black"
    ),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    plot.title=element_text(colour="black", size=10),
    axis.text.x=element_text(angle = -30,hjust = -0.05,
                             colour="black", size=7),
    axis.text.y=element_text(#face = "bold",
      colour="black", size=12),
    panel.background=element_blank())
ggsave('protein_class2_enrich.pdf',width = 5,height = )
#plot(class$class,class$value)
