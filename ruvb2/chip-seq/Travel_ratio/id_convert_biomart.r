library('biomaRt')
library("curl")
library(data.table)
library(dplyr)
library("EnsDb.Mmusculus.v79")
options(stringsAsFactors = F)
setwd('/home/boyuanli/bashscript/bin/ruvb2/chip-seq/Travel_ratio')
system('mysql --user=genome -N --host=genome-mysql.cse.ucsc.edu -A -D mm10 -e "select name,name2 from refGene" > Refseq2Gene.txt')
convert_list <- fread('Refseq2Gene.txt',col.names = c('gene_name','symbol'))
ucsc_gene <- '/DATA/work/lbyybl/genomes/mm10/DNA_elements_made_by_Boyuan/mm10_gene.bed'
ucsc_gene <- fread(ucsc_gene,stringsAsFactors = F,header = F)
ucsc_gene <- ucsc_gene[,1:6]
colnames(ucsc_gene) <- c('chr','st','en','gene_name','score','strand')
merge<- merge(ucsc_gene,convert_list,by='gene_name')

length(unique(merge$gene_name))
length(unique(ucsc_gene$gene_name))
length(unique(convert_list$gene_name))
diff <- setdiff(unique(ucsc_gene$gene_name),unique(merge$gene_name))
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
mms_mrna<- getBM(attributes=c('refseq_mrna','external_gene_name'),
                 filters = 'refseq_mrna', values = diff, mart = mart)
colnames(mms_mrna) <- c('gene_name','symbol')
convert_list <- rbind(convert_list,mms_mrna)
merge<- merge(ucsc_gene,convert_list,by='gene_name')

length(unique(merge$gene_name))
length(unique(ucsc_gene$gene_name))
length(unique(convert_list$gene_name))
diff <- setdiff(unique(ucsc_gene$gene_name),unique(merge$gene_name))
convert_list2 <- fread('ncbiRefSeq.txt',col.names = c('gene_name','symbol'))
convert_list2$gene_name <- apply(convert_list2[,1], 1, function(nc){
  strsplit(nc,'\\.')[[1]][[1]]
})
diff_convert <- convert_list2 %>%
  dplyr::filter(gene_name %in% diff)
convert_list <- rbind(convert_list,diff_convert)
merge<- merge(ucsc_gene,convert_list,by='gene_name')

length(unique(merge$gene_name))
length(unique(ucsc_gene$gene_name))
length(unique(convert_list$gene_name))
diff <- setdiff(unique(ucsc_gene$gene_name),unique(merge$gene_name))
diff <- data.frame('gene_name'=diff,
                   'symbol'=c('Tpd52','Cypt10','Cypt9','Cypt7','Cypt8','Celf4'))
convert_list <- rbind(convert_list,diff)
merge<- merge(ucsc_gene,convert_list,by='gene_name')
nrow(unique(merge))

fwrite(unique(merge),'id_convert_gene.bed',sep = '\t',col.names = F)
merge2 <- unique(merge)

#---select unique gene and the longest 要循环至少3 次；
merge <- fread('id_convert_gene.bed',stringsAsFactors = F,
               col.names = c('gene_name','chr','st','en','score','strand','symbol'))
merge <- merge %>%
  arrange(symbol)
merge$length <- abs(merge$en-merge$st)
merge <- unique(merge)
for (i in 1:(nrow(merge)-1)){
  if (merge$symbol[i]==merge$symbol[i+1]){
    if (merge$length[i]==merge$length[i+1]){
      merge[i+1,1:6]<- merge[i,1:6]
    }else if (merge$length[i]>merge$length[i+1]){
      merge[i+1,c(1:6,8)]<- merge[i,c(1:6,8)]
    }else if (merge$length[i]<merge$length[i+1]){
      merge[i,c(1:6,8)]<- merge[i+1,c(1:6,8)]
    }
  }
}
nrow(unique(merge))
final_merge <- unique(merge)
merge <- merge[,c(2:4,7,5,6,1,8)]
merge <- merge %>%
  arrange(chr,st,en)
fwrite(unique(merge),'long_unique_gene.bed',sep = '\t',col.names = F)
#fwrite(diff,'diff.txt')
# # listMarts()
# # usemart <- useMart('ENSEMBL_MART_MOUSE')
# # database <- listDatasets(usemart)
# # mart <- useDataset('maj_gene_ensembl',usemart)
# mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
# 
# a <- listAttributes(mart)
# b <- listFilters(mart)
# # my_ensembl_gene_id <- read.table('protien_name.txt',header = T,stringsAsFactors = F)
# ucsc_mrna <- ucsc_gene %>%
#   filter('NM' == substr(gene_name,1,2))
# ucsc_ncrna <- ucsc_gene %>%
#   filter('NR' == substr(gene_name,1,2))
# mms_mrna<- getBM(attributes=c('refseq_mrna','external_gene_name'),
#                     filters = 'refseq_mrna', values = unique(ucsc_mrna$gene_name), mart = mart)
# mms_ncrna<- getBM(attributes=c('refseq_ncrna','external_gene_name'),
#                     filters = 'refseq_ncrna', values = unique(ucsc_ncrna$gene_name), mart = mart)
# head(mms_symbols)
# test <- getBM(attributes=c('external_gene_name','name_1006','definition_1006'),
#       filters = 'uniprotswissprot', values = my_ensembl_gene_id[1:5,], mart = mart)
# fwrite(mms_symbols,'id_trans.tsv',sep = '\t')
# 
# #BiocManager::install("EnsDb.Mmusculus.v79")
# 
# keytypes(EnsDb.Mmusculus.v79)
# select(EnsDb.Mmusculus.v79, key=unique(ucsc_mrna$gene_name)[1:30], 
#                columns=c("SEQNAME", "SYMBOL"), 
#                keytype="UNIPROTID")
