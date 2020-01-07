library('biomaRt')
library("curl")
library(dplyr)
# setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/Ruvbl1/sample20190509/graph/subtract_heatmap/Ruvb1_2_overlap_motif/graph')
# mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
# my_ensembl_gene_id <- read.table('chrM_gene.txt')
# mms_symbols<- getBM(attributes=c('ensembl_gene_id','entrezgene_id','external_gene_name','goslim_goa_description',"description"),
#                     filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
# head(mms_symbols)
# fwrite(mms_symbols,'id_trans.tsv',sep = '\t')

setwd('/DATA/work/lbyybl/wh/ruvb2/chip-ms')
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
a <- listAttributes(mart)
b <- listFilters(mart)
my_ensembl_gene_id <- read.table('protien_name.txt',header = T,stringsAsFactors = F)
mms_symbols<- getBM(attributes=c('uniprotswissprot','external_gene_name',"mgi_description",'go_id'),
                    filters = 'uniprotswissprot', values = my_ensembl_gene_id, mart = mart)
head(mms_symbols)
test <- getBM(attributes=c('external_gene_name','name_1006','definition_1006'),
      filters = 'uniprotswissprot', values = my_ensembl_gene_id[1:5,], mart = mart)
fwrite(mms_symbols,'id_trans.tsv',sep = '\t')

# chipms <- read.csv('WH_20190620_R1_R2_ChIP-MS.csv',header = T,stringsAsFactors = F)
# 
# merge <- merge(chipms,mms_symbols,by.x='Accession',by.y='uniprotswissprot',all.x=T)
# names(merge)[26:30] <- c('rep1','rep2','rep3','gfp','input')
# merge[is.na(merge)] <- 0
# merge <- merge %>%
#   mutate(ipvsgfp=(rep1+rep2+rep3+1)/((gfp+1)*3),ipvsinput=(rep1+rep2+rep3+1)/((input+1)*3)) %>%
#   mutate(enrich=ipvsgfp+ipvsinput)
# merge <- merge %>%
#   arrange(-enrich,-(rep1+rep2+rep3))
# fwrite(merge,'chipms_trans_id.csv',sep = ',')
# 
# sub_merge <- merge
# sub_merge <- merge %>%
# #  filter(rep1+rep2+rep3>10) %>%
#   mutate(sum=(rep1+rep2+rep3))
# sub_merge$class <- 'no'
# for (i in 1:nrow(sub_merge)){
#   if (sub_merge$Accession[i] %in% c("Q6ZPV2","P60122","Q9WTM5","Q80US4","Q8R2S9","Q9Z2N8","Q99PT3","Q8BHA0","Q3U1J1","Q6PIJ4","Q99L90")){
#     sub_merge$class[i] <- 'ino'
#   }else if(sub_merge$Accession[i] %in% c("Q8CHK4","P60122","Q9WTM5","Q62481","Q9Z2N8","Q9CR11","P60710","Q9JI44","Q80YV3","P60762","Q2VPQ9","Q9DAT2","Q8C9X6","Q8VEK6","Q8R3B7","Q8CHI8")){
#     sub_merge$class[i] <- 'tip'
#   }else if(sub_merge$Accession[i] %in% c("SRCAP","P60122","Q9WTM5","Q62481","Q9Z2N8","Q9CR11","P60710","Q9JI44","Q9D864")){
#     sub_merge$class[i] <- 'srcap'
#   } else if (sub_merge$Accession[i] %in% c("Q9CQJ2","P60122","Q9WTM5","Q9D706")){
#     sub_merge$class[i] <- 'r2tp'
#   }
# }
# sub_merge$sz <- as.numeric(0.1)
# for (i in 1:nrow(sub_merge)){
#   if (sub_merge$class[i]!='no'){
#     sub_merge$sz[i] <- as.numeric(0.5)
#   }
# }
# 
# sub_merge$num <- 1:nrow(sub_merge)
# ggplot(sub_merge,aes(x=num,y=enrich,color=class,size=sz))+geom_point()+
#   scale_x_continuous(limits = c(1,1500))+
#   scale_y_continuous(limits = c(1,3))
# ggsave('distribu3.pdf',width = 6,height = 5)
