library('biomaRt')
library("curl")
library(dplyr)

setwd('/DATA/work/lbyybl/wh/ruvb2/chip-ms')
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
a <- listAttributes(mart)
b <- listFilters(mart)
my_ensembl_gene_id <- read.table('protien_name.txt',header = T,stringsAsFactors = F)
mms_symbols<- getBM(attributes=c('uniprotswissprot','external_gene_name',"mgi_description"),
                    filters = 'uniprotswissprot', values = my_ensembl_gene_id, mart = mart)
head(mms_symbols)
test <- getBM(attributes=c('external_gene_name','name_1006','definition_1006'),
              filters = 'uniprotswissprot', values = my_ensembl_gene_id[1:5,], mart = mart)
fwrite(mms_symbols,'id_trans.tsv',sep = '\t')

