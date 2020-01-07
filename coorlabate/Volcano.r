rm(list = ls())
setwd('/DATA2/work/lbyybl/coorlaborate/YB/YB_RNA_seq/total/output/graph/38nodaoxvs38jiadox/contact_valo')
dir()
library(dplyr)
library(data.table)
gene_test <- read.table('contact_gene.tsv',stringsAsFactors = F, header = F)
colnames(gene_test) <- c('chr','st','en','gene','baseMean','log2FoldChange',
                         'lfcSE','stat','pvalue','padj')
head(gene_test)
id2gene <- read.table('/DATA/work/lbyybl/genomes/hg19/id_convert_gene.bed')#,col.names = c('gene','symbol'))
id2gene <- id2gene[,c(1,7)]
colnames(id2gene) <- c('gene','symbol')
gene_test <- merge(gene_test,id2gene,by='gene')
gene_test$sig <- 'no'
gene_test$change <- 'no'
gene_test[is.na(gene_test)] <- 1
for (i in 1:nrow(gene_test)){
  if (gene_test$padj[i] < 0.05){
    gene_test$sig[i] <- 'Yes'
  }
}

for (i in 1:nrow(gene_test)){
  if (gene_test$log2FoldChange[i] > 1){
    gene_test$change[i] <- 'up'
  } else if (gene_test$log2FoldChange[i] < -1) {
    gene_test$change[i] <- 'down'
  }
}

gene_test$ch_si <- paste0(gene_test$sig, gene_test$change)
for (i in 1:nrow(gene_test)){
  if (gene_test$sig[i] =='no'){
    gene_test$ch_si[i] <-  'no'
  }
}

for (i in 1:nrow(gene_test)){
  if (gene_test$change[i] =='no'){
    gene_test$ch_si[i] <-  'no'
  }
}

library(ggplot2)
label_file <- read.table('label.txt')
colnames(label_file) <- c('gene','chr','st','en','baseMean','log2FoldChange',
                           'lfcSE','stat','pvalue','padj','symbol','sig','change','ch_si')
ggplot(gene_test,aes(x=log2FoldChange,y=-log10(padj),color=ch_si)) + geom_point(size=2) +
  labs(x="log2 Fold Change",y="-log (adjusted p-value)",color="")+
  geom_vline(xintercept = c(-log(2),log(2)),linetype='dashed',
             color='darkslateblue',size=1)+
  geom_hline(yintercept = -log10(0.05),linetype='dashed',
             color='darkslateblue',size=1)+
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(0,20))+
  scale_color_manual(values=c("grey", "#005792",  "#ca3e47"), name="", 
                     label=c("unchanged (31750)", "down (358)", "up (300)"))+
  theme(plot.title=element_text(colour="black", size=14, hjust=0.5),
        axis.line.x = element_blank(), axis.line.y = element_blank(),
        legend.text = element_text(size=15),
        axis.text=element_text(colour="black", size=17), axis.title = element_text(size=20),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.background=element_blank())+
  geom_label(data=label_file,aes(x=log2FoldChange,y=-log10(padj),
                                 label=symbol),color='black',size=2)

ggsave('contact_gene.pdf',width = 7, height = 5)


change_gene <- gene_test %>%
  filter(ch_si!='no')
up_gene <- change_gene %>% 
  filter(ch_si=='Yesup')
down_gene <- change_gene %>%
  filter(ch_si=='Yesdown')

fwrite(unique(up_gene),'contact_up.tsv',sep = '\t', col.names = F)
fwrite(unique(down_gene),'contact_down.tsv',sep = '\t', col.names = F)
