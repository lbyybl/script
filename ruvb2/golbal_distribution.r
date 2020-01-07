#--- sts the distribution of ruvb2
library(data.table)
library(dplyr)
library(ggplot2)
require("RColorBrewer")
library(plotrix)
#setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/bw_bam/Ruvbl2/merge/graph/distribution/pie_plot')
#setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/Ruvbl1/sample20190509/graph/subtract_heatmap/distribution')
setwd('/WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl/R1_2_overlap/grahp/pie')
distribu_file <- fread('all_sts.tsv',col.names = c('chr','st','en','strand','element'))
distribu_file$class <- 'no'
for (i in 1:nrow(distribu_file)){
  if (stringr::str_detect(distribu_file$element[i],'promo')){
    distribu_file$class[i] <- 'Promoter'
  }else if (stringr::str_detect(distribu_file$element[i],'Super_enhancer')){
    distribu_file$class[i] <- 'Super_enhancer'
  }else if (stringr::str_detect(distribu_file$element[i],'Enhancer')){
    distribu_file$class[i] <- 'Enhancer'
  }else if (stringr::str_detect(distribu_file$element[i],'Insulator')){
    distribu_file$class[i] <- 'Insulator'
  }else if (stringr::str_detect(distribu_file$element[i],'genebody')){
    distribu_file$class[i] <- 'genebody'
  }else if (stringr::str_detect(distribu_file$element[i],'None_gene')){
    distribu_file$class[i] <- 'None_gene'
  }
}
# distribu_file %>%
#   filter(class=='no')
# distribu_file <- distribu_file %>%
#   mutate(V10=substring(V8,2,4))
sts <- distribu_file %>%
  filter(class!='no') %>%
  group_by(class) %>%
  summarise(value=n(),length=sum(en-st)) 
sts$per <- paste0(round(sts$value/sum(sts$value)*100,2),"%")
sts$enrich <- sts$value/sts$length
sts <- sts %>%
  arrange(-value) 
pdf('Ruvbl1_2overlap_ele.pdf')
pie(sts$value/sum(sts$value),labels = sts$per,radius=0.7,
    cex=0.8,col=sort(brewer.pal(6,"Set1")),
    border="white")
legend("bottomright",legend=sts$class,cex=0.6,bty="n",
       fill=sort(brewer.pal(6,"Set1")))
dev.off()
pdf('Ruvbl1_2overlap_ele3d.pdf')
pie3D(sts$value/sum(sts$value),labels = sts$per,radius=0.7,
      labelcex=0.8,col=sort(brewer.pal(6,"Set1")),
      border="gray",theta = 1.2)
dev.off()


#--------------------------------- promoter distribution
promoter_dis <- fread('four_gene_prmototer_sts.tsv',
                      col.names = c('chr','st','en','chr2','st2','en2','strand','element','ov_len'))
promoter_dis$class <- 'no'
for (i in 1:nrow(promoter_dis)){
  if (stringr::str_detect(promoter_dis$element[i],'Activ')){
    promoter_dis$class[i] <- 'Active'
  }else if (stringr::str_detect(promoter_dis$element[i],'Initiat')){
    promoter_dis$class[i] <- 'Initiated'
  }else if (stringr::str_detect(promoter_dis$element[i],'Bivalent')){
    promoter_dis$class[i] <- 'Bivalent'
  }else if (stringr::str_detect(promoter_dis$element[i],'Silent')){
    promoter_dis$class[i] <- 'Silent'
  }else if (stringr::str_detect(promoter_dis$element[i],'Unclassified')){
    promoter_dis$class[i] <- 'Others'
  }
}
# distribu_file %>%
#   filter(class=='no')
# distribu_file <- distribu_file %>%
#   mutate(V10=substring(V8,2,4))
psts <- promoter_dis %>%
  filter(class!='no') %>%
  group_by(class) %>%
  summarise(value=n(),length=sum(en-st)) 
psts$per <- paste0(round(psts$value/sum(psts$value)*100,2),"%")
psts$enrich <- psts$value/psts$length
psts <- psts %>%
  arrange(-value) 
pdf('Ruvbl1_2overlap_pro.pdf')
pie(psts$value/sum(psts$value),labels = psts$per,radius=0.7,
    cex=0.8,col=sort(brewer.pal(6,"Set1")),
    border="white")
legend("bottomright",legend=psts$class,cex=0.6,bty="n",
       fill=sort(brewer.pal(6,"Set1")))
dev.off()
pdf('Ruvbl1_2overlap_pro3d.pdf')
pie3D(psts$value/sum(psts$value),labels = psts$per,radius=0.7,
      labelcex=0.8,col=sort(brewer.pal(6,"Set1")),
      border="gray",theta = 1.2)
dev.off()
# blank_theme <- theme_minimal()+
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     panel.border = element_blank(),
#     panel.grid=element_blank(),
#     axis.ticks = element_blank(),
#     plot.title=element_text(size=14, face="bold")
#   )
# ggplot(arrange(sts, class),aes(x='',y=value,fill=class)) + geom_col(width = 1,stat = "identity") + coord_polar("y", start=0) +
#   blank_theme + scale_fill_brewer(palette="Dark2") +geom_text(aes(label=sts$per))+
#   theme(axis.text.x=element_blank())+
#   scale_fill_discrete(limits=sts$class)
# ggsave('elements_distribute.pdf',width = 8,height = 5)
# ggsave('four_promoter_rick_distribute.pdf',width = 8, height = 5)

# ggplot(sts,aes(x='',y=enrich,fill=V10)) + geom_col(width = 1) + coord_polar("y", start=0) +
#   blank_theme + scale_fill_brewer(palette="Dark2") +
#   theme(axis.text.x=element_blank())
# ggsave('elements_enrich.pdf',width = 8, height = 5)
# ggsave('promoter_enrich.pdf',width = 8, height = 5)