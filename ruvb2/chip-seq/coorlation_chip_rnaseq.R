#--- to draw the coorlation between chip-seq and RNA-seq
# 想看的是Ruvbl1/2 binding与基因表达的关系，所以是将Ruvbl1/2 分成表达量高中低
# 看基因表达怎么变
library(data.table)
library(dplyr)
library(ggplot2)
setwd('/WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl/R1_2_overlap/grahp/gene_coorlation')

exon <- fread('exon.bed',col.names = c('chr','st','en','gene','nu','strand'))
exon$len <- abs(exon$en-exon$st)
exon <- exon %>%
  dplyr::select(gene,len) %>%
  group_by(gene) %>%
  summarise(length=sum(len))
length(unique(exon$gene))
nrow(exon)

over_gene <- fread('overlap_gene.bed',col.names = c('chr','st','en','gene'))
over_gene <- over_gene %>%
  select(gene)
# RNA_seq <- fread('RNA-seq_count.txt',stringsAsFactors = F,
#                  col.names = c('gene','sympol','stat','rep1','rep2'))
# RNA_seq$rna <- RNA_seq$rep1 + RNA_seq$rep2
# RNA_seq <- RNA_seq %>%
#   select(gene,rna)
RNA_seq <- fread('Rna_count.tsv',col.names = c('gene','rna'),stringsAsFactors = F)
Ruvbl1 <- fread('R1_1ChIP_FKDL190735073_merge.tsv',
                col.names = c('gene','Ruvbl1'),stringsAsFactors = F)
Ruvbl2 <- fread('Ruvb2_34_merge.tsv',col.names = c('gene','Ruvbl2'),stringsAsFactors = F)
H3k9me3 <- fread('H3K9me3_unique_uniqe.tsv',col.names = c('gene','H3k9me3'),stringsAsFactors = F)
H3k27ac <- fread('H3k27ac_rmdup_uniqe.tsv',col.names = c('gene','H3k27ac'),stringsAsFactors = F)
join <- dplyr::full_join(RNA_seq,Ruvbl1,by='gene')
join <- dplyr::full_join(join,Ruvbl2,by='gene')
join <- dplyr::full_join(join,H3k9me3,by='gene')
join <- dplyr::full_join(join,H3k27ac,by='gene')
join <- dplyr::full_join(join,exon,by='gene')
join <- dplyr::left_join(over_gene,join,by='gene')
nrow(join)
head(join)
join[is.na(join)] <- 0
sum(join$Ruvbl2>10)
join <- join %>%
  filter(Ruvbl2>10 & Ruvbl1 > 10)

RPKM <- function(cuout,length,total_count){
  
  rpkm.count <- cuout*1000000000/(as.numeric(length)*total_count)
  return((rpkm.count))
}
RPKM2 <- function(cuout,length,total_count){
  
  rpkm.count <- cuout*1000000000/(4000*total_count)
  return((rpkm.count))
}
sum_rna <- sum(join$rna)
sum_R1 <- sum(join$Ruvbl1)
sum_R2 <- sum(join$Ruvbl2)
sum_h3k9 <- sum(join$H3k9me3)
sum_h3k27 <- sum(join$H3k27ac)
join.rpkm<-join
# notice !!!! this is log transformed rpkm
for (i in 1:nrow(join)){
  join.rpkm$rna[i] <- RPKM(join$rna[i],join$length[i],sum_rna)
  join.rpkm$Ruvbl1[i] <- RPKM2(join$Ruvbl1[i],join$length[i],sum_R1)
  join.rpkm$Ruvbl2[i] <- RPKM2(join$Ruvbl2[i],join$length[i],sum_R2)
  join.rpkm$H3k9me3[i] <- RPKM2(join$H3k9me3[i],join$length[i],sum_h3k9)
  join.rpkm$H3k27ac[i] <- RPKM2(join$H3k27ac[i],join$length[i],sum_h3k27)
}
join.rpkm <- join.rpkm %>%
  filter(length>0 & Ruvbl2>0 )
# filter_3p <- function(data,i){
#   data <- data[order(data[i]),]
#   n=nrow(data)
#   p=3*n/100
#   data <- data[p:(n-p),]
#   return(data)
# }
# join <- join %>%
#   filter(rna>0 & Ruvbl1 > 0 & Ruvbl2 > 0 & H3k9me3 >0 )
# for (i in 2:ncol(join.rpkm)){
#   join.rpkm <- filter_3p(join.rpkm,i)
# }

# join <- join %>%
#   filter(rna>0)

# join.rpkm <- join.rpkm %>%
#   arrange(rna)
# join.rpkm$rank <- 1:nrow(join.rpkm)

n=nrow(join.rpkm)
f=round(n/3)
join.rpkm <- join.rpkm %>%
  arrange(Ruvbl1)

mean(join.rpkm$rna[1:f])
mean(join.rpkm$rna[f:(2*f)])
mean(join.rpkm$rna[(2*f):n])

mean(join.rpkm$Ruvbl2[1:f])
mean(join.rpkm$Ruvbl2[f:(2*f)])
mean(join.rpkm$Ruvbl2[(2*f):n])

join.rpkm$R1_rank <- 1:nrow(join.rpkm)

join.rpkm$R1_rank[1:f] <- 'low'
join.rpkm$R1_rank[f:(2*f)] <- 'media'
join.rpkm$R1_rank[(2*f):n] <- 'high'

join.rpkm <- join.rpkm %>%
  arrange(Ruvbl2)
join.rpkm$R2_rank <- 1:nrow(join.rpkm)

join.rpkm$R2_rank[1:f] <- 'low'
join.rpkm$R2_rank[f:(2*f)] <- 'media'
join.rpkm$R2_rank[(2*f):n] <- 'high'

join.rpkm <- join.rpkm %>%
  arrange(H3k27ac)
join.rpkm$k27_rank <- 1:nrow(join.rpkm)

join.rpkm$k27_rank[1:f] <- 'low'
join.rpkm$k27_rank[f:(2*f)] <- 'media'
join.rpkm$k27_rank[(2*f):n] <- 'high'

join.rpkm <- join.rpkm %>%
  arrange(H3k9me3)
join.rpkm$k9_rank <- 1:nrow(join.rpkm)
join.rpkm$k9_rank[1:f] <- 'low'
join.rpkm$k9_rank[f:(2*f)] <- 'media'
join.rpkm$k9_rank[(2*f):n] <- 'high'

join.rpkm2 <- join.rpkm %>%
  mutate(Ruvbl1='Ruvbl1',Ruvbl2='Ruvbl2',H3k27ac='H3k27ac',H3k9me3='H3k9me3')

join.rpkm2 <- join.rpkm2[,c(2,8:11)]
library(reshape2)
join.rpkm2 <- melt(join.rpkm2,id.vars = 'rna')

ggplot(join.rpkm2,aes(x=value,y=rna,fill=variable))+ geom_boxplot(outlier.color = NA)+
  scale_y_log10() +
  scale_x_discrete(limits=c('low','media','high'))+
  theme(plot.title=element_text(size=40,hjust=0.5),
        panel.border = element_rect(colour="black",  fill=NA, size=1),
        axis.title=element_text(size=10),
        axis.text=element_text(size=10),
        legend.title=element_text(size=25),
        panel.background=element_blank(),
        legend.text=element_text(size=22))+
  ylab("RNA expression")+xlab("ChIP-seq binding desinty")

ggsave('Chip-seq-rna-cor-boxplot.pdf')


#-----------------------使用分位数

n=nrow(join.rpkm)
f=round(n/3)
join.rpkm <- join.rpkm %>%
  arrange(Ruvbl1)

mean(join.rpkm$rna[1:f])
mean(join.rpkm$rna[f:(2*f)])
mean(join.rpkm$rna[(2*f):n])

sum(join.rpkm$Ruvbl1 < f)
mean(join.rpkm$Ruvbl2[f:(2*f)])
mean(join.rpkm$Ruvbl2[(2*f):n])

join.rpkm$R1_rank <- 1:nrow(join.rpkm)
range(join.rpkm$Ruvbl1)
f=(31.77048-0)/3
join.rpkm$R1_rank <- 'high'
join.rpkm$R1_rank[join.rpkm$Ruvbl1 < (2*f)] <- 'media'
join.rpkm$R1_rank[join.rpkm$Ruvbl1 < f] <- 'low'



join.rpkm <- join.rpkm %>%
  arrange(Ruvbl2)
join.rpkm$R2_rank <- 1:nrow(join.rpkm)

join.rpkm$R2_rank[1:f] <- 'low'
join.rpkm$R2_rank[f:(2*f)] <- 'media'
join.rpkm$R2_rank[(2*f):n] <- 'high'

join.rpkm <- join.rpkm %>%
  arrange(H3k27ac)
join.rpkm$k27_rank <- 1:nrow(join.rpkm)

join.rpkm$k27_rank[1:f] <- 'low'
join.rpkm$k27_rank[f:(2*f)] <- 'media'
join.rpkm$k27_rank[(2*f):n] <- 'high'

join.rpkm <- join.rpkm %>%
  arrange(H3k9me3)
join.rpkm$k9_rank <- 1:nrow(join.rpkm)
join.rpkm$k9_rank[1:f] <- 'low'
join.rpkm$k9_rank[f:(2*f)] <- 'media'
join.rpkm$k9_rank[(2*f):n] <- 'high'

join.rpkm2 <- join.rpkm %>%
  mutate(Ruvbl1='Ruvbl1',Ruvbl2='Ruvbl2',H3k27ac='H3k27ac',H3k9me3='H3k9me3')

join.rpkm2 <- join.rpkm2[,c(2,8:11)]
library(reshape2)
join.rpkm2 <- melt(join.rpkm2,id.vars = 'rna')

pdf('scater_cor.pdf')
plot(join.rpkm2)
dev.off()
cor(join.rpkm$rna,join.rpkm$Ruvbl1)
cor(join.rpkm$rna,join.rpkm$Ruvbl2)
cor(join.rpkm$rna,join.rpkm$H3k9me3)
plot(join.rpkm$Ruvbl1,join.rpkm$Ruvbl2)
cor(join.rpkm$rank,join.rpkm$R1_rank)
cor(join.rpkm$R1_rank,join.rpkm$R2_rank)
plot(join.rpkm$R1_rank,join.rpkm$R2_rank)
ggplot(join.rpkm,aes(rank,Ruvbl1))+ geom_smooth()+
  geom_smooth(aes(rank,Ruvbl2),color='red')+
  geom_smooth(aes(rank,H3k27ac),color='yellow')+
  geom_smooth(aes(rank,H3k9me3),color='green')  +
  theme(plot.title=element_text(size=40,hjust=0.5),
        panel.border = element_rect(colour="black",  fill=NA, size=1),
        axis.title=element_text(size=10),
        axis.text=element_text(size=10),
        legend.title=element_text(size=25),
        panel.background=element_blank(),
        legend.text=element_text(size=22))+
  ylab("ChIP-seq density")+xlab("Genes sorted by expression")
ggsave("chip-rna_cor1.pdf")
cor(join.rpkm$rna,join.rpkm$Ruvbl1)
cor(join.rpkm$rna,join.rpkm$Ruvbl2)
cor(join.rpkm$rna,join.rpkm$H3k9me3)
cor(join.rpkm$rna,join.rpkm$H3k27ac)
# ggplot(join.rpkm,aes(rna,Ruvbl1)) + geom_smooth() + geom_point() +
#   geom_smooth(aes(rna,Ruvbl2),color='red') +
#   geom_smooth(aes(rna,H3k9me3),color='green')# +
#   scale_x_continuous(limits=c(0,10))+scale_y_continuous(limits=c(0,20))
# 
#   lm.fit=lm(rna~Ruvbl1,data = join)
#   summary(lm.fit)
#   lm.fit$coefficients
#   coef(lm.fit)
#   confint(lm.fit)
#   # predict(lm.fit,data.frame(lstat=(c(5,10,15))),interval = 'confidence')
#   # predict(lm.fit,data.frame(lstat=(c(5,10,15))),interval = 'prediction')
#   plot(join$rna,join$Ruvbl1)
#   abline(lm.fit)
#   abline(lm.fit,lwd=3,col='red')
#   plot(Boston$lstat,Boston$medv,pch=20)
#   plot(Boston$lstat,Boston$medv,pch='+')
#   plot(1:20,1:20,pch=1:20)
#   #image(1:20,1:20,matrix(1:20))
#   par(mfrow=c(2,2))
#   plot(lm.fit)
#   predict(lm.fit)
#   plot(predict(lm.fit),residuals(lm.fit))
#   plot(predict(lm.fit),rstudent(lm.fit))
#   plot(hatvalues(lm.fit))
#   which.max(hatvalues(lm.fit))
#   #--- 多元线性回归
#   lm.fit=lm(medv~lstat+age,data = Boston)
#   summary(lm.fit)
#   par(mfrow=c(2,2))
#   plot(lm.fit)
#   lm.fit=lm(medv~.,data = Boston)
#   library(car)
#   vif(lm.fit)
#   lm.fit=lm(medv~.-age,data = Boston)
#   lm.fit1=update(lm.fit,~.-age)
#   #---添加交互项
#   lm.fit=lm(medv~lstat*age,data=Boston)
#   summary(lm.fit)
#   plot(lm.fit)
#   #--- 预测变量的非线性变换
#   lm.fit2=lm(medv~lstat+I(lstat^2),data=Boston)
#   summary(lm.fit2)
#   lm.fit=lm(medv~lstat,data=Boston)
#   anova(lm.fit,lm.fit2)
#   plot(lm.fit2)
#   lm.fit5=lm(medv~poly(lstat,5),data=Boston)
#   summary(lm.fit5)
  