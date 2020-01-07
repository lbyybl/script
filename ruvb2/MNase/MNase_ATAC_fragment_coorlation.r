# setwd('/WORK/lbyybl/WH/rvb/cor_uniq/MNase/graph/TR_change')
# data <- fread('Pol2_tr.bed')
# diff <- data.frame('dox05'=rep(0,nrow(data)),
#                    'dox1'=rep(0,nrow(data)))
# diff$dox05 <- data$V5-data$V4
# diff$dox1 <- data$V6-data$V4
# diff_matrix <- as.matrix(diff)
# rownames(diff_matrix) <- data$V7
# mean(abs(diff_matrix))
# a <- pheatmap::pheatmap(diff_matrix[abs(diff_matrix[,'dox05'])>7 & abs(diff_matrix[,'dox1'])>7,],scale = 'column')

#---ATAC-seq  质控
#-相关性
setwd('/WORK/lbyybl/WH/rvb/ATAC/sample20191014')
atac <- fread('scores_per_bin10K.tab',skip = 1, col.names = c('chr','st','en',
                                                              'dox1','dox2','h051',
                                                              'h052','h11','h12'))
pdf('atac_dox_cor.pdf',width = 4.3,height = 5)
smoothScatter(log2(atac$dox1+1),log2(atac$dox2+1),xlim = c(0,6),ylim = c(0,6), 
              xlab = 'log2(Dox rep1 density)', ylab = 'log2(Dox rep2 density)')
text(1.5,5.5,paste0('Pearson coorlation = ',round(cor(log2(atac$dox1+1),log2(atac$dox2+1)),3)),cex=0.3)
dev.off()
pdf('atac_h05_cor.pdf',width = 4.3,height = 5)
smoothScatter(log2(atac$h051+1),log2(atac$h052+1),xlim = c(0,6),ylim = c(0,6), xlab = 'log2(Degron 0.5h rep1 density)', ylab = 'log2(Degron 0.5h rep2 density)')
text(1.5,5.5,paste0('Pearson coorlation = ',round(cor(log2(atac$h051+1),log2(atac$h052+1)),3)),cex=0.3)
dev.off()
pdf('atac_h1_cor.pdf',width = 4.3,height = 5)
smoothScatter(log2(atac$h11+1),log2(atac$h12+1),xlim = c(0,6),ylim = c(0,6), xlab = 'log2(Degron 1h rep1 density)', ylab = 'log2(Degron 1h rep2 density)')
text(1.5,5.5,paste0('Pearson coorlation = ',round(cor(log2(atac$h11+1),log2(atac$h12+1)),3)),cex=0.3)
dev.off()

#---- MNase-seq
#-相关性
setwd('/WORK/lbyybl/WH/rvb/cor_uniq/uniq_corbw/MNase')
atac <- fread('scores_per_bin10K.tab',skip = 1, col.names = c('chr','st','en',
                                                              'dox1','dox2','h051',
                                                              'h052','h11','h12'))
pdf('atac_dox_cor.pdf',width = 4.3,height = 5)
smoothScatter(log2(atac$dox1+1),log2(atac$dox2+1),xlim = c(0,1.7),ylim = c(0,1.7), 
              xlab = 'log2(Dox rep1 density)', ylab = 'log2(Dox rep2 density)')
text(0.3,1.5,paste0('Pearson coorlation = ',round(cor(log2(atac$dox1+1),log2(atac$dox2+1)),3)),cex=0.3)
dev.off()
pdf('atac_h05_cor.pdf',width = 4.3,height = 5)
smoothScatter(log2(atac$h051+1),log2(atac$h052+1),xlim = c(0,1.7),ylim = c(0,1.7), xlab = 'log2(Degron 0.5h rep1 density)', ylab = 'log2(Degron 0.5h rep2 density)')
text(0.3,1.5,paste0('Pearson coorlation = ',round(cor(log2(atac$h051+1),log2(atac$h052+1)),3)),cex=0.3)
dev.off()
pdf('atac_h1_cor.pdf',width = 4.3,height = 5)
smoothScatter(log2(atac$h11+1),log2(atac$h12+1),xlim = c(0,1.7),ylim = c(0,1.7), xlab = 'log2(Degron 1h rep1 density)', ylab = 'log2(Degron 1h rep2 density)')
text(0.3,1.5,paste0('Pearson coorlation = ',round(cor(log2(atac$h11+1),log2(atac$h12+1)),3)),cex=0.3)
dev.off()
#--- freagment length distribution
setwd('/WORK/lbyybl/WH/rvb/cor_uniq/MNase/unique')
read_data <- function(file_name){
  freg <- fread(file_name)
  names(freg) <- 'len'
  freg_freq <- freg %>%
    dplyr::group_by(len) %>%
    summarise(n=n())
  freg_freq$freq <- freg_freq$n/nrow(freg)
  return(freg_freq)
}
dox_file <- read_data('R2DoxMNase_merge_freg.txt')
h05_file <- read_data('R205hMNase_merge_freg.txt')
h1_file <- read_data('R21hMNase_merge_freg.txt')
plot(dox_file$len,dox_file$n,type='l',ylab='frequency',xlab='frengment length')
#axis(1,0:1000)
label<-data.frame(x=c(130,290,460,200,390,567),y=rep(0,6))
ggplot(dox_file,aes(x=len,y=freq)) + geom_line(color='red') + scale_x_continuous(limits=c(0,270))+
  geom_line(data=h05_file,aes(x=len,y=freq),color='forestgreen')+
  geom_line(data=h1_file,aes(x=len,y=freq),color='black')+
  # geom_vline(xintercept =c(130,290,460),color='blue')+
  # geom_vline(xintercept = c(200,390,567),color='red')+xlab('frengment length')+
  ylab('frequency')+ #geom_text(data=label,aes(x,y,label=x),size=2)+
  theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(), panel.background=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2))
ggsave('length_distribtion.pdf',width = 5,height = 4)

dox1 <- read_data('R2DoxMNase-1_FKDL190762899-1a_rmdup_uniqe_freg.txt')
dox2 <- read_data('R2DoxMNase-2_FKDL190762900-1a_rmdup_uniqe_freg.txt')
h051 <- read_data('R205hMNase-1_FKDL190762901-1a_rmdup_uniqe_freg.txt')
h052 <- read_data('R205hMNase-2_FKDL190762902-1a_rmdup_uniqe_freg.txt')
h11<- read_data('R21hMNase-1_FKDL190762903-1a_rmdup_uniqe_freg.txt')
h12 <- read_data('R21hMNase-2_FKDL190762904-1a_rmdup_uniqe_freg.txt')
ggplot(dox1,aes(x=len,y=freq)) + geom_line(color='red') + scale_x_continuous(limits=c(0,270))+
  geom_line(data=dox2,aes(x=len,y=freq),color='palevioletred1')+
  geom_line(data=h051,aes(x=len,y=freq),color='forestgreen')+
  geom_line(data=h052,aes(x=len,y=freq),color='mediumseagreen')+
  geom_line(data=h11,aes(x=len,y=freq),color='black')+
  geom_line(data=h12,aes(x=len,y=freq),color='gray19')+
  # geom_vline(xintercept =c(130,290,460),color='blue')+
  # geom_vline(xintercept = c(200,390,567),color='red')+xlab('frengment length')+
  ylab('frequency')+ #geom_text(data=label,aes(x,y,label=x),size=2)+
  theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(), panel.background=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2))
ggsave('rep_length_distribtion.pdf',width = 5,height = 4)
