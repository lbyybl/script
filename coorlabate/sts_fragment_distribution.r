# #---统计杨博的文章中的不同 酶切位点的片段的distribution；
# setwd('/DATA2/work/lbyybl/coorlaborate/YB/fragemnt_distribution')
# read_data <- function(file_name){
#   freg <- fread(file_name,col.names = c('chr','st','en','id','sc','strand'))
#   freg$len <- freg$en-freg$st
#   freg_freq <- freg %>%
#     dplyr::group_by(len) %>%
#     summarise(n=n())
#   freg_freq$freq <- freg_freq$n/nrow(freg)
#   return(freg_freq)
# }
# ALU1 <- read_data('hg19_HBV_ALU1.bed')
# EcoRI <- read_data('hg19_HBV_EcoRI.bed')
# MboI <- read_data('hg19_HBV_MboI.bed')
# NdeI <- read_data('hg19_HBV_NdeI.bed')
# XhoI <- read_data('hg19_HBV_XhoI.bed')
# plot(MboI$len,MboI$n,type='l',ylab='frequency',xlab='frengment length',xlim=c(1,1000))
# #axis(1,0:1000)
# 
# ggplot(ALU1,aes(x=len,y=freq)) + geom_line(color='red') + scale_x_continuous(limits=c(0,1000))+
#   geom_line(data=EcoRI,aes(x=len,y=freq),color='palevioletred1')+
#   geom_line(data=MboI,aes(x=len,y=freq),color='forestgreen')+
#   geom_line(data=NdeI,aes(x=len,y=freq),color='mediumseagreen')+
#   geom_line(data=XhoI,aes(x=len,y=freq),color='black')+
#   # geom_vline(xintercept =c(130,290,460),color='blue')+
#   # geom_vline(xintercept = c(200,390,567),color='red')+xlab('frengment length')+
#   ylab('frequency')+ #geom_text(data=label,aes(x,y,label=x),size=2)+
#   theme(plot.title=element_text(hjust=0.5),
#         panel.grid=element_blank(), panel.background=element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.2))
# ggsave('enzyme_length_distribtion.pdf',width = 5,height = 4)

#---统计杨博的文章中的不同 酶切位点的片段的distribution；
setwd('/DATA2/work/lbyybl/coorlaborate/YB/fragemnt_distribution')
read_data <- function(file_name){
  freg <- fread(file_name,col.names = c('chr','st','en','id','sc','strand'))
  freg$len <- freg$en-freg$st

  return(freg)
}
ALU1 <- read_data('hg19_HBV_ALU1.bed')
EcoRI <- read_data('hg19_HBV_EcoRI.bed')
MboI <- read_data('hg19_HBV_MboI.bed')
NdeI <- read_data('hg19_HBV_NdeI.bed')
XhoII <- read_data('hg19_HBV_XhoII.bed')
median(ALU1$len)
median(EcoRI$len)
median(MboI$len)
median(NdeI$len)
median(XhoII$len)
#plot(MboI$len,MboI$n,type='l',ylab='frequency',xlab='frengment length',xlim=c(1,1000))
#axis(1,0:1000)

ggplot(ALU1,aes(x=len)) + geom_freqpoly(color='red',binwidth = 10) + 
  scale_x_continuous(limits=c(0,3000))+
  geom_freqpoly(data=EcoRI,aes(x=len),color='purple',binwidth = 10)+
  geom_freqpoly(data=MboI,aes(x=len),color='yellow',binwidth = 10)+
  geom_freqpoly(data=NdeI,aes(x=len),color='blue',binwidth = 10)+
  geom_freqpoly(data=XhoII,aes(x=len),color='black',binwidth = 10)+
  # geom_vline(xintercept =c(130,290,460),color='blue')+
  # geom_vline(xintercept = c(200,390,567),color='red')+xlab('frengment length')+
  ylab('frequency')+ #geom_text(data=label,aes(x,y,label=x),size=2)+
  theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(), panel.background=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2))
ggsave('enzyme_length_distribtion.pdf',width = 5,height = 4)
