#---draw traveling ratio curve
# setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/bw_bam/Ruvbl2/merge/bam/Travelrato')
# setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/bw_bam/Pol2/merge/bam/Travelrato')
setwd('/WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl/Pol2/graph/travel')
library(data.table)
library(stringr)
library(readr)
library(dplyr)
library(ggplot2)
read_data <- function(file){
  data <- readr::read_tsv(file, col_names = F)
  data$pro <- data$X4
  data$body <- data$X4
  for (i in 1:nrow(data)){
    data$pro[i] <- str_split(data$X4[i],",")[[1]][1]
    data$body[i] <- str_split(data$X4[i],",")[[1]][2]
  }
  data$pro <- as.numeric(data$pro)
  data$body <- as.numeric(data$body)
  data <- data %>%
    dplyr::filter(body>0)
  #data$pro[1]/data$body[1]
  data$ratio <- data$pro/data$body
  data <- data[order(data$ratio),]
  return(data)
}

rmdup <- function(data){
  data <- data %>% arrange(X1,X2,X3)
  n<-nrow(data)
  for (i in 1:n){
    if (i+1 <= n){
      if (sum(data[i,c(1:3,6)]==data[i+1,c(1:3,6)])==4){
        data[i+1,5] <- data[i,5]
        
      }
    }
    
  }
  return(unique(data))
}

for (i in c('R205hPolII','R21hPolII','R2doxPolII')){
  data <- read_data(paste0(i,'_allgene.tr.bed'))
  data <- rmdup(data)
  assign(paste0(i,'_all'),data)
  # data <- read_data(paste0(i,'_unbindgene.tr.bed'))
  # data <- rmdup(data)
  # assign(paste0(i,'_unbind'),data)
  # data <- read_data(paste0(i,'_bindgene.tr.bed'))
  # data <- rmdup(data)
  # assign(paste0(i,'_bind'),data)
}

# ko_data <- read_data('Pol2_R2_degron.tr.pol2.bed')
# ko_data <- rmdup(ko_data)
# wt_data <- read_data('Pol2_R2_WT.tr.pol2.bed')
# wt_data <- rmdup((wt_data))
ggplot(R2doxPolII_all,aes(x=log(ratio+1),y=cumsum(..count..)/nrow(R2doxPolII_all)))+ 
  geom_freqpoly(bins=nrow(R2doxPolII_all),color='red') + geom_hline(yintercept = 0,linetype = "dashed") + geom_hline(yintercept = 1,linetype = "dashed")+
  # geom_freqpoly(data = R205hPolII_all, aes(x=log(ratio+1),y=cumsum(..count..)/nrow(R205hPolII_all)),bins=nrow(R205hPolII_all))+
  # geom_freqpoly(data = R21hPolII_all, aes(x=log(ratio+1),y=cumsum(..count..)/nrow(R21hPolII_all)),bins=nrow(R21hPolII_all),color='blue')+
  theme(plot.title=element_text(size=40,hjust=0.5), 
        panel.border = element_rect(colour="black",  fill=NA, size=1),
        axis.title=element_text(size=30),
        axis.text=element_text(size=28),
        legend.title=element_text(size=25),
        panel.background=element_blank(),
        legend.text=element_text(size=22))+
  ylab("Percent of Pol II bound gene")+xlab("Traveling ratio")+
  scale_y_continuous(labels=scaleFUN)+
  scale_x_continuous(labels=scaleX)
scaleFUN <- function(x) sprintf("%.0f", x*100)
scaleX <- function(x) sprintf("%.0f", 10^x)
ggsave('Pol2_travel_ratio.pdf',width = 12,height = 6)
ruvbl2 <- read_data('Ruvb2.tr.pol2.bed')
ruvbl2 <- rmdup(ruvbl2)
ggplot(ruvbl2,aes(x=log(ratio+1),y=cumsum(..count..)/nrow(ruvbl2)))+ 
  geom_freqpoly(bins=nrow(ruvbl2)) + geom_hline(yintercept = 0,linetype = "dashed") +
  geom_hline(yintercept = 1,linetype = "dashed")+
  theme(plot.title=element_text(size=40,hjust=0.5), 
        panel.border = element_rect(colour="black",  fill=NA, size=1),
        axis.title=element_text(size=30),
        axis.text=element_text(size=28),
        legend.title=element_text(size=25),
        panel.background=element_blank(),
        legend.text=element_text(size=22))+
  ylab("Percent of Pol II bound gene")+xlab("Traveling ratio")+
  scale_y_continuous(labels=scaleFUN)+
  scale_x_continuous(labels=scaleX)
scaleFUN <- function(x) sprintf("%.0f", x*100)
scaleX <- function(x) sprintf("%.0f", 10^x)
ggsave('RUVBL2_travel_ratio.pdf',width = 12,height = 6)
#data <- data[-which(data$ratio==0),]

#colnames(data)<- c('chr','st','TR','id','score','strand')
head(data)

# data$num <- 1
# data <- unique(data)
# data <- readr::read_tsv('Pol2_R2_WT.tr.pol2.bed',col_names = F)
# data$pro <- data$X4
# data$body <- data$X4
# for (i in 1:nrow(data)){
#   data$pro[i] <- str_split(data$X4[i],",")[[1]][1]
#   data$body[i] <- str_split(data$X4[i],",")[[1]][2]
# }
# data$pro <- as.numeric(data$pro)
# data$body <- as.numeric(data$body)
# #data$pro[1]/data$body[1]
# data$ratio <- data$pro/data$body
# data <- data[order(data$ratio),]
#data$num <- 1:nrow(data)