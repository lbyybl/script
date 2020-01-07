##--- 将找到的GNG比例用图展示，打算用饼图和累积分布图展示
setwd("/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/graph")
library(data.table)
library(ggplot2)
library(dplyr)
library(scales)
file1 <- "pdf"
file2 <- "jpeg"
CTCF <- "/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/CTCF/Gene_nongeneinteraction.bed"
YY1 <- "/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/YY1/Gene_nongeneinteraction.bed"
pol1 <- "/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/pol1/Gene_nongeneinteraction.bed"
pol2 <- "/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/pol2/Gene_nongeneinteraction.bed"
pol3 <- "/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/pol3/Gene_nongeneinteraction.bed"
pol1vspol2 <- "/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/pol1vspol2/Gene_nongeneinteraction.bed"
pol1vspol3 <- "/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/pol1vspol3/Gene_nongeneinteraction.bed"
pol2vspol3 <- "/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/pol2vspol3/Gene_nongeneinteraction.bed"
#--- 将三列式数据转成一列
transformat <- function(file, sample_name){
  data <- fread(file)
  data$GG[data$GG!=""] <- "GG"
  data$GN[data$GN!=""] <- "GN"
  data$NN[data$NN!=""] <- "NN"
  data <- data.frame("class"=c(data$GG[data$GG!=""],data$GN[data$GN!=""],data$NN[data$NN!=""]))
  data$sample <- sample_name
  return(data)
}
#data <- CTCFbed
total_data <- data.frame('class'=NULL,
                         'sample'=NULL)
sample_list <- c("pol1","pol2","pol3","CTCF","YY1","pol1vspol2","pol1vspol3","pol2vspol3")

for (i in sample_list){
  assign(paste0(i,'_trans'),transformat(get(i),i))
  total_data <- bind_rows(total_data,get(paste0(i,'_trans')))
}

# pol1_trans <- transformat(pol1,"pol1") 
# pol2_trans <- transformat(pol2,"pol2") 
# pol3_trans <- transformat(pol3,"pol3") 
# CTCF_trans <- transformat(ctcf,"CTCF") 
# YY1_trans <- transformat(yy1,"YY1") 
# pol1vspol2_trans <- transformat(pol1vspol2,"pol1vspol2") 
# 
# total_data <- rbind(pol1_trans,pol2_trans,pol3_trans,CTCF_trans,YY1_trans,pol1vspol2_trans)
#--- 画bar plot
drawbarplot <- function(data){
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      #panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold"),
      axis.text.x=element_text(angle=30, hjust=1.2,vjust = 1.2, #face = "bold",
                                 colour="black", size=12),
      axis.text.y=element_text(#face = "bold",
        colour="black", size=12),
      panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    )  
  ggplot(total_data,aes(sample,fill=class)) + geom_bar() + blank_theme +
    scale_x_discrete(limits=c('pol1vspol2','pol1vspol3','pol2vspol3','pol1',
                              'pol2','pol3','CTCF','YY1'))
}
savefile2 <- function(name){
  ggsave(paste0(name,".jpeg"), width = 5,height = 5,dpi = 1000)
  ggsave(paste0(name,".pdf"), width = 5,height = 5,dpi = 1000)
}
drawbarplot(total_data)

savefile2("total_data")

#test <- total_data %>% group_by(class,sample) %>% summarise(num=n())


#--- draw pie plot  
#--- trans data format
trans_pie_formt <- function(file,name){
  data <- fread(file)
  data$GG[data$GG!=""] <- "GG"
  rowGG <- length(data$GG[data$GG!=""])
  data$GN[data$GN!=""] <- "GN"
  rowGN <- length(data$GN[data$GN!=""])
  data$NN[data$NN!=""] <- "NN"
  rowNN <- length(data$NN[data$NN!=""])
  data <- data.frame("class"=c("GG","GN","NN"),
                     "value"=c(rowGG,rowGN,rowNN))
  data$sample <- name
  return(data)
}

total_pie_data <- data.frame('class'=NULL,
                             'value'=NULL,
                             'sample'=NULL)
#sample_list <- c("pol1","pol2","pol3","CTCF","YY1","pol1vspol2")

for (i in sample_list){
  assign(paste0(i,'_pietrans'),trans_pie_formt(get(i),i))
  total_pie_data <- bind_rows(total_pie_data,get(paste0(i,'_pietrans')))
}
# pol1_pietrans <- trans_pie_formt(pol1,"pol1")
# pol2_pietrans <- trans_pie_formt(pol2,"pol2")
# pol3_pietrans <- trans_pie_formt(pol3,"pol3")
# CTCF_pietrans <- trans_pie_formt(ctcf,"CTCF")
# YY1_pietrans <- trans_pie_formt(yy1,"YY1")
# pol1vspol2bed_pietrans <- trans_pie_formt(pol1vspol2,"kpnawvspol2")
# total_pie_data <- rbind(pol1_pietrans,pol2_pietrans,pol3_pietrans,CTCF_pietrans,YY1_pietrans,pol1vspol2bed_pietrans)  
#--- draw plot
#ggplot(pol1_pietrans,aes(x="",y=value,fill=class)) + geom_bar(width = 1) + coord_polar("y")

drawpieplot <- function(data){
  bp<- ggplot(data, aes(x="", y=value, fill=class))+
    geom_bar(width = 1, stat = "identity")
  pie <- bp + coord_polar("y", start=0)
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold")
    )  
  pie + scale_fill_brewer(palette="Dark2") + blank_theme +
    theme(axis.text.x=element_blank()) +
    geom_text(aes(y =(sum(value)-cumsum(value))+ sum(value)/7,label = percent(value/sum(value))), size=5)
}
savefile <- function(name){
  ggsave(paste0(name,".jpeg"), width = 5,height = 5,dpi = 1000)
  ggsave(paste0(name,".pdf"), width = 5,height = 5,dpi = 1000)
}
for (i in sample_list){
  drawpieplot(get(paste0(i,'_pietrans')))
  savefile(i)
}
total_pie_data %>% 
  group_by(sample) %>%
  summarise(n=sum(value))

# drawpieplot(pol1_pietrans)
# savefile("pol1")
# drawpieplot(pol2_pietrans)
# savefile("pol2")
# drawpieplot(pol3_pietrans)
# savefile("pol3")
# drawpieplot(CTCF_pietrans)
# savefile("CTCF")
# drawpieplot(YY1_pietrans)
# savefile("YY1")
# drawpieplot(pol1vspol2bed_pietrans)
# savefile("kpnawvspol2bed")

#----------------
# sample region
#----------------
setwd('/home/boyuanli/bashscript/bin/HiCcompare/sample/sample_pair_result')
#sample100 <- 'sample100sample_pair/sample100sample_pair.txt'
summ_data <- function(file){
  data <- fread(file,col.names = c('GG','GN','NN'))
  NN_num <- sum(data$NN!='')
  GN_num <- sum(data$GN!='')
  GG_num <- sum(data$GG!='')
  data_summ <- data.frame('GG'=GG_num,
                          'GN'=GN_num,
                          'NN'=NN_num)
  return(data_summ)
}

readsamplefile <- function(name){
  assign(paste0('sample',name),paste0('sample',name,'sample_pair/sample',name,'sample_pair.txt'))
  data_summ <- summ_data(get(paste0('sample',name)))
  return(data_summ)
}
total_data_summ <- data.frame('GG'=NULL,
                              'GN'=NULL,
                              'NN'=NULL)
for (i in 1:100){
  assign(paste0('sample',i,'_summ'),readsamplefile(i))
  total_data_summ <- bind_rows(total_data_summ,get(paste0('sample',i,'_summ')))
}
  

# mean_data_summ <- data.frame('GG'=mean(total_data_summ$GG),
#                              'GN'=mean(total_data_summ$GN),
#                              'NN'=mean(total_data_summ$NN))

mean_data_summ <- data.frame('class'=c('GG','GN','NN'),
                             'num'=c(mean(total_data_summ$GG),
                                     mean(total_data_summ$GN),
                                     mean(total_data_summ$NN)))

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )  
ggplot(mean_data_summ,aes(x='',y=num,fill=class)) + geom_col() +
  geom_text(aes(x='',y = (sum(num)-cumsum(num))+ sum(num)/7,label = percent(num/sum(num))), size=5) +
  #scale_fill_discrete(limits=c('GG','GN','NN')) + 
  coord_polar("y", start=0) + blank_theme +
  theme(axis.text.x=element_blank()) + scale_fill_brewer(palette="Dark2")
savefile('random')
# ggplot(mean_data_summ,aes(x='',y=num,fill=class)) + geom_col() + 
#   geom_text(aes(x='',y=(sum(num)-cumsum(num))+ sum(num)/7,label=percent(num/sum(num)))) + 
#   scale_fill_discrete(limits=c('NN','GN','GG'),
#                       guide = guide_legend(reverse = TRUE)) + 
#   coord_polar("y", start=0) + blank_theme + 
#   theme(axis.text.x = element_blank())
  


  