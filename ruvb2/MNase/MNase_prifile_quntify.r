#---- MNase profile 
setwd('/WORK/lbyybl/WH/rvb/cor_uniq/MNase/graph/TSS_profile/deeptools/merge/')
trans_df <- function(data,name){
  data <- data %>%
    summarise_all(mean)
  data <- t(data)
  colnames(data) <- "value"
  data <- as.data.frame(data)
  data$sample <- name
  data$coor <- ((0:399)*10-2000)
  return(data)
}
rbind_df <- function(data1,data2,data3){
  med <- median(c(data1$value,data2$value,data3$value))
  data1$value <- (data1$value/median(data1$value))*med
  data2$value <- (data2$value/median(data2$value))*med
  data3$value <- (data3$value/median(data3$value))*med
  merge <- rbind(data1,data2)
  merge<-rbind(merge,data3)
  return(merge)
}
tss_profile <- fread('dox1h_top3000.gz',skip = 1,stringsAsFactors = F)
names(tss_profile)
dox <- tss_profile[,7:406]
h05 <- tss_profile[,407:806]
h1 <- tss_profile[,807:1206]
h05_mean<- trans_df(h05,'h05')
h1_mean <-trans_df(h1,'h1')
dox_mean <- trans_df(dox,'dox')
plot(dox_mean$coor,dox_mean$value,type='l')

merge <- rbind(dox_mean,h05_mean)
merge <- rbind(merge,h1_mean)
#merge<-rbind_df(dox_mean,h05_mean,h1_mean)

ggplot(merge,aes(coor,value,color=sample))+geom_line()+
  scale_x_continuous(limits=c(-1000,1000))+xlab('')+
  ylab('Read density(normalized)')+ 
  #scale_y_continuous(limits = c(0.4,1),expand=c(0,0))+
  theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(), panel.background=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2))#+
  # geom_vline(xintercept =c(130,330,520,700,880,-210,-400,-590,-780,-960),color='blue')+
  # geom_vline(xintercept = c(-305,-495,-680,-860,230,420,610,790),color='red')
ggsave('dox1h_top3000.pdf',width=7,height = 3)
tss_min <- c(-305,-495,-680,-860,230,420,610,790)
tss_max <- c(130,330,520,700,880,-210,-400,-590,-780,-960)
ctcf_max <-c(130,330,520,700,880,-130,-330,-520,-700,-880)
ctcf_min <- c(-230,-420,-610,-790,230,420,610,790)
#plot(1:400,dox_mean[1,],type='l')
#---quantify
max_local <- function(data,loc){
  data <- data %>%
    filter(coor>(loc-50) & coor <(loc+50))
  return(max(data$value))
}
min_local <- function(data,loc){
  data <- data %>%
    filter(coor>(loc-50) & coor <(loc+50))
  return(min(data$value))
}
find_real_coor <- function(data,coor_df,type){
  if (type=='max'){
    for (i in 1:nrow(coor_df)){
      va <- max_local(data,coor_df$cand[i])
      coor_df$value[i]<- va
      loc<-which(data$value==va)
      coor_df$real_coor[i]<-data$coor[loc]
    }
  }else if (type=='min'){
    for (i in 1:nrow(coor_df)){
      va <- min_local(data,coor_df$cand[i])
      coor_df$value[i]<- va
      loc<-which(data$value==va)
      coor_df$real_coor[i]<-data$coor[loc]
    }
  }
  return(coor_df)
}
find_value <- function(data,coor){
  data2 <- data.frame('coor'=coor,
                      'value'=0)
  for (i in 1:nrow(data2)){
    loc <- which(data$coor==data2$coor[i])
    data2$value[i] <- data$value[loc]
  }
  data2<-data2%>%
    arrange(data2$coor)
  return(data2)
}
#----find coor to quantify
all <- dox_mean
all$value <- dox_mean$value+h05_mean$value+h1_mean$value
cand_max <- data.frame('cand'=c(130,330,520,700,880,-210,-400,-590,-780,-960),
                       'value'=0,'real_coor'=0)
real_max <- find_real_coor(dox_mean,cand_max,'max')

cand_min <- data.frame('cand'=c(-305,-495,-680,-860,230,420,610,790),
                       'value'=0,'real_coor'=0)
real_min <- find_real_coor(dox_mean,cand_min,'min')

dox_max <- find_value(dox_mean,real_max$real_coor)
dox_min <- find_value(dox_mean,real_min$real_coor)
h05_max <- find_value(h05_mean,real_max$real_coor)
h05_min <- find_value(h05_mean,real_min$real_coor)
h1_max <- find_value(h1_mean,real_max$real_coor)
h1_min <- find_value(h1_mean,real_min$real_coor)
find_change <- function(dox_max,h05_max,h1_max){
  names(dox_max) <- c('coor','dox')
  names(h05_max) <- c('coor','h05')
  names(h1_max) <- c('coor','h1')
  merge <- merge(dox_max,h05_max)
  merge_max <- merge(merge,h1_max)
  merge_max$dox_05 <- merge_max$h05/merge_max$dox
  merge_max$dox_1 <- merge_max$h1/merge_max$dox
  merge_max$h05_1 <- merge_max$h1/merge_max$h05
  return(merge_max)
}

merge_max <- find_change(dox_max,h05_max,h1_max)
merge_min <- find_change(dox_min,h05_min,h1_min)
fwrite(merge_max,'tss_unbinding_max.csv')
fwrite(merge_min,'tss_unbinding_min.csv')

# c <- rbind(dox_max,dox_min)
# ggplot(c,aes(coor,value))+geom_line()+geom_line(data = dox_mean,aes(coor,value))
