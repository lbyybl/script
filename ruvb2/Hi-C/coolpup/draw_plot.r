library(ggplot2)
setwd('/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/loop/APA')
wt <- 'wt_25k.txt'
ko <- 'ko_25k.txt'
#file <- read.table('wt_25000-25.0K_over_tad_10-shifts.np.txt')
data <- read.table('wt_25k.txt')
APA_cal <- function(file){
  data <- read.table(file)
  APA <- as.numeric(data[11,11])/mean(c(data[17:21,17],data[17:21,18],data[17:21,19],data[17:21,20],data[17:21,21]))
  return(APA)
}
APA_cal(wt)
APA_cal(ko)
read_data <- function(data){
  file <- read.table(data)
  long_file <- melt(file)
  df <- data.frame("bin2"=rep(1:nrow(file),each=nrow(file)),
                   'bin1'=rep(1:nrow(file),nrow(file)),
                   'int'=long_file$value)
  return(df)
}
wt_file <- read_data(wt_file)
ko_file <- read_data(ko_file)
max_va <- (max(log(wt_file$int + 1),log(ko_file$int + 1)))*1.01
min_va <- (min(log(wt_file$int + 1),log(ko_file$int + 1)))*0.99
# file <- read.table('test')
# long_file <- melt(file)
# #x <- rester(nrow=nrow(dense), ncol=nrow(dense))
# #x[] <- long_file$value
# df <- data.frame("bin2"=rep(1:nrow(file),each=nrow(file)),
#                         'bin1'=rep(1:nrow(file),nrow(file)),
#                         'int'=long_file$value)
# # df <- data.frame("bin1"=rep(1:nrow(file),each=nrow(file)),
#                  'bin2'=rep(1:nrow(file),nrow(file)),
#                  'int'=0)

# for (i in 1:nrow(file)){
#   for (j in 1:ncol(file)){
#     df[(i-1)*nrow(file)+j,3] <- file[i,j]
#   }
# }
ggplot(data=wt_file, aes(x=bin1, y=bin2, fill=log(int+1)))+geom_raster()+
  #ggtitle(paste0("Averaged Interaction within TAD\n(Ruvbl2 ", arg, " merge)"))+
  scale_fill_distiller(palette = "RdYlBu",direction = -1 ,values=c(0,0.45,0.5,1),limit=c(min_va,max_va)
  )+labs(fill="log(interaction+1)")+
  scale_x_continuous(expand=c(0, 0))+
  scale_y_reverse(expand=c(0, 0))+
  #scale_fill_continuous()+
  theme(plot.title=element_text(size=40,hjust=0.5), panel.border = element_rect(colour="black",  fill=NA, size=1), 
        axis.title=element_text(size=30),axis.text=element_text(size=28),legend.title=element_text(size=25),legend.text=element_text(size=22))+
  xlab("-100kb 0 100kb")+ylab("-100kb 0 100kb")
ggsave('wt_APA.pdf',width = 7,height = 4)

ggplot(data=ko_file, aes(x=bin1, y=bin2, fill=log(int+1)))+geom_raster()+
  #ggtitle(paste0("Averaged Interaction within TAD\n(Ruvbl2 ", arg, " merge)"))+
  scale_fill_distiller(palette = "RdYlBu",direction = -1 ,values=c(0,0.45,0.5,1),limit=c(min_va,max_va)
  )+labs(fill="log(interaction+1)")+
  scale_x_continuous(expand=c(0, 0))+
  scale_y_reverse(expand=c(0, 0))+
  #scale_fill_continuous()+
  theme(plot.title=element_text(size=40,hjust=0.5), panel.border = element_rect(colour="black",  fill=NA, size=1), 
        axis.title=element_text(size=30),axis.text=element_text(size=28),legend.title=element_text(size=25),legend.text=element_text(size=22))+
  xlab("-100kb 0 100kb")+ylab("-100kb 0 100kb")
ggsave('ko_APA.pdf',width = 7,height = 4)
diff <- ko_file
diff$int <- ko_file$int - wt_file$int
ggplot(data=diff, aes(x=bin1, y=bin2, fill=int))+geom_raster()+
  #ggtitle(paste0("Averaged Interaction within TAD\n(Ruvbl2 ", arg, " merge)"))+
  scale_fill_distiller(palette = "RdYlBu",direction = -1 ,values=c(0,0.45,0.5,1)#,limit=c(min_va,max_va)
  )+labs(fill="log(interaction+1)")+
  scale_x_continuous(expand=c(0, 0))+
  scale_y_reverse(expand=c(0, 0))+
  #scale_fill_continuous()+
  theme(plot.title=element_text(size=40,hjust=0.5), panel.border = element_rect(colour="black",  fill=NA, size=1), 
        axis.title=element_text(size=30),axis.text=element_text(size=28),legend.title=element_text(size=25),legend.text=element_text(size=22))+
  xlab("-100kb 0 100kb")+ylab("-100kb 0 100kb")
ggsave('diff_APA.pdf',width = 7,height = 4)
