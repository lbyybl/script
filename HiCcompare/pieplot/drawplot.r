# setwd("/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/test")
# library(data.table)
# library(dplyr)
# library(stringi)
# interafile <- fread("interaction.bed")
# View(interafile)
# interafile <- interafile[,c(1:6,23,28)]
# colnames(interafile) <- c("chr1","start1","end1","chr2","start2","end2","region1","region2")
# interafile$region1 <- str_sub(interafile$region1,1,4)
# interafile$region2 <- str_sub(interafile$region2,1,4)
# interafile$region2[interafile$region2 != "NONG"] <- "GERE"
# interafile <- interafile %>%
#   mutate(regiont = paste0(chr1,"-",start1,"-",end1,"-",chr2,"-",start2,"-",end2),interaction=paste0(region1,"-",region2))
# fwrite(interafile,file = "interafile.bed", sep = "\t")

##--- 将找到的GNG比例用图展示，打算用饼图和累积分布图展示
setwd("/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/graph")
library(data.table)
library(ggplot2)
library(scales)
file1 <- "pdf"
file2 <- "jpeg"
ctcf <- "/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/CTCF/Gene_nongeneinteraction.bed"
ctcfbed <- fread(ctcf)
yy1 <- "/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/YY1/Gene_nongeneinteraction.bed"
yy1bed <- fread(yy1)
pol1 <- "/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/pol1/Gene_nongeneinteraction.bed"
pol1bed <- fread(pol1)
pol2 <- "/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/pol2/Gene_nongeneinteraction.bed"
pol2bed <- fread(pol2)
pol3 <- "/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/pol3/Gene_nongeneinteraction.bed"
pol3bed <- fread(pol3)
pol1vspol2 <- "/DATA/work/lbyybl/ypjiang/extract_matrix/Hiccompare/pol1vspol2/Gene_nongeneinteraction.bed"
pol1vspol2bed <- fread(pol1vspol2)

#--- 将三列式数据转成一列
transformat <- function(data, sample_name){
  data$GG[data$GG!=""] <- "GG"
  data$GN[data$GN!=""] <- "GN"
  data$NN[data$NN!=""] <- "NN"
  data <- data.frame("class"=c(data$GG[data$GG!=""],data$GN[data$GN!=""],data$NN[data$NN!=""]))
  data$sample <- sample_name
  return(data)
}
#data <- CTCFbed
pol1_trans <- transformat(pol1bed,"pol1") #2964
pol2_trans <- transformat(pol2bed,"pol2") #2604
pol3_trans <- transformat(pol3bed,"pol3") #3078
CTCF_trans <- transformat(ctcfbed,"CTCF") #2751
YY1_trans <- transformat(yy1bed,"YY1") # 2311
pol1vspol2_trans <- transformat(pol1vspol2bed,"pol1vspol2") # 2008
total_data <- rbind(pol1_trans,pol2_trans,pol3_trans,CTCF_trans,YY1_trans,pol1vspol2_trans)
#--- 画bar plot
drawbarplot <- function(data){
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold")
    )  
  ggplot(total_data,aes(sample,fill=class)) + geom_bar() + blank_theme 
}
savefile2 <- function(name){
  ggsave(paste0(name,".jpeg"), width = 5,height = 5,dpi = 1000)
  ggsave(paste0(name,".pdf"), width = 5,height = 5,dpi = 1000)
}
drawbarplot(total_data)
#p + geom_text(aes(y = value,label = percent(value/sum(value))), size=5)
savefile2("total_data")
# ggplot(pol1_trans,aes(class)) + geom_bar()+coord_polar()
# ggplot(total_data,aes(sample,fill=class)) + geom_bar() + blank_theme
#--- draw pie plot  
  #--- trans data format
  trans_pie_formt <- function(data,name){
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

  pol1_pietrans <- trans_pie_formt(pol1bed,"pol1")
  pol2_pietrans <- trans_pie_formt(pol2bed,"pol2")
  pol3_pietrans <- trans_pie_formt(pol3bed,"pol3")
  CTCF_pietrans <- trans_pie_formt(ctcfbed,"CTCF")
  YY1_pietrans <- trans_pie_formt(yy1bed,"YY1")
  pol1vspol2bed_pietrans <- trans_pie_formt(pol1vspol2bed,"pol1vspol2")
  total_pie_data <- rbind(pol1_pietrans,pol2_pietrans,pol3_pietrans,CTCF_pietrans,YY1_pietrans,pol1vspol2bed_pietrans)  
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
      geom_text(aes(y = (value)*5/5,label = percent(value/sum(value))), size=5)
  }
  savefile <- function(name){
    ggsave(paste0(name,".jpeg"), width = 5,height = 5,dpi = 1000)
    ggsave(paste0(name,".pdf"), width = 5,height = 5,dpi = 1000)
  }
  drawpieplot(pol1_pietrans)
  savefile("pol1")
  drawpieplot(pol2_pietrans)
  savefile("pol2")
  drawpieplot(pol3_pietrans)
  savefile("pol3")
  drawpieplot(CTCF_pietrans)
  savefile("CTCF")
  drawpieplot(YY1_pietrans)
  savefile("YY1")
  drawpieplot(pol1vspol2bed_pietrans)
  savefile("kpnawvspol2bed")
#  drawpieplot(pol1_pietrans)

 
  
  
  
  
  
  # bp<- ggplot(pol1_pietrans, aes(x="", y=value, fill=class))+
  #   geom_bar(width = 1, stat = "identity")
  # pie <- bp + coord_polar("y", start=0)
  # pie
  # pie + scale_fill_brewer(palette="Dark2")
  # blank_theme <- theme_minimal()+
  #   theme(
  #     axis.title.x = element_blank(),
  #     axis.title.y = element_blank(),
  #     panel.border = element_blank(),
  #     panel.grid=element_blank(),
  #     axis.ticks = element_blank(),
  #     plot.title=element_text(size=14, face="bold")
  #   )  
  # 
  # pie + scale_fill_brewer(palette="Dark2") + blank_theme +
  #   theme(axis.text.x=element_blank()) +
  #   geom_text(aes(y = value + c(0, cumsum(value)[-length(value)]), 
  #                 label = percent(value/sum(value))), size=5)
  # pie + scale_fill_brewer(palette="Dark2") + blank_theme +
  #   theme(axis.text.x=element_blank()) +
  #   geom_text(aes(y = value*4/5,label = percent(value/sum(value))), size=5) #+
  # scale_y_discrete(limit=c("GG","GN","NN"))
  # 
  
  
  
