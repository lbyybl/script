setwd("/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/MAplot/Hiccups")
#setwd("/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/MAplot/Hiccups")
# 读取表达量的表格
options(stringsAsFactors = FALSE)
library(edgeR)
library(ggplot2)
library(data.table)
library(dplyr)
#--- read 4 file and merge to one file
wt1 <- 'pol2_h3k27_hichip_WT_1_loop.bed'
wt2 <- 'pol2_h3k27_hichip_WT_2_loop.bed'
ko1 <- 'pol2_h3k27_hichip_KO_1_loop.bed'
ko2 <- 'pol2_h3k27_hichip_KO_2_loop.bed'
# wt1 <- 'pol2_oceanc_WT_1_loop.bed'
# wt2 <- 'pol2_oceanc_WT_2_loop.bed'
# ko1 <- 'pol2_oceanc_KO_1_loop.bed'
# ko2 <- 'pol2_oceanc_KO_2_loop.bed'
readfile <- function(file,name){
  data <- read.table(file,header = FALSE)
  colnames(data) <- c("loc",name)
  return(data)
}
wt1_file <- readfile(wt1,"wt1")
wt2_file <- readfile(wt2,"wt2")
ko1_file <- readfile(ko1,"ko1")
ko2_file <- readfile(ko2,"ko2")

wt_merge <- merge(wt1_file,wt2_file,by="loc", sort=TRUE)#,all=TRUE)
ko_merge <- merge(ko1_file,ko2_file,by="loc", sort=TRUE)#,all=TRUE)
all_merge <- unique(merge(wt_merge,ko_merge,by="loc"))#,all=TRUE))
all_merge[is.na(all_merge)] <- 0
counts <- as.matrix(all_merge[,2:5])
rownames(counts) <- as.vector(all_merge$loc)
counts <- as.data.frame(counts)
# counts <- read.table(
#   "endgr_input.bed",
#   header=F,
#   sep="\t",
#   row.names=1,
#   comment.char="",
#   check.names=F)
#colnames(counts) <- c("ko","wt")
counts$ko1 <- round(counts$ko1*100)
counts$ko2 <- round(counts$ko2*100)
counts$wt1 <- round(counts$wt1*100)
counts$wt2 <- round(counts$wt2*100)
# 设置样本分组
groups <- factor(c(1,1,2,2))

# 构建edgeR中的对象
y <- DGEList(counts=counts,group=groups,remove.zeros=TRUE)

y <- calcNormFactors(y)
design <- model.matrix(~groups)
y <- estimateDisp(y,design)
et <- exactTest(y)
res <- et$table
#write.table(res, "edgeR.xls", header = T, col.names = NA, sep = "\t" )
ma_input <- res
ma_input$logCPM <- log2((counts$ko1+counts$ko2)/100*(counts$wt1+counts$wt2)/100)/2
ma_input$logFC <- log2((counts$ko1+counts$ko2)/(counts$wt1+counts$wt2))
# ma_input <- data.frame("logCPM"=log((counts$ko1+counts$ko2)/100*(counts$wt1+counts$wt2)/100)/2,
#                        "logFC"=res$logFC,
#                        "PValue"=res$PValue )
ma_input$logFC <- ma_input$logFC
ma_input$change <- ma_input$logFC
ma_input[which(ma_input$logFC < 0),]$change <- "down"
ma_input[which(ma_input$logFC > 0),]$change <- "up"
ma_input[which(ma_input$logFC >= 0 & ma_input$logFC <= 0),]$change <- "con"
ma_input$PValue <- as.numeric(ma_input$PValue)
ma_input$FDR <- p.adjust(ma_input$PValue,method="fdr",n=length(ma_input$PValue))
ma_input$sig <- ma_input$FDR
ma_input[which(ma_input$FDR < 0.1),]$sig <- "s"
ma_input[which(ma_input$FDR >= 0.1 ),]$sig <- "ns"
ma_input$color <- paste0(ma_input$change,ma_input$sig)
#ifelse(ma_input$PValue < 0.0001,ma_input$PValue <- "s",ma_input$PValue <- "ns")
#ma_input[which(ma_input$PValue < 0.0001),]$PValue <- "s"
#ma_input[which(ma_input$PValue != 's'),]$PValue <- "ns"
#plotMA(ma_input,colNonSig="black",colSig="red")
ma_input %>%
  group_by(change) %>%
  summarise(class_n=n())

#sub_input <- ma_input[1:1000,]
ggplot(ma_input,aes(x=logCPM,y=logFC,color=color))+ geom_point(size=0.1,alpha=0.7) +
  geom_hline(aes(yintercept = 0),colour="black",size=0.2,linetype="dashed")+
  #geom_hline(aes(yintercept = log2(2)),color="cyan3",size=0.2,linetype="dashed")+
  #geom_hline(aes(yintercept = log2(0.5)),color="cyan3",size=0.2,linetype="dashed")+
  labs(x="log2 (Degron*untreated)/2",y="log2 Degron/untreated FC",color="") +
  scale_y_continuous(expand=c(0,0),limits = c(-4,4)) +
  #scale_color_manual(values = c("grey67","blue","blue","red","red"))+
  scale_color_manual(values = c("red","red","red","red","red"))+
  #scale_color_manual(values = c("grey67","grey67","blue","red"))+
  scale_x_continuous(limits = c(1.5,6))+
  theme(
    legend.key.width = unit(0.1,"line"),
    legend.key.height = unit(1,"cm"),
    legend.text = element_text(
      size = 10,
      hjust = 0,
      #face = "italic",
      colour = "black",
      angle = 0
    ),
    
    plot.title=element_text(colour="black", size=10),
    axis.line.x = element_blank(),#element_line(color="black", size = 0.5),
    axis.line.y = element_blank(), #element_line(color="black", size = 0.5),
    axis.text.x=element_text(angle=0, hjust=0.5,vjust = 0,
                             colour="black", size=12),
    axis.text.y=element_text(#face = "bold",
      colour="black", size=12),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    panel.background=element_blank()
  )

filename="pol2_hichip_hicups"
fwrite(ma_input,paste0(filename,".bed"),sep = "\t",row.names = TRUE)
ggsave(paste0(filename,".jpeg"),width =4, height = 4 )
ggsave(paste0(filename,".pdf"),width =4, height = 4 )

