rm(list=ls())

library(plyr)
library(data.table)
library(ggplot2)
source("/home/khlun/R_function/multiplot.R")


narm<-function(dataframe){
  loc<-which(is.na(dataframe[,3])==T)
  if(length(loc)!=0){
    dataframe<-dataframe[-loc,]
  }
  dataframe
}

setwd("/home/khlun/HMG_Info")


seq_tpm<-fread("../gtex/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct",skip=2,header=T,sep="\t",quote="")
names(seq_tpm)<-gsub("\\.","_",names(seq_tpm))

minit<-seq_tpm[1:20,1:20]
View(minit)

# subset HMG gene data

gene_match<-"(^TRMT112|^FARSA|^RRM1|^ATP2A2|^RUVBL2|^DNM1L|^IARS|^DYNC1H1|^SF3B3|^PRMT5|^ALDOA|^SNRPD3|^BUB3|^GNB2L1|^PFAS|^ATIC|^ATP6V1A|^PCNA|^VARS|^PRPF19)"
gene_sub_group<-""

pick <- seq_tpm[grep(paste0(gene_match,gene_sub_group),seq_tpm$Description, ignore.case = T),]

View(pick[grep("AS[0-9]$", invert = T, pick$Description), 1:20])

# sum(seq_tpm[,2])
# mean(seq_tpm[,2])

# HMGN3-AS1--HMGN3 antisense RNA 1

View(pick[-grep("(P[0-9]$|AS[0-9])",pick$Description),1:20])
pick_main<-pick[-grep("(P[0-9]$|AS[0-9])",pick$Description),]
pick_main<-pick_main[order(pick_main$Description),]

# read annotation

anno<-fread("/home/khlun/gtex/GTEx_v7_Annotations_SampleAttributesDS.txt",sep="\t",fill=T,quote="",header=T)

ID_by_tissue<-split(ID_all<-data.frame(ID=anno$SAMPID,tissue=anno$SMTSD, stringsAsFactors = F),ID_all$tissue)
names(ID_by_tissue)

# grep(as.character(ID_by_tissue[[1]]$ID)[1],names(pick_main))


# transverse HMG data to join with ID

map_data<-data.frame(t(pick_main))
# dim(map_data)
map_data<-data.frame(cbind(rownames(map_data),map_data), stringsAsFactors = F)
names(map_data)<-c("ID",as.matrix(map_data[2,-1]))
map_data<-map_data[-c(1,2),]
map_data$ID <- as.character(map_data$ID)
map_data[1:10,]

# join ID annotation with tpm data in list

framelist<-lapply(ID_by_tissue, as.data.frame, stringsAsFactors=F)
map_Data<-lapply(framelist, join, map_data, by="ID")
tpm_by_tissue<-lapply(map_Data, narm)

rm(framelist,map_data,map_Data)

# no ID match in tpm data for ID[[23]]
#View(tpm_by_tissue[[23]][,-c(1,2)])
names(tpm_by_tissue)[23]
# "Cells - Leukemia cell line (CML)"


# calculate mean expr of HMGs for each tissue
expr_mean<-data.frame()
for(i in c(1:length(tpm_by_tissue))[-23]){
  tissue_mean<-t(colMeans(apply(as.matrix(tpm_by_tissue[[i]][,-c(1,2)]),2,as.numeric)))
  expr_mean<-rbind(expr_mean,cbind(names(tpm_by_tissue)[i],tissue_mean))
}

rname<-expr_mean[,1]
expr_mean<-expr_mean[,-1]
expr_mean<-as.data.frame(apply(as.matrix(expr_mean),2,as.numeric))
rownames(expr_mean)<-rname

rm(rname,tissue_mean,i)

tpm_thrs<-10

expr_over_thrs<-expr_mean-tpm_thrs
expr_over_thrs[expr_over_thrs<0]<-0
expr_over_thrs[expr_over_thrs>0]<-1
count<-data.frame(Gene=names(expr_over_thrs),percent=colSums(expr_over_thrs)/nrow(expr_over_thrs)*100)

rm(expr_over_thrs)
count<-count[order(-count$percent),]

# plot expr data

expr_graph1<-ggplot(count,aes(x=reorder(count$Gene,-count$percent),y=count$percent))+
  geom_bar(stat="identity",width=.7,fill="steelblue3",alpha=0.9)+
  xlab("HMG Gene Subtype")+ylab("Expression Percentage (%)")+
  theme(axis.text=element_text(size=14),axis.text.x=element_text(angle=30,vjust=0.9,hjust=0.8),
        axis.title = element_text(size=15,face="bold"))
expr_graph1

# png("HMG_expr_tpm.png",width=930,height=550)
# expr_graph1
# dev.off()

# seq_rpkm6<-read.gct("GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct")
# 
# anno6<-read.table("GTEx_Data_V6_Annotations_SampleAttributesDS.txt",sep="\t",fill=T,quote="",header=T)
# 
# ID_by_tissue6<-split(ID_data6<-data.frame(ID=anno6$SAMPID,tissue=anno6$SMTSD),tissue)
# names(ID_by_tissue6)


############# CRISPR score for cell essentiality

cs<-read.csv("table_s3.csv")
cs_gene<-cs[grep(paste0(gene_match,gene_sub_group),cs$Gene, ignore.case = T),]


## 1 ############ KBM7 line

cs_gene$pcut<-0
cs_gene$pcut[cs_gene$KBM7.adjusted.p.value<0.05]<-1
cs_gene$essential<-0
cs_gene$essential[cs_gene$pcut==1 & cs_gene$KBM7.CS<(-1)]<-1
cs_gene<-data.frame(Gene=cs_gene$Gene,CS_KBM7=cs_gene$KBM7.CS,p_value=cs_gene$KBM7.adjusted.p.value,
                   pcut=cs_gene$pcut,essential=cs_gene$essential)


# CswithCons<-join(cs_gene, cons)
# 
# Sum<-join(count,CswithCons)

Sum <- join(count, cs_gene)


# plot CS score distribution

cs$row_num<-c(1:nrow(cs))
csplot<-cs[order(cs$KBM7.CS),]
csplot<-join(csplot,data.frame(Gene=cs_gene$Gene,essential=cs_gene$essential,HMG=1))
csplot[is.na(csplot)==T]<-0

c_score<-ggplot(csplot,aes(x=as.numeric(reorder(csplot$row_num,csplot$KBM7.CS)),y=csplot$KBM7.CS,group=1))+
  geom_line(alpha=0.7)+geom_point(aes(alpha=factor(csplot$HMG),color=factor(csplot$essential)))+
  scale_alpha_manual(values=c(0,0.8),guide=F)+
  scale_color_manual(values=c("#888888","brown2"),labels=c("non-essential","essential"))+
  scale_x_continuous(breaks=seq(0,max(csplot$row_num),5000))+
  scale_y_continuous(breaks=seq(floor(min(csplot$KBM7.CS)),max(csplot$KBM7.CS),1))+
  theme(legend.title=element_blank(),legend.text=element_text(size=11, color="black"),
        plot.title=element_text(hjust=0.5), panel.background = element_blank(),
        axis.text = element_text(color="black", size=12), axis.title=element_text(size=14,face="bold"),
        panel.border = element_rect(color="black", fill=NA, size=0.8), legend.key = element_blank())+
  xlab("genes ranked by CRISPR score")+ylab("CRISPR score")+ggtitle("Cell Essentiality")#+
  annotate("segment", x=2000, xend=800, y=-2.9, yend=-2.38,colour="brown2")+
  annotate("text",label="ZNF131",x=3500,y=-3.1,colour="brown2",size=4)+
  annotate("segment", x=6000, xend=4200, y=-1.3, yend=-0.6,colour="#888888")+
  annotate("text",label="ZBTB14",x=7500,y=-1.5,colour="#888888",size=4)+
  annotate("segment", x=13500, xend=11200, y=-1, yend=-0.16,colour="#888888")+
  annotate("text",label="ZBTB21",x=13700,y=-1.2,colour="#888888",size=4)

c_score

# png("crispr_score.png",width=630,height=430)
# c_score
# dev.off()

# ## 2 ######### cs score summary of all lines
# 
# trans<-function(arr){
#   arr<-arr[-2]
#   cell<-c("KBM7","K562","Jiyoye","Raji")
#   newframe<-data.frame(Gene=arr[1],cs_score=as.numeric(arr[seq(2,8,2)]),p_value=as.numeric(arr[seq(3,9,2)]),cell=cell)
# }
# 
# 
# cs_trans<-data.frame()
# for(i in 1:nrow(cs_gene)){
#   cs_trans<-rbind(cs_trans, trans(cs_gene[i,]))
# }
# 
# cs_trans$pcut_05<-"Insignificant"
# cs_trans$pcut_05[which(cs_trans$p_value<0.05)]<-as.character(cs_trans$cell[which(cs_trans$p_value<0.05)])
# 
# cs_trans$pcut_10<-"Insignificant"
# cs_trans$pcut_10[cs_trans$p_value<0.10]<-as.character(cs_trans$cell[cs_trans$p_value<0.10])
# 
# cs_trans<-cs_trans[-which(cs_trans$cell=="K562"),]
# 
# gene_order_by_expr<-c("HMGB1","HMGB2","HMGN1","HMGN2","HMGN3","HMGN4","HMGA1","HMG20B",
#                       "HMGXB3","HMG20A","HMGXB4","HMGN5","HMGB3","HMGA2","HMGB4")
# 
# 
# 
# cs_point<-ggplot(data=cs_trans,aes(x=Gene,y=cs_score,color=pcut_05))+
#   geom_point()+scale_color_manual(values=c("#999999","brown2","chocolate1","darkorchid4"))+
#   xlab("HMG Gene Subtype")+ylab("CRISPR Score")+
#   scale_x_discrete(limits=gene_order_by_expr)+scale_y_continuous(breaks=seq(-2,0.5,0.5))+
#   theme(axis.text=element_text(size=14),axis.text.x=element_text(angle=30,vjust=0.9,hjust=0.8),
#         legend.title=element_blank(),axis.title = element_text(size=15,face="bold"),
#         legend.text=element_text(size=13))+
#   geom_hline(yintercept=-1,color="#999999")
# # annotate("text",label="essential",x=Inf,y=-1.25,colour="#999999",size=5,hjust=1.1)+
# # annotate("text",label="non essential",x=Inf,y=-0.75,colour="#999999",size=5,hjust=1.07)
# cs_point
# png("crispr_score_modified.png",width=630,height=430)
# cs_point
# dev.off()
# 



# plot expr data with cell essentiality labeled

expr_graph2<-ggplot(Sum,aes(x=reorder(Sum$Gene,-Sum$percent),y=Sum$percent,fill=factor(Sum$essential)))+
  geom_bar(stat="identity",width=.6,alpha=0.9)+
  xlab("Genes")+ylab("Expression Percentage (%)")+
  scale_y_continuous(expand=c(0,0), limits = c(0, 105))+
  scale_fill_manual(values=c("#999999","brown2"),labels=c("non-essential","essential"))+
  theme(axis.text=element_text(size=10, color="black"),axis.text.x=element_text(angle=30,vjust=0.9,hjust=0.8),
        legend.title=element_blank(),legend.text=element_text(size=10, color="black"),
        axis.title=element_text(size=12,face="bold"), panel.background = element_blank(), 
        panel.border = element_rect(color="black", fill=NA, size=0.8))
expr_graph2

# png("gene_expr_with_CS.png",width=930,height=550)
# expr_graph2
# dev.off()


# cons_graph<-ggplot(Sum,aes(x=reorder(Sum$Gene,-Sum$percent),y=Sum$mean_cons_score,group=1))+
#   geom_line(alpha=0.7)+geom_point(aes(color=factor(Sum$essential)))+
#   scale_color_manual(values=c("#999999","brown2"),labels=c("non-essential","essential"))+
#   scale_y_continuous(breaks=seq(floor(min(Sum$mean_cons_score)),max(Sum$mean_cons_score),0.05))+
#   theme(axis.text.x=element_blank(),axis.text.y=element_text(size=14),
#         axis.title.x=element_blank(),axis.title.y=element_text(size=15,face="bold"),
#         axis.ticks.x=element_blank(),panel.grid.minor=element_blank(),
#         legend.title=element_blank(),legend.text=element_text(size=11))+
#   ylab("phastCons score")
# cons_graph
# 
# 
# layout<-matrix(c(rep(1,3),rep(2,6)),nrow=3,byrow=T)
# multiplot(cons_graph,expr_graph2,layout=layout)
# 
# 
# png("HMG_expr_with_cons.png",width=930,height=750)
# multiplot(cons_graph,expr_graph2,layout=layout)
# dev.off()
# 
# 
# 




