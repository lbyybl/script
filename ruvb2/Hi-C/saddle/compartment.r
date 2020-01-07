#--- change the ob/ex to dense and change store it in the raster
setwd('/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/merge/saddle_plot/ko_plot/_ploter/500000')
setwd('/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/merge/saddle_plot/wt_plot/_ploter/500000')
library(tidyverse)
library(raster)
#-- read the ob/ex file
for (i in 1:19){
  file_name <- paste0("chr",i,".oe.resort")
  assign(paste0("chr",i,"_oe"),read_csv(file_name, col_names = T))
  data <- get(paste0("chr", i, "_oe"))
  colnames(data) <- c("bin1", "bin2", "interaction")
  #data <- data[-which(data$bin1==data$bin2),,] # remove self interaction
  assign(paste0("chr", i, "_oe"),data)
}

date()
#--- rm the line with non-zero value less than 20
remove_sparse_line <- function(data){
  data_rm_sparse <- data %>%
    group_by(bin1) %>%
    filter(sum(interaction!=0)>20) 
  bin_sect <- intersect(unique(data_rm_sparse$bin2), unique(data_rm_sparse$bin1))
  data_rm_sparse <- data_rm_sparse %>%
    filter(bin2 %in% bin_sect)
  return(data_rm_sparse)
}

for (i in 1:19){
  data <- get(paste0("chr",i,"_oe"))
  data_rm_sparse <- data %>%
    group_by(bin1) %>%
    filter(sum(interaction!=0)>20) 
  bin_sect <- intersect(unique(data_rm_sparse$bin2), unique(data_rm_sparse$bin1))
  data_rm_sparse <- data_rm_sparse %>%
    filter(bin2 %in% bin_sect)
  assign(paste0("chr",i,"_dense_rm_zero"),data_rm_sparse)
}

date()
#-- rm the blank line

remove_blank_line <- function(data){
  bin_min <- min(data$bin1,data$bin2)
  bin_max <- max(data$bin1,data$bin2)
  data <- data[order(data$bin1,data$bin2),,]
  bin_number <- length(unique(data$bin1))
  data_rm_blank <- data.frame(bin1=rep(1:bin_number, each=bin_number), 
                              bin2=rep(1:bin_number, bin_number))
  data <- data[order(data$bin1,data$bin2),,]
  data_rm_blank$interaction <- data$interaction
  return(data_rm_blank)
}

for (i in 1:19){
  data <- get(paste0("chr", i, "_dense_rm_zero"))
  assign(paste0("chr",i,"_dense_rm_blank"), remove_blank_line(data))
}

date()
transform_resolution <- function(dense, res){
  r <- raster(nrow=sqrt(nrow(dense)), ncol=sqrt(nrow(dense)))
  dense <- dense[order(dense$bin1,dense$bin2),,]
  r[] <- dense$interaction
  s <- raster(nrow=res, ncol=res)
  s <- resample(r, s, method="bilinear")
  newdense <- data.frame(bin1=rep(1:res, each=res), bin2=rep(1:res, res), interaction=s[])
  return(newdense)
}

#--- change the resolution to 80 or 122

reso_res <- 300

for (i in 1:19){
  data <- get(paste0("chr",i,"_dense_rm_blank"))
  assign(paste0("chr",i,"_cR_oe"),transform_resolution(data,reso_res))
}
date()

#-----get expect

data_oe <- data.frame(bin1=chr1_cR_oe$bin1, bin2=chr1_cR_oe$bin2,
                      interaction=c(rep(0,nrow(chr1_cR_oe))))
for (i in 1:19){
  data <- get(paste0("chr",i,"_cR_oe"))
  data_oe$interaction <- data_oe$interaction + data$interaction
}
data_oe$interaction <- data_oe$interaction/19
drawheatmap <- function(data){
  ggplot(data=data, aes(x=bin1, y=bin2, fill = log(interaction+1))) + geom_raster() +
    #scale_fill_gradient2(low="dodgerblue2", mid="white", high="orangered2",midpoint=2) +
    scale_fill_gradientn(colors = c("navy","royalblue4", "white", "orangered", "orangered2"),
                         values = c(0,0.2,0.3,0.4,1), limits=c(0.25,1.75)
    ) +
    # scale_fill_manual(breaks=c("0","0.5","1.5","7.5"),
    #                   values=c("dodgerblue", "white", "orange", "orangered")) +
    scale_y_reverse(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    theme(plot.title=element_text(hjust=0.5),
          #panel.border = element_rect(colour="black",  fill=NA, size=1),
          panel.grid=element_blank(), panel.background=element_blank()) +
    xlab("bin 25k") +
    ylab("bin 25k") + theme(axis.line=element_blank(),
                            axis.text.x=element_blank(),
                            axis.text.y=element_blank(),
                            axis.ticks=element_blank(),
                            # axis.title.x=element_blank(),
                            # axis.title.y=element_blank(),
                            #legend.position="none",
                            panel.background=element_blank(),
                            panel.border=element_blank(),
                            panel.grid.major=element_blank(),
                            panel.grid.minor=element_blank(),
                            plot.background=element_blank())
}
drawheatmap(chr6_cR_oe)
drawheatmap(chr6_oe)
drawheatmap(chr6_dense_rm_zero)
drawheatmap(data_oe)
setwd('/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/merge/saddle_plot/')
ggsave('R2_ko_500k_saddle.pdf',width = 7.5, height = 6)
wt_data_oe <- data_oe
ko_data_oe <- data_oe
diff_data_oe <- ko_data_oe
diff_data_oe$interaction <- ko_data_oe$interaction-wt_data_oe$interaction
draw_diff <- function(data){
  ggplot(data=data, aes(x=bin1, y=bin2, fill = interaction)) + geom_raster() +
    #scale_fill_gradient2(low="dodgerblue2", mid="white", high="orangered2",midpoint=2) +
    scale_fill_gradientn(colors = c("navy","royalblue4", "white", "orangered", "orangered2"),
                         values = c(0,0.45,0.5,0.65,1),limits=c(-3.3,3.3)
    ) +
    # scale_fill_manual(breaks=c("0","0.5","1.5","7.5"),
    #                   values=c("dodgerblue", "white", "orange", "orangered")) +
    scale_y_reverse(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    theme(plot.title=element_text(hjust=0.5),
          #panel.border = element_rect(colour="black",  fill=NA, size=1),
          panel.grid=element_blank(), panel.background=element_blank()) +
    xlab("bin 25k") +
    ylab("bin 25k") + theme(axis.line=element_blank(),
                            axis.text.x=element_blank(),
                            axis.text.y=element_blank(),
                            axis.ticks=element_blank(),
                            # axis.title.x=element_blank(),
                            # axis.title.y=element_blank(),
                            #legend.position="none",
                            panel.background=element_blank(),
                            panel.border=element_blank(),
                            panel.grid.major=element_blank(),
                            panel.grid.minor=element_blank(),
                            plot.background=element_blank())
}
ggsave('R2_diff_500k_saddle.pdf',width = 7.5, height = 6)
#--- compare the difference of two method and why difference
#--- for eg1
#--- for old
setwd("/DATA/work/lbyybl/wangcl/hic/KPNA2/trim_linker/DoStreAna/Compartment/Resulation500k/ploter/kpna2_saddle/erro_detect/zscore/old")
old_chr6zscorefile <- read_tsv("sample1018_ko1_500000_iced_chr6_dense.addedHeaders.zScore.compartments")
View(old_chr6zscorefile)
nrow(old_chr6zscorefile) # 300
old_chr6zscorefile <- select(old_chr6zscorefile, start, end, eigen1,index)
old_chr6zscorefile$eigen1[is.na(old_chr6zscorefile$eigen1)] <- 0 
ggplot(old_chr6zscorefile,aes(index, eigen1)) + geom_line(color="blue")

#--- for new
setwd("/DATA/work/lbyybl/wangcl/hic/KPNA2/trim_linker/DoStreAna/Compartment/Resulation500k/ploter/kpna2_saddle/erro_detect/zscore/observed")
new_chr6zscorefile <- read_tsv("chr6_sample_all_ko_merge_observed_dense.addedHeaders.zScore.compartments")  
View(new_chr6zscorefile) 
nrow(new_chr6zscorefile) # 300
new_chr6zscorefile <- select(new_chr6zscorefile,start, end, eigen1,index)
new_chr6zscorefile$eigen1[is.na(new_chr6zscorefile$eigen1)] <- 0
ggplot(new_chr6zscorefile,aes(index,eigen1)) + geom_line(color="red")
#--- draw the eg1 for two file
ggplot(old_chr6zscorefile,aes(index, eigen1)) + geom_line(color="blue") + 
  geom_line(data=new_chr6zscorefile,aes(index,eigen1),color="red")
#--- combine towfiles eg1
combine_eg1_file <- merge(old_chr6zscorefile,new_chr6zscorefile,by=c("start","end"),all=T)
View(combine_eg1_file)
colnames(combine_eg1_file)<-c("start","end","eigen_old","index_old","eigen_new","index_new")
ggplot(combine_eg1_file,aes(index_old,eigen_old)) + geom_line(color="blue") +
  geom_line(aes(index_old,eigen_new_mod),color="red")
combine_eg1_file$eigen_new_mod <- combine_eg1_file$eigen_new
cor(combine_eg1_file$eigen_old,combine_eg1_file$eigen_new_mod) # 0.9896837
cor(combine_eg1_file$eigen_old,combine_eg1_file$eigen_new) # 0.833433
#--- for sort file
#--- for the corresponding number file old
setwd("/DATA/work/lbyybl/wangcl/hic/KPNA2/trim_linker/DoStreAna/Compartment/Resulation500k/ploter/kpna2_saddle/erro_detect/ploter/old")
old_cornumber <- read_tsv("chr6_correspondingnumber_",col_names = F)
colnames(old_cornumber) <- c("old_bin","new_bin")
View(old_cornumber)
ggplot(old_cornumber,aes(new_bin,old_bin)) + geom_point(color="blue")
#-- for new
setwd("/DATA/work/lbyybl/wangcl/hic/KPNA2/trim_linker/DoStreAna/Compartment/Resulation500k/ploter/kpna2_saddle/erro_detect/ploter/observed")
new_cornumber <- read_tsv("chr6_correspondingnumber_",col_names = F)
colnames(new_cornumber) <- c("old_bin","new_bin")
View(new_cornumber)
#--
ggplot(old_cornumber,aes(old_bin,new_bin)) + geom_line(color="blue") + 
  geom_line(data=new_cornumber, aes(old_bin+1,new_bin),color="red") 

#--- try a new sor method

#--- read the chr6 file
setwd("/DATA/work/lbyybl/wangcl/hic/KPNA2/trim_linker/DoStreAna/Compartment/Resulation500k/ploter/kpna2_saddle/erro_detect/ploter/observed")
chr6_test <- read_tsv("chr6_ko.oe",col_names = F)
colnames(chr6_test) <- c("bin1","bin2","interaction")
View(chr6_test)
chr6_test_bin <- chr6_test %>%
  mutate(bin1=bin1/500000 +1 ,bin2=bin2/500000 +1)

egivector_vec <- combine_eg1_file$eigen_new

chr6_test_bin_dense <- sparse2dense(chr6_test_bin,1,300)
chr6_observed_bin_dense <- chr6_test_bin_dense
oe_ob <- merge(chr6_observed_bin_dense,chr6_test_bin_dense,by=c("bin1","bin2"),all=T)
is.na(oe_ob$interaction.x)
oe_ob <- mutate(oe_ob, expected=interaction.x/interaction.y)
oe_ob$expected[is.na(oe_ob$expected)] <- 0
oe_ob <- select(oe_ob, bin1, bin2, expected)
colnames(oe_ob) <- c("bin1","bin2","interaction")
drawheatmap(chr6_test_bin_dense)
chr6_test_bin_egi <- chr6_test_bin_dense #&>&
#mutate(bin1_eign=egivector_vec[chr6_test_bin_dense$bin1], bin2_eign=egivector_vec[chr6_test_bin_dense$bin2])
chr6_test_bin_egi$bin1_egi <- egivector_vec[chr6_test_bin_egi$bin1]
chr6_test_bin_egi$bin2_egi <- egivector_vec[chr6_test_bin_egi$bin2]
chr6_test_bin_egi_sort <- chr6_test_bin_egi[order(chr6_test_bin_egi$bin1_egi,chr6_test_bin_egi$bin2_egi),,]
View(chr6_test_bin_egi_sort)
new_matrix_sort <- data.frame(bin1=c(rep(1:300,each=300)),
                              bin2=c(rep(1:300,300))) 
new_matrix_sort$interaction <- chr6_test_bin_egi_sort$interaction
new_matrix_sort_rm_zero <- remove_sparse_line(new_matrix_sort)
new_matrix_sort_rm_blank <- remove_blank_line(new_matrix_sort_rm_zero)
drawheatmap(new_matrix_sort_rm_blank)
