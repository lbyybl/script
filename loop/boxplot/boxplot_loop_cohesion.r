#---boxplot for with cohesion and no cohesion loop
#setwd('/DATA2/work/lbyybl/cohesion/peak')
pol1 <- '/DATA2/work/lbyybl/cohesion/peak/pol1_hic'
pol2 <- '/DATA2/work/lbyybl/cohesion/peak/pol2_hic'
pol3 <- '/DATA2/work/lbyybl/cohesion/peak/pol3_hic'
library(ggplot2)
library(data.table)
library(dplyr)
library(ggpubr)
sample_list <- c('no_cohe','tow_cohe','one_cohe')

readfile <- function(name){
  data <- fread(name,col.names = c('loc','interaction'))
  return(data)
}

readpathfile <- function(prefix,path){
  setwd(path)
  all_data <- data.frame('loc'=NULL,
                         'interaction'=NULL,
                         'class'=NULL,
                         'sample'=NULL)
  for (i in sample_list){
    assign(paste0(prefix,i,'_wt'),readfile(paste0(prefix,i,'_wt_loop.bed')))
    data <- get(paste0(prefix,i,'_wt'))
    data$class <- 'wt'
    data$sample <- i
    all_data <- bind_rows(all_data,data)
    assign(paste0(prefix,i,'_wt'),data)
    assign(paste0(prefix,i,'_ko'),readfile(paste0(prefix,i,'_ko_loop.bed')))
    data <- get(paste0(prefix,i,'_ko'))
    data$class <- 'ko'
    data$sample <- i
    assign(paste0(prefix,i,'_ko'),data)
    all_data <- bind_rows(all_data,data)
    
  }
  return(all_data)
}
pol1_data <- readpathfile(prefix = 'pol1_',path = pol1)
pol2_data <- readpathfile(prefix = 'pol2_',path = pol2)
pol3_data <- readpathfile(prefix = 'pol3_',path = pol3)

drawboxplot <- function(all_data){
  ggplot(all_data,aes(x=sample,y=interaction,fill=class))+geom_boxplot(outlier.color = NA)+
    scale_y_continuous(limits = c(0,50)) +
    theme(plot.title=element_text(hjust=0.5),
          panel.grid=element_blank(), panel.background=element_blank()) +
    theme(
      title = element_blank() ,
      legend.key.width = unit(2,"line"),
      legend.key.height = unit(0.5,"cm"),
      legend.text = element_text(
        size = 7,
        hjust = 1.2,
        face = "italic",
        colour = "black"
      ),
      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
      plot.title=element_text(colour="black", size=10),
      axis.text.x=element_text(#face = "bold",
        colour="black", size=12),
      axis.text.y=element_text(#face = "bold",
        colour="black", size=12),
      panel.background=element_blank(),
      plot.background=element_blank()) +
    labs(fill=NULL )
}

drawboxplot(pol1_data)
ggsave('pol1_interaction_change_and_cohesion.pdf',width = 5,height = 5)
drawboxplot(pol2_data)
ggsave('pol2_interaction_change_and_cohesion.pdf',width = 5,height = 5)
drawboxplot(pol3_data)
ggsave('pol3_interaction_change_and_cohesion.pdf',width = 5,height = 5)


#---- t-test for three data
# 1. no

dotest <- function(data){
  data_no_wt <- data %>% filter(class == 'wt' & sample == 'no_cohe')
  data_one_wt <- data %>% filter(class == 'wt' & sample == 'one_cohe')
  data_both_wt <- data %>% filter(class == 'wt' & sample == 'tow_cohe')
  data_no_ko <- data %>% filter(class == 'ko' & sample == 'no_cohe')
  data_one_ko <- data %>% filter(class == 'ko' & sample == 'one_cohe')
  data_both_ko <- data %>% filter(class == 'ko' & sample == 'tow_cohe')
  
  no_pvalue <- t.test(data_no_wt$interaction,data_no_ko$interaction)$p.value
  one_pvalue <- t.test(data_one_wt$interaction,data_one_ko$interaction)$p.value
  two_pvalue <- t.test(data_both_wt$interaction,data_both_ko$interaction)$p.value
  print(paste0('no_P => ',no_pvalue))
  print(paste0('one_P => ',one_pvalue))
  print(paste0('tow_P => ',two_pvalue))
}
dotest(pol1_data)
dotest(pol2_data)
dotest(pol3_data)
#my_comparisons <- list(c("no_cohe","one_cohe"),c("one_cohe","two_cohe"),c("no_cohe","two_cohe"))

