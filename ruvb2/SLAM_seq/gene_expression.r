setwd('/DATA/work/lbyybl/wh/ruvb2/SLAM/slam_output/count_2TC')

read_file <- function(read1,read2){
  R1_file <- fread(read1,skip = 2)
  R1_file_sub <- R1_file %>%
    select(Name, ReadCount,TcReadCount) 
  R2_file <- fread(read2,skip = 2)
  R2_file_sub <- R2_file %>%
    select(Name, ReadCount,TcReadCount) 
  
  merge_file <- merge(R1_file_sub,R2_file_sub,by='Name')
  
  merge_file <- merge_file %>%
    mutate(total=ReadCount.x+ReadCount.y,new=TcReadCount.x+TcReadCount.y) %>%
    select(Name,total,new) %>%
    mutate(rate=new/total)
  
}

WT_R1 <- 'L10-R2-NT+IAA-SLAM_FKDL171549250-1A_1.fq_slamdunk_mapped_filtered_tcount.tsv'
WT_R2 <- 'L10-R2-NT+IAA-SLAM_FKDL171549250-1A_2.fq_slamdunk_mapped_filtered_tcount.tsv'

WT_slam <- read_file(WT_R1,WT_R2)

KO_R1 <- 'L12-R2-IAA+IAA-SLAM_FKDL171549252-1A_1.fq_slamdunk_mapped_filtered_tcount.tsv'
KO_R2 <- 'L12-R2-IAA+IAA-SLAM_FKDL171549252-1A_2.fq_slamdunk_mapped_filtered_tcount.tsv'

KO_slam <- read_file(KO_R1,KO_R2)

merge_data <- merge(WT_slam,KO_slam,by='Name')
colnames(merge_data) <- c('Name','wt_total','wt_new','wt_rate','ko_total','ko_new','ko_rate')
merge_data <- merge_data %>%
  filter(wt_rate + ko_rate>0) %>%
  mutate(wtvsko=wt_rate/ko_rate,wtplusko=wt_rate+ko_rate)

ggplot(merge_data,aes(x=log(wt_rate),y=log(ko_rate),color=wtplusko)) + geom_point(alpha=0.1,shape=15) +
  #scale_color_gradientn(colors = c("Wheat","Khaki","Gold","Orange","Coral","Salmon","OrangeRed","Red", "MediumVioletRed","DeepPink")) +
  scale_color_gradientn(colors = 'red','red')+
  geom_abline(intercept = 0, slope = 1,linetype="dashed", size=0.2) +
  theme(plot.title=element_text(hjust=0.5)) +
  labs(x="wt APA",y="ko APA",color="") +
  # scale_y_continuous(expand=c(0,0),limits=c(v_min,v_max)) +
  # scale_x_continuous(expand=c(0,0),limits=c(v_min,v_max))+
  scale_y_continuous(expand=c(0,0),limits = c(-6,0)) +
  scale_x_continuous(expand=c(0,0),limits=c(-6,0))+
  # scale_color_continuous(limits=c(-750,300))
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
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.text.x=element_text(angle=60, hjust=0,vjust = 0,
                             colour="black", size=12),
    axis.text.y=element_text(#face = "bold",
      colour="black", size=12),
    panel.background=element_blank()
  )
ggsave('SLAM_seq.pdf',width = 5.8,height = 5)
down <- sum(merge_data$wt_rate>merge_data$ko_rate)
up <- sum(merge_data$wt_rate<merge_data$ko_rate)

gene_exp <- data.frame('class'=c('up','down'),
                       'value'=c(up,down))

ggplot(gene_exp,aes(x=class,y=value)) + geom_col(fill='blue') +
  #scale_color_gradientn(colors = c("Wheat","Khaki","Gold","Orange","Coral","Salmon","OrangeRed","Red", "MediumVioletRed","DeepPink")) +
  #scale_color_gradientn(colors = 'red','red')+
  geom_abline(intercept = 0, slope = 1,linetype="dashed", size=0.2) +
  theme(plot.title=element_text(hjust=0.5)) +
  labs(x="wt APA",y="ko APA",color="") +
  # scale_y_continuous(expand=c(0,0),limits=c(v_min,v_max)) +
  # scale_x_continuous(expand=c(0,0),limits=c(v_min,v_max))+
  #scale_y_continuous(expand=c(0,0),limits = c(-6,0)) +
  #scale_x_continuous(expand=c(0,0),limits=c(-6,0))+
  # scale_color_continuous(limits=c(-750,300))
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
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.text.x=element_text(angle=60, hjust=0,vjust = 0,
                             colour="black", size=12),
    axis.text.y=element_text(#face = "bold",
      colour="black", size=12),
    panel.background=element_blank()
  )
ggsave('SLAM_gene_exp.pdf',width = 3,height = 5)
#ggplot(merge_data,aes(x=wtplusko,y=log(wtvsko))) + geom_point(alpha=0.1)

wc_r1 <- 'L11-R2-NT-IAA-SLAM_FKDL171549251-1A_1.fq_slamdunk_mapped_filtered_tcount.tsv'
wc_r2 <- 'L11-R2-NT-IAA-SLAM_FKDL171549251-1A_2.fq_slamdunk_mapped_filtered_tcount.tsv'
WT_control <- read_file(wc_r1,wc_r2)
kc_r1 <- 'L13-R2-IAA-IAA-SLAM_FKDL171549253-1A_1.fq_slamdunk_mapped_filtered_tcount.tsv'
kc_r2 <- 'L13-R2-IAA-IAA-SLAM_FKDL171549253-1A_2.fq_slamdunk_mapped_filtered_tcount.tsv'
KO_control <- read_file(kc_r1,kc_r2)

quality_data <- data.frame('wt_control'=sum(WT_control$new),
                           'wt_slam'=sum(WT_slam$new),
                           'ko_control'=sum(KO_control$new),
                           'ko_slam'=sum(KO_slam$new))
quality_data <- data.frame('class'=c('wt_control','wt_slam','ko_control','ko_slam'),
                           'value'=c(sum(WT_control$new)/sum(WT_control$total),
                                     sum(WT_slam$new)/sum(WT_slam$total),
                                     sum(KO_control$new)/sum(KO_control$total),
                                     sum(KO_slam$new)/sum(KO_slam$total)))
ggplot(quality_data,aes(x=class,y=value)) + geom_col() +
  #scale_color_gradientn(colors = c("Wheat","Khaki","Gold","Orange","Coral","Salmon","OrangeRed","Red", "MediumVioletRed","DeepPink")) +
  #scale_color_gradientn(colors = 'red','red')+
  geom_abline(intercept = 0, slope = 1,linetype="dashed", size=0.2) +
  theme(plot.title=element_text(hjust=0.5)) +
  labs(x="wt APA",y="ko APA",color="") +
  # scale_y_continuous(expand=c(0,0),limits=c(v_min,v_max)) +
  # scale_x_continuous(expand=c(0,0),limits=c(v_min,v_max))+
  #scale_y_continuous(expand=c(0,0),limits = c(-6,0)) +
  #scale_x_continuous(expand=c(0,0),limits=c(-6,0))+
  # scale_color_continuous(limits=c(-750,300))
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
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.text.x=element_text(angle=60, hjust=0,vjust = 0,
                             colour="black", size=12),
    axis.text.y=element_text(#face = "bold",
      colour="black", size=12),
    panel.background=element_blank()
  )
ggsave('SLAM_quality_ratio.pdf',width = 5,height = 5)