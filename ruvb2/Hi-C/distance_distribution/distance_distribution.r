setwd('/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/loop/distance_distribute')
wt_file <- '/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/hic_file/RUVB2_hic_merge_wt_sub.allValidPairs'
ko_file <- '/DATA/work/lbyybl/wh/ruvb2/Hi-C/sample0227/hic_file/RUVB2_hic_merge_ko_sub.allValidPairs'

wt_data <- fread(wt_file)
wt_data <- wt_data[,2:6]
colnames(wt_data) <- c('chr1','loc1','sd1','chr2','loc2')
ko_data <- fread(ko_file)
ko_data <- ko_data[,2:6]
colnames(ko_data) <- c('chr1','loc1','sd1','chr2','loc2')
chr_filter <- function(data,chr){
  chr_data <- data %>%
    filter(chr1==chr)
  chr_data <- chr_data %>%
    mutate(dis=abs(loc2-loc1)/10^6)
  return(chr_data)
}
wt_chr1 <- chr_filter(wt_data,'chr2')
ko_chr1 <- chr_filter(ko_data,'chr2')
ggplot(data=ko_chr1, aes(x=dis))+
  geom_line(stat = "density", adjust=1.5, size=0.7)+
  scale_x_log10()+scale_y_log10()+
  xlab("Distance (Kb)")+
  scale_color_manual(values=c("black", "darkgreen", "tomato", "dodgerblue3"))+
  theme(element_blank(), panel.background = element_blank(),
        panel.border = element_rect(size=0.8, fill=NA),
        axis.text = element_text(color="black"),
        axis.title = element_text(color="black"))