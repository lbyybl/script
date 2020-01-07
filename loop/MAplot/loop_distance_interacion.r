#--- find the loop interaction distribution along the distance
#--- 1. 读取wt和ko，2. 看距离分布画成freq图
oceanc_wt <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/Contact_strength/hiccups/E_E_int_distribution/pol2_ocean_c_wt_sub_loop.bed'
oceanc_ko <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/Contact_strength/hiccups/E_E_int_distribution/pol2_ocean_c_ko_sub_loop.bed'

hichip_wt <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/Contact_strength/hiccups/E_E_int_distributinon/pol2_hichip_wt_sub_loop.bed'
hichip_ko <- '/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/Contact_strength/hiccups/E_E_int_distributinon/pol2_hichip_ko_sub_loop.bed'
library(dplyr)
library(data.table)
library(ggplot2)
readdata <- function(file){
  data <- fread(file, col.names = c('chr1','st1','st2','pets'))
  data$distance <- abs(data$st2 - data$st1)
  return(data)
}
# readdata <- function(wt,ko){
#   wt_file <- fread(wt, col.names = c('chr1','st1','st2','pets'))
#   wt_file$class <- 'wt'
#   ko_file <- fread(ko, col.names = c('chr1','st1','st2','pets'))
#   ko_file$class <- 'ko'
#   data <- bind_rows(wt_file, ko_file)
#   return(data)
# }

oceanc_wt_file <- readdata(oceanc_wt)
oceanc_ko_file <- readdata(oceanc_ko)
hichip_wt_file <- readdata(hichip_wt)
hichip_ko_file <- readdata(hichip_ko)
distance_distribute <- function(data){
  data <- data %>%
    group_by(distance) %>%
    summarise(int=mean(pets))
}

oceanc_wt_stas <- distance_distribute(oceanc_wt_file)
oceanc_ko_stas <- distance_distribute(oceanc_ko_file)
hichip_wt_stas <- distance_distribute(hichip_wt_file)
hichip_ko_stas <- distance_distribute(hichip_ko_file)

# ggplot(oceanc_wt_stas,aes(distance,int)) + geom_line(color='red') +
#   geom_line(data=oceanc_ko_stas,aes(distance,int),color='blue')
# 
# ggplot(hichip_wt_stas,aes(distance,int)) + geom_line(color='red') +
#   geom_line(data=hichip_ko_stas,aes(distance,int),color='blue')

hichip_data <- merge(hichip_wt_file,hichip_ko_file,by=c('chr1','st1','st2'))
colnames(hichip_data) <- c('chr1','st1','st2','pets','distance','pets_2','distance2')
hichip_data$distance2 <- "duan"
hichip_data[which(hichip_data$distance>=300000),]$distance2 <- "long"
hichip_data[which(hichip_data$distance<300000),]$distance2 <- "short"

oceanc_data <- merge(oceanc_wt_file,oceanc_ko_file,by=c('chr1','st1','st2'))
colnames(oceanc_data) <- c('chr1','st1','st2','pets','distance','pets_2','distance2')

ggplot(hichip_data,aes(x=pets, y=pets_2, color=log10(distance)
)) +
  scale_color_gradientn(colors = c("Wheat","Khaki","Gold","Orange","Coral","Salmon","OrangeRed","Red", "MediumVioletRed","DeepPink"),
                        values = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)#,
                        #limits=c(-750,300)
  ) +
  geom_point(shape=15) +
  geom_abline(intercept = 0, slope = 1,linetype="dashed", size=0.2) +
  theme(plot.title=element_text(hjust=0.5)) +
  labs(x="wt loop PETs",y="ko loop PETs",color="") +
  # scale_y_continuous(expand=c(0,0),limits=c(v_min,v_max)) +
  # scale_x_continuous(expand=c(0,0),limits=c(v_min,v_max))+
  scale_y_continuous(expand=c(0,0),limits=c(0,220)) +
  scale_x_continuous(expand=c(0,0),limits=c(0,220))+
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

ggplot(oceanc_data,aes(x=pets, y=pets_2, color=log10(distance)
)) +
  scale_color_gradientn(colors = c("Wheat","Khaki","Gold","Orange","Coral","Salmon","OrangeRed","Red", "MediumVioletRed","DeepPink"),
                        values = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)#,
                        #limits=c(-750,300)
  ) +
  geom_point(shape=15) +
  geom_abline(intercept = 0, slope = 1,linetype="dashed", size=0.2) +
  theme(plot.title=element_text(hjust=0.5)) +
  labs(x="wt loop PETs",y="ko loop PETs",color="") +
  # scale_y_continuous(expand=c(0,0),limits=c(v_min,v_max)) +
  # scale_x_continuous(expand=c(0,0),limits=c(v_min,v_max))+
  scale_y_continuous(expand=c(0,0),limits=c(0,90)) +
  scale_x_continuous(expand=c(0,0),limits=c(0,90))+
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


