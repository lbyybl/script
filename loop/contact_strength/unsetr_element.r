# this script is used to stastic the overlap between ocean and hichip
# for hiccups ; and due to the big intervel for loop anchor, many elements
# share the same loop; so this script will also stastic the shared loop 
# overlap

library(UpSetR)
library(stringr)
library(dplyr)
#--- for elements shared loops
setwd('/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/Contact_strength/hiccups')

got_file_list <- function(prefix){
  unset_list<-list()
  i=1
  for (file in system('ls *.bed',intern = T)){
    #print(file)
    read_data <- read.table(file, col.names = c('chr1','st1','en1','chr2','st2','en2','d1','d2'))
    read_data <- read_data %>%
      mutate(loop=str_c(chr1,st1,en1,chr2,st2,en2,sep = "-"))
    assign(paste0(prefix,strsplit(file,"\\.")[[1]][1]), read_data$loop)
    unset_list[[i]] <- get(paste0(prefix,strsplit(file,"\\.")[[1]][1]))
    names(unset_list)[i] <- paste0(prefix,strsplit(file,"\\.")[[1]][1])
    i=i+1
  }
  return(unset_list)
}
unset_list <- got_file_list('')
jpeg('hichip.jpeg',width = 1000, height = 1000)
upset(fromList(unset_list),nsets=9, order.by = "freq", mb.ratio = c(0.3, 0.7))
dev.off()
#-- oceanc elements overlap
setwd('/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/Contact_strength/hiccups')
unset_list<-list()
i=1
for (file in system('ls *.bed',intern = T)){
  #print(file)
  read_data <- read.table(file, col.names = c('chr1','st1','en1','chr2','st2','en2','d1','d2'))
  read_data <- read_data %>%
    mutate(loop=str_c(chr1,st1,en1,chr2,st2,en2,sep = "-"))
  assign(strsplit(file,"\\.")[[1]][1], read_data$loop)
  unset_list[[i]] <- get(strsplit(file,"\\.")[[1]][1])
  names(unset_list)[i] <- strsplit(file,"\\.")[[1]][1]
  i=i+1
}
jpeg('oceanc.jpeg',width = 1000, height = 1000)
upset(fromList(unset_list),nsets=, order.by = "freq", mb.ratio = c(0.3, 0.7))
dev.off()

#--- hichip and oceanc overlap

setwd('/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/Contact_strength/hiccups')
unset_list<-list()
i=1
for (file in system('ls *.bed',intern = T)){
  #print(file)
  read_data <- read.table(file, col.names = c('chr1','st1','en1','chr2','st2','en2','d1','d2'))
  read_data <- read_data %>%
    mutate(loop=str_c(chr1,st1,en1,chr2,st2,en2,sep = "-"))
  assign(strsplit(file,"\\.")[[1]][1], read_data$loop)
  unset_list[[i]] <- get(strsplit(file,"\\.")[[1]][1])
  names(unset_list)[i] <- paste0('hichip',strsplit(file,"\\.")[[1]][1])
  i=i+1
}
#-- oceanc 
setwd('/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/Contact_strength/hiccups')

for (file in system('ls *.bed',intern = T)){
  #print(file)
  read_data <- read.table(file, col.names = c('chr1','st1','en1','chr2','st2','en2','d1','d2'))
  read_data <- read_data %>%
    mutate(loop=str_c(chr1,st1,en1,chr2,st2,en2,sep = "-"))
  assign(strsplit(file,"\\.")[[1]][1], read_data$loop)
  unset_list[[i]] <- get(strsplit(file,"\\.")[[1]][1])
  names(unset_list)[i] <- paste0('oceanc',strsplit(file,"\\.")[[1]][1])
  i=i+1
}
jpeg('hichip_oceanc.jpeg',width = 1000, height = 1000)
upset(fromList(unset_list),nsets=18, order.by = "freq", mb.ratio = c(0.3, 0.7))
dev.off()
#--- hichip and oceanc overlap total
setwd('/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_h3k27ac_hichip_loop/allvalidpairs/Contact_strength/hiccups')
unset_list<-list()
i=1
for (file in system('ls ko_hichp.bedpe wt_hichp.bedpe',intern = T)){
  #print(file)
  read_data <- read.table(file, col.names = c('chr1','st1','en1','chr2','st2','en2','d1','d2'))
  read_data <- read_data %>%
    mutate(loop=str_c(chr1,st1,en1,chr2,st2,en2,sep = "-"))
  assign(strsplit(file,"\\.")[[1]][1], read_data$loop)
  unset_list[[i]] <- get(strsplit(file,"\\.")[[1]][1])
  names(unset_list)[i] <- paste0('hichip',strsplit(file,"\\.")[[1]][1])
  i=i+1
}
#-- oceanc 
setwd('/DATA/work/lbyybl/ypjiang/geneloop/hichip_oceanc/pol2_ocean_c/allvalidpairs/Contact_strength/hiccups')

for (file in system('ls ko_oceanc.bedpe wt_oceanc.bedpe',intern = T)){
  #print(file)
  read_data <- read.table(file, col.names = c('chr1','st1','en1','chr2','st2','en2','d1','d2'))
  read_data <- read_data %>%
    mutate(loop=str_c(chr1,st1,en1,chr2,st2,en2,sep = "-"))
  assign(strsplit(file,"\\.")[[1]][1], read_data$loop)
  unset_list[[i]] <- get(strsplit(file,"\\.")[[1]][1])
  names(unset_list)[i] <- paste0('oceanc',strsplit(file,"\\.")[[1]][1])
  i=i+1
}
jpeg('hichip_oceanc.jpeg',width = 1000, height = 1000)
upset(fromList(unset_list),nsets=4, order.by = "freq", mb.ratio = c(0.3, 0.7))
dev.off()




