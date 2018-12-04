# 2018/11/7
# Boyuan_li
# it's used to trans the file output from find_region.sh
command=matrix(c("Input_file","f",1,"character",

                                "Output_file","o",1,"character",

                 "help","h",0,"logical"),byrow=T,ncol=4)

args=getopt::getopt(command)

if (!is.null(args$help) || is.null(args$Input_file) || is.null(args$Output_file)) {

    cat(paste(getopt::getopt(command, usage = T), "\n"))

    q()

}


library(data.table)
library(dplyr)
library(stringr)
interafile <- fread(args$Input_file)
#View(interafile)
interafile <- interafile[,c(1:6,23,28)]
colnames(interafile) <- c("chr1","start1","end1","chr2","start2","end2","region1","region2")
interafile$region1 <- str_sub(interafile$region1,1,4)
interafile$region2 <- str_sub(interafile$region2,1,4)
interafile$region2[interafile$region2 != "NONG"] <- "GERE"
interafile$region1[interafile$region1 != "NONG"] <- "GERE"
interafile <- interafile %>%
  mutate(regiont = paste0(chr1,"-",start1,"-",end1,"-",chr2,"-",start2,"-",end2),interaction=paste0(region1,"-",region2))
fwrite(interafile,file = args$Output_file, sep = "\t")
