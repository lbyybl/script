#!/usr/bin/R
# 2018-9-7
# this is used to sum the interaction, only for these data

command=matrix(c("Input","i",1,"character",
				 "Output","o",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)

args=getopt::getopt(command)

if (!is.null(args$help) || is.null(args$Input)  || is.null(args$Output)) {

    cat(paste(getopt::getopt(command, usage = T), "\n"))

    q()

}

#suppressMessages(datatest <- readr::read_tsv(args$Input,col_names = c("ID1","ID2","interaction","chr1","chr2")))    
#sum(datatest$interaction)

suppressMessages(library(dplyr))
uniqw <- read.table(args$Input)
colnames(uniqw) <- c("ID1","ID2","NUM")
uniqw <- mutate(uniqw, ID3=paste0(ID1,ID2))
total <- summarise(group_by(uniqw,ID3),sum(NUM))
uniquniq <- select(uniqw, ID1,ID2,ID3)
uniquniq <- unique(uniquniq)
mergeuniq <- merge(y = total, x = uniquniq,by.x = "ID3",by.y = "ID3")
colnames(mergeuniq) <- c("ID3","ID1","ID2","NUM")
mergeuniq <- select(mergeuniq, ID1,ID2,NUM)
write.table(mergeuniq,args$Output,col.names = F, row.names = F,quote = F)