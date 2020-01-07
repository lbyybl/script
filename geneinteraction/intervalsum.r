#/bin/R
# 2018-9-7
# this is used to sum the interaction, only for these data
#library(getopt)
command=matrix(c("Input","i",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)

args=getopt::getopt(command)

if (!is.null(args$help) || is.null(args$Input) ) {

    cat(paste(getopt::getopt(command, usage = T), "\n"))

    q()

}

suppressMessages(datatest <- readr::read_tsv(args$Input,col_names = c("Chr","st","en","strand","ID")))    

sum(as.numeric(datatest$en-datatest$st))
#sum(as.numeric(abs(datatest$en-datatest$st)))
