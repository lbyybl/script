#!/bin/R
# 2018/7/23
# This scropt is used to sort the matrix
library(getopt)
command=matrix(c("Input_compartment","c",1,"character",
                 "Output_file_name","o",1,"character",
			"Input_matrix_file","m",1,"character",
			"Resulotion","r",1,"integer",
                 "help","h",0,"logical"),byrow=T,ncol=4)

args=getopt(command)

if (!is.null(args$help) || is.null(args$Input_compartment) || is.null(args$Output) || is.null(args$Input_matrix_file) || is.null(args$Resulotion) ) {

    cat(paste(getopt(command, usage = T), "\n"))

    q()

}
library(tidyverse)


sparse2dense<-function(mat, st, en){
  colnames(mat) <- c("bin1", "bin2", "interaction")
  dense <- data.frame(bin1=rep(st:en, each=(en-st+1)), 
                      bin2=rep(st:en, en-st+1))
  sparse <- mat
  sparse <- sparse[order(sparse$bin1,sparse$bin2),,]
  sparse <- unique(sparse)
  dense <- merge(sparse, dense, by=c("bin1", "bin2"), all=T)
  dense$interaction[is.na(dense$interaction)] <- 0
  lower <- dense[order(dense$bin2), ]
  dense$interaction <- apply(cbind(dense$interaction, lower$interaction), 1, max)
  dense$chr <- unique(dense$chr)[1]
  return(dense)
}


  #--- for new
  new_chr6zscorefile <- read_tsv(args$Input_compartment, col_names = T)  
  nrow(new_chr6zscorefile) # 300
  new_chr6zscorefile <- select(new_chr6zscorefile,start, end, eigen1,index)
  new_chr6zscorefile$eigen1[is.na(new_chr6zscorefile$eigen1)] <- 0
  
#--- try a new sor method
  
  #--- read the chr6 file
    chr6_test <- read_tsv(args$Input_matrix_file, col_names = F)
    colnames(chr6_test) <- c("bin1","bin2","interaction")
    chr6_test_bin <- chr6_test %>%
      mutate(bin1=bin1/args$Resulotion +1 ,bin2=bin2/args$Resulotion +1)
    
    egivector_vec <- new_chr6zscorefile$eigen1

    st_n <- min(chr6_test_bin$bin1, chr6_test_bin$bin2)
    en_n <- max(chr6_test_bin$bin1, chr6_test_bin$bin2)
    chr6_test_bin_dense <- sparse2dense(chr6_test_bin,st_n,en_n)
    chr6_test_bin_egi <- chr6_test_bin_dense 
    chr6_test_bin_egi$bin1_egi <- egivector_vec[chr6_test_bin_egi$bin1]
    chr6_test_bin_egi$bin2_egi <- egivector_vec[chr6_test_bin_egi$bin2]
    chr6_test_bin_egi_sort <- chr6_test_bin_egi[order(chr6_test_bin_egi$bin1_egi,chr6_test_bin_egi$bin2_egi),,]
    new_matrix_sort <- data.frame(bin1=c(rep(st_n:en_n,each=(en_n-st_n+1))),
                                  bin2=c(rep(st_n:en_n,en_n-st_n+1))) 
    new_matrix_sort$interaction <- chr6_test_bin_egi_sort$interaction
    write.csv(new_matrix_sort, args$Output, row.names = F, col.names = T)
    


