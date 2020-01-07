#!/bin/bash
set -exu
# this script is used to call compartment form the result producted by hic-pro;
# it need the script of matrixhic-resorter
# Boyuan-Li
# 2018-5-28


function helps
{
        echo ""
        echo -e "Usage: $0 hicpro_result/raw/bed hicpro/ice/matrix output_dir genome_version(such as mm9,mm10,hg18 et al.) genome_size_file_used_by_hic_pro"
   

}

while getopts "h" optionName
do
        case "$optionName" in
                
                h)
                        helps
                        exit 0
                        ;;
        esac
done


#source ~/.bashrc
hicpro_dir=$(dirname $(which HiC-Pro))
cworldscript_dir=$(dirname $(which addMatrixHeaders.pl))

Input_bedfile=$(readlink -e $1)
Input_matrixfile=$(readlink -e $2)
Output_dir=$(pwd)/$3
Genomenane=$4   # $4 is the genome version, you should give;
Genomesizefile=$(readlink -e $5)
filename=${Input_matrixfile##*/}
Input_dir=${Input_matrixfile%/*}
prefix=${filename%%.*}
mkdir -p $Output_dir
cd $Output_dir
#####################################################################################################################################
# trans the spare matrix into dense matrix; $1 is the bed file in the hicrop-result/raw; and $2 is the matrix file in the hicpro-result/ICE
#######################################################################################################################################
$hicpro_dir/utils/sparseToDense.py -b $Input_bedfile $Input_matrixfile --perchr

ls $Input_dir/${prefix}_chr* | parallel mv {} ./   # mv the dense matrix produced by sparseToDense.py to Output_dir
# for example $hicpro_dir/utils/sparseToDense.py -b ./trimlinker4_resullt/hic_results/matrix/sample1124_wt/raw/100000/sample1124_wt_100000_abs.bed ./trimlinker4_resullt/hic_results/matrix/sample1124_wt/iced/100000/sample1124_wt_100000_iced.matrix

####################################################################################################################################
# add header for the dense matrix
#####################################################################################################################################
# make header for the dense matrix

awk -v Genomenane="$Genomenane" '{printf "%i|%s|%s:%i-%i\n",$4,Genomenane,$1,$2,$3}'  $Input_bedfile >> ${prefix}_header # make header file with the input bed file

chrname=$(cut -f 1 <(egrep -v "chrM" $Genomesizefile) | xargs | sed 's/ /\\n/g')     # obtain the chr name expect the chrM and chrMT
echo -e "$chrname" | parallel -j1 egrep "{}:" ${prefix}_header '>>' ${prefix}_{}_header  # got the header file per chr
#egrep "chrM" ${prefix}_header '>>' ${prefix}_chrY_header
rm ${prefix}_header # i think no use
# example: awk '{printf "%i|mm10|%s:%i-%i\n",$4,$1,$2,$3}' ./hicpro_latest_test/hic_results/matrix/dixon_2M/raw/1000000/dixon_2M_1000000_abs.bed >> header

# add header for the dense matrix

echo -e "$chrname" | parallel -j10 perl -I $cworldscript_dir/../../lib/ $cworldscript_dir/addMatrixHeaders.pl -i ${filename%.*}_{}_dense.matrix --xhf ${prefix}_{}_header --yhf ${prefix}_{}_header


#example: perl -I $cworldscript_dir/../../lib/ $cworldscript_dir/addMatrixHeaders.pl -i /DATA/work/lbyybl/wangcl/hic/KPNA2/trim_linker/sample1124_wt_100000_iced_dense.matrix --xhf /DATA/work/lbyybl/wangcl/hic/KPNA2/trim_linker/header --yhf /DATA/work/lbyybl/wangcl/hic/KPNA2/trim_linker/header


#########################################################################################################
# call compartment with cworld
##########################################################################################################

#change 0.0 into nan

for i in $(echo -e $chrname); \
do \
sed -e '/^##/d' -e 's/\t0.0\t/\tnan\t/g' -e 's/\t0.0\t/\tnan\t/g' -e 's/\t0.0$/\tnan/g' <(less ${filename%.*}_${i}_dense.addedHeaders.matrix.gz) \
>> ${filename%.*}_${i}_dense.addedHeaders.matrix ; \
done
# calculate the eg1 for header matrix
echo -e "$chrname" | parallel -j10 perl -I $cworldscript_dir/../../lib/ $cworldscript_dir/matrix2compartment.pl -i ${filename%.*}_{}_dense.addedHeaders.matrix
# example: perl -I $cworldscript_dir/../../lib/ $cworldscript_dir/matrix2compartment.pl -i sample1124_wt_100000_iced_dense.addedHeaders.matrix
###############################################################################################################################################
# regot the resort matrix with the script: matrixhic-resorter
###############################################################################################################################################
# get the correspording bed file 
#
# resort 
#echo -e "$chrname" | parallel -j10 /home/boyuanli/bashscript/bin/matrixhic-resorter ${filename%.*}_{}_dense.addedHeaders.zScore.compartments \
#${prefix}_{}.bed 

# to got the bed(no chrM/chrMT); spare matrix (no chrM/chrMT); and mergeed eg1 file

# bed(no chrM/chrMT);
echo -e "$chrname" | parallel -j1 egrep "{}" $Input_bedfile '>>' ${prefix}_no_chrM.bed # to get an no chrM bed file

# spare matrix (no chrM/chrMT);
chrM_correspording_binnumber=$(grep "chrM" $Input_bedfile | cut -f 4 | xargs | sed -e 's/ /\\s\|\\s/g' -e 's/^/\\s/g' -e 's/$/\\s/g')
chrM_correspording_binnumber2=$(grep "chrM" $Input_bedfile | cut -f 4 | xargs | sed -e 's/ /\\s\|^/g' -e 's/^/\^/g' -e 's/$/\\s/g')
egrep -v "$chrM_correspording_binnumber|$chrM_correspording_binnumber2" $Input_matrixfile >> ${prefix}_no_chrM.matrix # to get a sparse matrix without chrM and chrMT

# mergeed eg1 file
echo -e "$chrname" | parallel -j1 cat ${filename%.*}_{}_dense.addedHeaders.zScore.compartments >> ${filename%.*}_merge_zSxore.compartments # merge the zScore.compartments file of all chr (except chrM and chrMT)

##########################################################################################################################################################
###got resort file
## after this step, you will got two files, one is the ${prefix}_no_chrM_resort.bed,another is the sparse matrix that have changed bin number :${prefix}_no_chrM_sparse.matrix
## and you can use it to draw heatmap (use python or Hicplotter)
##########################################################################################################################################################
matrixhic-resorter ${filename%.*}_merge_zSxore.compartments \
${prefix}_no_chrM.bed ${prefix}_no_chrM.matrix ${prefix}_no_chrM_sparse.matrix # change the bin number according the eg1

### remove some file that don't need at all.
echo -e "$chrname" | parallel -j10 rm ${filename%.*}_{}_dense.matrix  ${prefix}_{}_header ## remove matrix and header in line 62

echo -e "$chrname" | parallel -j10 rm ${filename%.*}_{}_dense.addedHeaders.matrix # rm the matrix in line 77

rm ${prefix}_no_chrM.bed ${prefix}_no_chrM.matrix ${filename%.*}_merge_zSxore.compartments # rm the file in line 94 99 102 

###finally you can you use ${prefix}_no_chrM_sparse.matrix to draw heatmap




