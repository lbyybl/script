#!/bin/bash

# Boyuan-Li

# Wed Jul 11 18:26:41 CST 2018

# this script is used to call compartment form the result producted by hic-pro;
# it need the script of matrixhic-resorter

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -b <hicpro_result_bed_file> -m <hicpro_result_ice_matrix_file> -o <Output_dire> -g <Genome_version> -s <Genome_size_file_ued_by_hicpro> "
	echo ""
	echo " -b STRING          [required] the bed file producted by hic-pro, in the hic_results/matrix/ and raw directory"
	echo ""
	echo " -m STRING          [required] the matrix file producted by hic-pro, in the hic_results/matrix/ and ice directory"
	echo ""
	echo " -o STRING          [required] the output directory"
	echo ""
	echo " -g STRING          [required] the Genome version you used"
	echo ""
	echo " -s STRING          [required] the Genome size file used by hic-pro"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "b:m:o:g:s:h" optionName
do
	case $optionName in
		b) hicpro_result_bed_file="$OPTARG";;
		m) hicpro_result_ice_matrix_file="$OPTARG";;
		o) Output_dire="$OPTARG";;
		g) Genome_version="$OPTARG";;
		s) Genome_size_file_ued_by_hicpro="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $hicpro_result_bed_file = "" ]]; then
	echo " -b the hicpro_result_bed_file file is needed "
	exit 1
elif [[ ! -f $hicpro_result_bed_file ]]; then
	echo "$hicpro_result_bed_file:   is not found"
	exit 2
fi


if [[ $hicpro_result_ice_matrix_file = "" ]]; then
	echo " -m the hicpro_result_ice_matrix_file file is needed "
	exit 1
elif [[ ! -f $hicpro_result_ice_matrix_file ]]; then
	echo "$hicpro_result_ice_matrix_file:   is not found"
	exit 2
fi


if [[ $Output_dire = "" ]]; then
	echo " -o the Output_dir directory is needed "
	exit 1
elif [[ ! -d $Output_dire ]]; then
	echo "$Output_dire:   is not found"
	echo "mkdir $Output_dire ......"
	mkdir -p $Output_dire
fi


if [[ $Genome_version = "" ]]; then
	echo " -g the Genome_version STRING is needed "
	exit 1
fi


if [[ $Genome_size_file_ued_by_hicpro = "" ]]; then
	echo " -s the Genome_size_file_ued_by_hicpro file is needed "
	exit 1
elif [[ ! -f $Genome_size_file_ued_by_hicpro ]]; then
	echo "$Genome_size_file_ued_by_hicpro:   is not found"
	exit 2
fi


	hicpro_dir=$(dirname $(which HiC-Pro))
	cworldscript_dir=$(dirname $(which addMatrixHeaders.pl))
	Input_bedfile=$(readlink -e $hicpro_result_bed_file)
	Input_matrixfile=$(readlink -e $hicpro_result_ice_matrix_file)
	Output_dir=$(readlink -e $Output_dire)
	Genomenane=$Genome_version   # $Genome_version is the genome version, you should give;
	Genomesizefile=$(readlink -e $Genome_size_file_ued_by_hicpro)
	filename=${Input_matrixfile##*/}
	Input_dir=${Input_matrixfile%/*}
	prefix=${filename%%.*}
	cd $Output_dir
#####################################################################################################################################
# trans the spare matrix into dense matrix; $hicpro_result_bed_file is the bed file in the hicrop-result/raw; and $hicpro_result_ice_matrix_file is the matrix file in the hicpro-result/ICE
#######################################################################################################################################
	$hicpro_dir/utils/sparseToDense.py -b $Input_bedfile $Input_matrixfile --perchr

	ls $Input_dir/${prefix}_chr* | parallel mv {} ./   # mv the dense matrix produced by sparseToDense.py to Output_dir
	# for example $hicpro_dir/utils/sparseToDense.py -b sample1124_wt_100000_abs.bed sample1124_wt_100000_iced.matrix

####################################################################################################################################
# add header for the dense matrix
#####################################################################################################################################
# make header for the dense matrix

	awk -v Genomenane="$Genomenane" '{printf "%i|%s|%s:%i-%i\n",$Genome_version,Genomenane,$hicpro_result_bed_file,$hicpro_result_ice_matrix_file,$Output_dire}'  $Input_bedfile >> ${prefix}_header # make header file with the input bed file

	chrname=$(cut -f 1 <(egrep -v "chrM" $Genomesizefile) | xargs | sed 's/ /\\n/g')     # obtain the chr name expect the chrM and chrMT
	echo -e "$chrname" | parallel -j1 egrep "{}:" ${prefix}_header '>>' ${prefix}_{}_header  # got the header file per chr


#egrep "chrM" ${prefix}_header '>>' ${prefix}_chrY_header
	rm ${prefix}_header # i think no use
	# example: awk '{printf "%i|mm10|%s:%i-%i\n",$Genome_version,$hicpro_result_bed_file,$hicpro_result_ice_matrix_file,$Output_dire}' ./hicpro_latest_test/hic_results/matrix/dixon_2M/raw/1000000/dixon_2M_1000000_abs.bed >> header

	# add header for the dense matrix

	echo -e "$chrname" | parallel -j10 perl -I $cworldscript_dir/../../lib/ $cworldscript_dir/addMatrixHeaders.pl -i ${filename%.*}_{}_dense.matrix --xhf ${prefix}_{}_header --yhf ${prefix}_{}_header


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
