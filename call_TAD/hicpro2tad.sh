#!/bin/bash

# Boyuan-Li

# Tue Nov 27 20:36:24 CST 2018

# It's used to call tad with the .matrix and .bed file produced by hicpro.

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -b <Bed_file_by_hicpro_corresponding_to_resoluton> -m <matrix_file_by_hicpro_corresponding_to_resoluton> -o <Output> -g <Genome_version> "
	echo ""
	echo " -b STRING          [required] Bed_file_by_hicpro_corresponding_to_resoluton"
	echo ""
	echo " -m STRING          [required] matrix_file_by_hicpro_corresponding_to_resoluton"
	echo ""
	echo " -o STRING          [required] The output directory"
	echo ""
	echo " -g STRING          [required] The genome version you used, such as mm10,mm9,hg19 et al."
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "b:m:o:g:h" optionName
do
	case $optionName in
		b) Bed_file_by_hicpro_corresponding_to_resoluton="$OPTARG";;
		m) matrix_file_by_hicpro_corresponding_to_resoluton="$OPTARG";;
		o) Output="$OPTARG";;
		g) Genome_version="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Bed_file_by_hicpro_corresponding_to_resoluton = "" ]]; then
	echo " -b the Bed_file_by_hicpro_corresponding_to_resoluton file is needed "
	exit 1
elif [[ ! -f $Bed_file_by_hicpro_corresponding_to_resoluton ]]; then
	echo "$Bed_file_by_hicpro_corresponding_to_resoluton:   is not found"
	exit 2
fi


if [[ $matrix_file_by_hicpro_corresponding_to_resoluton = "" ]]; then
	echo " -m the matrix_file_by_hicpro_corresponding_to_resoluton file is needed "
	exit 1
elif [[ ! -f $matrix_file_by_hicpro_corresponding_to_resoluton ]]; then
	echo "$matrix_file_by_hicpro_corresponding_to_resoluton:   is not found"
	exit 2
fi

if [ ! -d $Output ]; then
	echo -e "$Output does not exists"
	echo -e "creating $Output "
	mkdir -p $Output
else 
	echo -e "$Output exists"
fi

if [[ $Genome_version = "" ]]; then
	echo " -g the Genome_version STRING is needed "
	exit 1
fi

	
	hicpro_dir=$(dirname $(which HiC-Pro))
	java_dir=/DATA/work/lbyybl/tools/juicer_fml/
	cworldscript_dir=$(dirname $(which addMatrixHeaders.pl))
	Input_bedfile=$(readlink -e $Bed_file_by_hicpro_corresponding_to_resoluton)
	Input_matrixfile=$(readlink -e $matrix_file_by_hicpro_corresponding_to_resoluton)
	Output_dir=$(readlink -e $Output)
	Genomename=$Genome_version   # $Genome_version is the genome version, you should give;
	filename=${Input_matrixfile##*/}
	Input_dir=${Input_matrixfile%/*}
	prefix=${filename%%.*}
	cd $Output_dir
	

#----- trans the sparse matrix into dense
	#--- get the bed file for every chr
	echo {{1..19},X,Y} | xargs -n1 | parallel  egrep -w \"chr{}\" ${Input_bedfile} ">>" chr{}_${prefix}.bed
	#--- trans the sparse to dense 
	#set +euo pipefail
	#ln -s ${Input_bedfile} ${Input_bedfile##*/}
	#ln -s ${Input_matrixfile} ${Input_matrixfile##*/}
	$hicpro_dir/utils/sparseToDense.py -b ${Input_bedfile} ${Input_matrixfile} -c 
	#set -euo pipefail
	ls $Input_dir/${prefix}_chr* | parallel mv {} ./
	
#--- add header for the dense matrix 
	#--- make header file
	echo {{1..19},X,Y} | xargs -n1 | parallel awk -v Genomename="$Genomename" \'{printf \"%i\|%s\|%s\:%i-%i\\n\",\$4,Genomename,\$1,\$2,\$3}\' chr{}_${prefix}.bed ">>" chr{}_${prefix}_header
	#--- test whether header's line is equal dense file's line 
	# for i in {{1..19},X,Y}; \
		# do \
		# Row_line=$(cat chr${i}_${prefix}_dense.matrix | wc -l)
		# Header_line=$(cat chr${i}_${prefix}_header | wc -l)
		# if [[ ${Row_line} != ${Header_line} ]]; then
			# sed -i '$d' chr${i}_${prefix}_header
		# fi
	# done	
	#--- add header
	echo {{1..19},X,Y} | xargs -n1 | parallel -j10 perl -I $cworldscript_dir/../../lib/ $cworldscript_dir/addMatrixHeaders.pl -i ${prefix}_chr{}_dense.matrix --xhf chr{}_${prefix}_header --yhf chr{}_${prefix}_header

	
#--- call compartment with cworld
	#--- change 0.0 into nan
	for i in {{1..19},X,Y}; \
		do \
		sed -e '/^##/d' -e 's/\t0.0\t/\tnan\t/g' -e 's/\t0.0\t/\tnan\t/g' -e 's/\t0.0$/\tnan/g' <(less ${prefix}_chr${i}_dense.addedHeaders.matrix.gz) \
		>> chr${i}_${prefix}_dense.addedHeaders.matrix ; \
	done
	# calculate the eg1 for header matrix
	echo {{1..19},X,Y} | xargs -n1 | parallel -j10 perl -I $cworldscript_dir/../../lib/ $cworldscript_dir/matrix2insulation.pl -i \
	chr{}_${prefix}_dense.addedHeaders.matrix --is 100000 --im sum 

#--- rm useless files
	rm -f *.{gz,boundaries.bed,log,pdf,boundaries,bedGraph,png,addedHeaders.matrix}
	rm -f *_dense.matrix
	rm -f chr*.bed
	echo {{1..19},X,Y} | xargs -n1 | parallel -j10 rm ${prefix}_chr{}_dense.matrix  chr{}_${prefix}_header ## remove matrix and header in line 121

	echo {{1..19},X,Y} | xargs -n1 | parallel -j10 rm ${prefix}_chr{}_dense.addedHeaders.matrix # rm the matrix in line 132
#--- rm file unless
	echo {{1..19},X,Y} | xargs -n1 | parallel rm chr{}_${prefix}.bed # chr{}_${prefix}_${Juicer_type}.matrix chr{}_${prefix}.${Juicer_type} # line 105 108 113


