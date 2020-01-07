#!/bin/bash

# Boyuan-Li

# Thu Jul 19 15:46:03 CST 2018

# It's used to call compartment with the .hic file produced by juicer.

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -b <Bed_file_by_hicpro_corresponding_to_resoluton> -j <Hic_file_by_juicer> -r <Resolution> -o <Output> -g <Genome_version> -t <calculate type>"
	echo ""
	echo " -b STRING          [required] Bed_file_by_hicpro_corresponding_to_resoluton"
	echo ""
	echo " -j STRING          [required] .hic file produced by juicer pre"
	echo ""
	echo " -r STRING          [required] bin size"
	echo ""
	echo " -o STRING          [required] The output directory"
	echo ""
	echo " -g STRING          [required] The genome version you used, such as mm10,mm9,hg19 et al."
	echo ""
	echo " -t STRING          [optional] the juicer observed or observed/expected"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
Juicer_type="observed"
while getopts "b:j:r:o:g:t:h" optionName
do
	case $optionName in
		b) Bed_file_by_hicpro_corresponding_to_resoluton="$OPTARG";;
		j) Hic_file_by_juicer="$OPTARG";;
		r) Resolution="$OPTARG";;
		o) Output="$OPTARG";;
		g) Genome_version="$OPTARG";;
		t) Juicer_type="$OPTARG";;
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


if [[ $Hic_file_by_juicer = "" ]]; then
	echo " -j the Hic_file_by_juicer file is needed "
	exit 1
elif [[ ! -f $Hic_file_by_juicer ]]; then
	echo "$Hic_file_by_juicer:   is not found"
	exit 2
fi


if [[ $Resolution = "" ]]; then
	echo " -r the Resolution STRING is needed "
	exit 1
fi


if [[ $Output = "" ]]; then
	echo " -o the Output directory is needed "
	exit 1
elif [[ ! -d $Output ]]; then
	echo "$Output:   is not found"
	exit 2
fi


if [[ $Genome_version = "" ]]; then
	echo " -g the Genome_version STRING is needed "
	exit 1
fi

	hicpro_dir=$(dirname $(which HiC-Pro))
	java_dir=/DATA/work/lbyybl/tools/juicer_fml/
	cworldscript_dir=$(dirname $(which addMatrixHeaders.pl))
	Input_bedfile=$(readlink -e $Bed_file_by_hicpro_corresponding_to_resoluton)
	Input_hicfile=$(readlink -e $Hic_file_by_juicer)
	Output_dir=$(readlink -e $Output)
	Genomename=$Genome_version   # $Genome_version is the genome version, you should give;
	filename=${Input_hicfile##*/}
	Input_dir=${Input_hicfile%/*}
	prefix=${filename%%.*}
	cd $Output_dir
	
#---- get the correspording resolution file from the .hic file, note -d can produced the dense file
	#--- sparse
	# echo {1..19} | xargs -n1 | parallel java -jar ${java_dir}/juicer_tools.1.8.9_jcuda.0.8.jar dump ${Juicer_type} \
	# KR ${Input_hicfile} {} {} BP $Resolution chr{}_${prefix}.${Juicer_type}
	#--- dense
	echo {1..19} | xargs -n1 | parallel java -jar ${java_dir}/juicer_tools.1.8.9_jcuda.0.8.jar dump ${Juicer_type} \
	KR ${Input_hicfile} {} {} BP $Resolution chr{}_${prefix}_${Juicer_type}_dense.matrix -d 
	#-- change the bin number
	# echo {1..19} | xargs -n1 | parallel awk -v res="${Resolution}" \'{printf \"%i \\t %i \\t %f \\n\",\$1/res+1,\$2/res+1,\$3}\' chr{}_${prefix}.${Juicer_type} '>>' chr{}_${prefix}_${Juicer_type}.matrix
	# sort -k 1n -k 2n chr{}_${prefix}_${Juicer_type}.matrix -o chr{}_${prefix}_${Juicer_type}.matrix
	
#----- trans the sparse matrix into dense
	#--- get the bed file for every chr
	echo {1..19} | xargs -n1 | parallel  egrep -w \"chr{}\" ${Input_bedfile} ">>" chr{}_${prefix}.bed
	#--- trans the sparse to dense 
	#set +euo pipefail
	# echo {1..19} | xargs -n1 | parallel $hicpro_dir/utils/sparseToDense.py -b chr{}_${prefix}.bed chr{}_${prefix}_${Juicer_type}.matrix 
	#set -euo pipefail
	#ls $Input_dir/${prefix}_chr* | parallel mv {} ./
	
#--- add header for the dense matrix 
	#--- make header file
	echo {1..19} | xargs -n1 | parallel awk -v Genomename="$Genomename" \'{printf \"%i\|%s\|%s\:%i-%i\\n\",\$4,Genomename,\$1,\$2,\$3}\' chr{}_${prefix}.bed ">>" chr{}_${prefix}_header
	#--- test whether header's line is equal dense file's line 
	for i in {1..19}; \
		do \
		Row_line=$(cat chr${i}_${prefix}_${Juicer_type}_dense.matrix | wc -l)
		Header_line=$(cat chr${i}_${prefix}_header | wc -l)
		if [[ ${Row_line} != ${Header_line} ]]; then
			st_line=$(echo ${Row_line}+1 | bc -l)
			sed -i "${st_line},$ d" chr${i}_${prefix}_header
		fi
	done	
	#--- add header
	echo {1..19} | xargs -n1 | parallel -j10 perl -I $cworldscript_dir/../../lib/ $cworldscript_dir/addMatrixHeaders.pl -i chr{}_${prefix}_${Juicer_type}_dense.matrix --xhf chr{}_${prefix}_header --yhf chr{}_${prefix}_header

	
#--- call compartment with cworld
	#--- change 0.0 into nan
	for i in {1..19}; \
		do \
		sed -e '/^##/d' -e 's/\t0.0\t/\tnan\t/g' -e 's/\t0.0\t/\tnan\t/g' -e 's/\t0.0$/\tnan/g' <(less chr${i}_${prefix}_${Juicer_type}_dense.addedHeaders.matrix.gz) \
		>> chr${i}_${prefix}_${Juicer_type}_dense.addedHeaders.matrix ; \
	done
	# calculate the eg1 for header matrix
	echo {1..19} | xargs -n1 | parallel -j10 perl -I $cworldscript_dir/../../lib/ $cworldscript_dir/matrix2compartment.pl -i \
	chr{}_${prefix}_${Juicer_type}_dense.addedHeaders.matrix

#--- rm useless files
	echo {1..19} | xargs -n1 | parallel -j10 rm chr{}_${prefix}_${Juicer_type}_dense.matrix  chr{}_${prefix}_header ## remove matrix and header in line 121

	echo {1..19} | xargs -n1 | parallel -j10 rm chr{}_${prefix}_${Juicer_type}_dense.addedHeaders.matrix # rm the matrix in line 132
#--- rm file unless
	echo {1..19} | xargs -n1 | parallel rm chr{}_${prefix}.bed # chr{}_${prefix}_${Juicer_type}.matrix chr{}_${prefix}.${Juicer_type} # line 105 108 113

