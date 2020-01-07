#!/bin/bash

# Boyuan-Li

# 2019/1/14



# it's used the sampled region to see if the pie plot is random

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -w <Wt_hic_file> -k <Ko_hic_file> -r <Resolution> -g <Genome_bed_file> -o <Output_file> "
	echo ""
	echo " -s STRING          [required] sample region file file"
	echo ""
	echo " -g STRING          [required] Gene_bed_file can be file merged by promoter enhancer gff file"
	echo ""
	echo " -o STRING          [required] Output_file name and the prefix also will be used and the dir will be used as output dir"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "s:g:o:h" optionName
do
	case $optionName in
		s) Sample_region_file="$OPTARG";;
		g) Genome_bed_file="$OPTARG";;
		o) Output_file="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Sample_region_file = "" ]]; then
	echo " -w the $Sample_region_file file is needed "
	exit 1
elif [[ ! -f $Sample_region_file ]]; then
	echo "$Sample_region_file:   is not found"
	exit 2
fi

if [[ $Genome_bed_file = "" ]]; then
	echo " -g the Genome_bed_file file is needed "
	exit 1
elif [[ ! -f $Genome_bed_file ]]; then
	echo "$Genome_bed_file:   is not found"
	exit 2
fi


if [[ $Output_file = "" ]]; then
        echo " -o the Output_file STRING is needed "
        exit 1
fi



#--- got the abs path
	diff_file=$(readlink -e ${Sample_region_file})
	Gene_file=$(readlink -e $Genome_bed_file)
	Output_file_name=${Output_file##*/}
	Output_dir="./"
	out_Prefix=${Output_file_name%.*}

	dir=$(echo -e "${Output_file}" | grep "/")
	if [[ $dir != "" && $dir != "./" ]]; then
			Output_dir=${Output_file%/*}
			if [[ ! -d $Output_dir ]]; then
				mkdir -p ${Output_dir}
			fi
	fi
	Output_dir=$(readlink -e ${Output_dir})

#--- run
	cd $Output_dir
	
	
	
	sed '/start1/ d' ${diff_file} | sort -k 1,1 -k 2n,2 -k 3n,3 -o ${diff_file}
	
	awk 'OFS="\t" {print $1,$2,$2,$1,$3,$3,"1","1","3","4","5","6","7","8","9","10","11","12"}' ${diff_file} >> ${diff_file##*/}_diff
	
	#--- find the interaction region interacting with the diff loop file
	
	/home/boyuanli/bashscript/bin/HiCcompare/find_region.sh -g $Genome_bed_file -d ${diff_file##*/}_diff -o ${out_Prefix}_region.bed
	
	#--- select the first 4 str from the inter region
	Rscript /home/boyuanli/bashscript/bin/HiCcompare/trans_str.R -f ${out_Prefix}_region.bed -o ${out_Prefix}_str.bed
	#--- got the region in different class
	python /home/boyuanli/bashscript/bin/HiCcompare/GNG.py -i ${out_Prefix}_str.bed -o ${Output_file_name}


	


