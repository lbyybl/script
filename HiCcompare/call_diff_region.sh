#!/bin/bash

# Boyuan-Li

# Thu Nov  8 22:36:41 CST 2018

# it's used to merge all step to a pipele, so that you can got the different diff region with the .hic file and genome region and resolution

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -w <Wt_hic_file> -k <Ko_hic_file> -r <Resolution> -g <Genome_bed_file> -o <Output_file> "
	echo ""
	echo " -w STRING          [required] WT .hic file"
	echo ""
	echo " -k STRING          [required] KO .hic file"
	echo ""
	echo " -r STRING          [required] The resolution to call difference"
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
while getopts "w:k:r:g:o:h" optionName
do
	case $optionName in
		w) Wt_hic_file="$OPTARG";;
		k) Ko_hic_file="$OPTARG";;
		r) Resolution="$OPTARG";;
		g) Genome_bed_file="$OPTARG";;
		o) Output_file="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Wt_hic_file = "" ]]; then
	echo " -w the Wt_hic_file file is needed "
	exit 1
elif [[ ! -f $Wt_hic_file ]]; then
	echo "$Wt_hic_file:   is not found"
	exit 2
fi


if [[ $Ko_hic_file = "" ]]; then
	echo " -k the Ko_hic_file file is needed "
	exit 1
elif [[ ! -f $Ko_hic_file ]]; then
	echo "$Ko_hic_file:   is not found"
	exit 2
fi


if [[ $Resolution = "" ]]; then
	echo " -r the Resolution STRING is needed "
	exit 1
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
	ko_file=$(readlink -e $Ko_hic_file)
	wt_file=$(readlink -e $Wt_hic_file)
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

#--- got the prefix
        filename=${Ko_hic_file##*/}
        prefix=$(echo $filename | sed -e "s/_ko_\|_KO_\|_Ko_\|.hic/_/g")	
#--- run
	cd $Output_dir
	#--- got the diff file and all diff file
	/home/boyuanli/bashscript/bin/HiCcompare/juicer2diffrence.sh -k $ko_file -w $wt_file -r $Resolution
	cat signifent${prefix}chr*.diff >> all_signifent_${prefix}.diff
	
	rm signifent${prefix}chr*.diff ${prefix}{Ko,Wt}*.observed ${prefix}chr*.diff
	sed '/start1/ d' all_signifent_${prefix}.diff | sort -k 1,1 -k 2n,2 -k 3n,3 -o all_signifent_${prefix}.diff
	#--- find the interaction region interacting with the diff loop file
	/home/boyuanli/bashscript/bin/HiCcompare/find_region.sh -g $Genome_bed_file -d all_signifent_${prefix}.diff -o ${out_Prefix}_region.bed
	#--- select the first 4 str from the inter region
	Rscript /home/boyuanli/bashscript/bin/HiCcompare/trans_str.R -f ${out_Prefix}_region.bed -o ${out_Prefix}_str.bed
	#--- got the region in different class
	python /home/boyuanli/bashscript/bin/HiCcompare/GNG.py -i ${out_Prefix}_str.bed -o ${Output_file_name}


	


