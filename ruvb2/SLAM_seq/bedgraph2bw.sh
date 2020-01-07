#!/bin/bash

# Boyuan_Li

# Fri Apr 12 15:56:07 2019

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -n <name> -o <output_name> -i <input_dir> "
	echo ""
	echo " -n	string       	[required] "
	echo ""
	echo " -o	string       	[required] "
	echo ""
	echo " -i	dir          	[required] "
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "n:o:i:h" optionName
do
	case $optionName in
		n) name="$OPTARG";;
		o) output_name="$OPTARG";;
		i) input_dir="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $name = "" ]]; then
	echo " the $name string is needed "
	exit 1
fi

if [[ $output_name = "" ]]; then
	echo " the $output_name string is needed "
	exit 1
fi

if [[ $input_dir = "" ]]; then
	echo "the $input_dir file is needed "
	exit 1
elif [[ ! -d $input_dir ]]; then
	 echo "$input_dir:   is not found"
	exit 2
fi
file1=${name}.fq_slamdunk_mapped_filtered_tcount_mins.bedgraph
file2=${name}.fq_slamdunk_mapped_filtered_tcount_plus.bedgraph

cat <(awk 'OFS="\t"{print $1,$2,$3,$4}' ${input_dir}/${file2}) \
 <(awk 'OFS="\t"{re=-$4;print $1,$2,$3,re}' ${input_dir}/${file1}) \
  | sort -k 1,1 -k 2n,2 -k 3n,3  | uniq >> ${output_name%.*}.bedgraph
  
bedGraphToBigWig ${output_name%.*}.bedgraph /home/boyuanli/tools/hic-pro/HiC-Pro_2.9.0/annotation/chrom_mm10.sizes ${output_name%.*}.bw
  