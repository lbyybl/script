#!/bin/bash

# Boyuan_Li

# Thu Nov  7 21:41:09 2019

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -b <name_sort_bam> -p <output_prefilx> "
	echo ""
	echo " -b	file         	[required] The input bam file should be sort by name"
	echo ""
	echo " -p	string       	[required] output prefix"
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "b:p:h" optionName
do
	case $optionName in
		b) name_sort_bam="$OPTARG";;
		p) output_prefilx="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $name_sort_bam = "" ]]; then
	echo "the $name_sort_bam file is needed "
	exit 1
elif [[ ! -f $name_sort_bam ]]; then
	echo "$name_sort_bam:   is not found"
	exit 2
fi

if [[ $output_prefilx = "" ]]; then
	echo " the $output_prefilx string is needed "
	exit 1
fi

bedtools bamtobed -i $name_sort_bam | pigz > ${output_prefilx}.bed.gz
/WORK/lbyybl/tools/NucTools/nucleosome_repeat_length.pl -m 10000000 -in ${output_prefilx}.bed.gz -out ${output_prefilx}.txt.gz
/usr/bin/Rscript /WORK/lbyybl/tools/NucTools/misc/plotNRL.R --input=${output_prefilx}.txt.gz  --out=${output_prefilx}.pdf --dir=./
