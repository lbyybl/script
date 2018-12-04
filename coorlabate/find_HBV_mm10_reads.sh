#!/bin/bash

# Boyuan-Li

# Sun Nov 11 21:48:18 CST 2018

# This file was used to find HBV mm10 mixture reads form hicpro results
# the pair.bam file was need in bowtie_results/bwt2/ dir or hic_results/data/(remove duplicated) 

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -b <BAM_file> -o <Output_file> "
	echo ""
	echo " -b STRING          [required] BAM file produced by HiCpro"
	echo ""
	echo " -o STRING          [required] the output file name samfile"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "b:o:h" optionName
do
	case $optionName in
		b) BAM_file="$OPTARG";;
		o) Output_file="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $BAM_file = "" ]]; then
	echo " -b the BAM_file file is needed "
	exit 1
elif [[ ! -f $BAM_file ]]; then
	echo "$BAM_file:   is not found"
	exit 2
fi


if [[ $Output_file = "" ]]; then
	echo " -o the Output_file STRING is needed "
	exit 1
fi

	samtools view $BAM_file -H >> $Output_file
	samtools view $BAM_file | awk '{if ($3!="HBV" && $7 == "HBV") print $0}' >> $Output_file
