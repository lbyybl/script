#!/bin/bash

# Boyuan-Li

# Thu Jun 21 17:05:08 CST 2018

# This script is used to find the intersection of pair end sequence;
# because maybe due to the cutadaptor maybe one reads in file1.fq, but not in file2.fq

#set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -a <File1> -b <File2> -p <Output_prefix> "
	echo ""
	echo " -a STRING          [required]  The fastq File1"
	echo ""
	echo " -b STRING          [required]  The fastq File2"
	echo ""
	echo " -p STRING          [required]  The prefix of output file; total three files intersection file; and two unintersection file"
	echo ""
	echo " -h                 help"
	echo -e "\tThis script is used to find the intersection of pair end sequence;\
	because maybe due to the cutadaptor maybe one reads in file1.fq, but not in file2.fq"
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "a:b:p:h" optionName
do
	case $optionName in
		a) File1="$OPTARG";;
		b) File2="$OPTARG";;
		p) Output_prefix="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $File1 = "" ]]; then
	echo " -a the File1 file is needed "
	exit 1
elif [[ ! -f $File1 ]]; then
	echo "$File1:   is not found"
	exit 2
fi


if [[ $File2 = "" ]]; then
	echo " -b the File2 file is needed "
	exit 1
elif [[ ! -f $File2 ]]; then
	echo "$File2:   is not found"
	exit 2
fi


if [[ $Output_prefix = "" ]]; then
	echo " -p the Output_prefix STRING is needed "
	exit 1
fi


function trans_format
{
	input_file=${1}
	output_file=${2}
	awk '{print $1" " $2 "\n" $3 "\n" $4 "\n" $5 }' ${input_file} >> ${output_file}
}

	cat ${File1} | xargs -l4 | sort -t ":" -k 5n -k 6n -k 7n -o ${Output_prefix}_sort1.fastq
	cat ${File2} | xargs -l4 | sort -t ":" -k 5n -k 6n -k 7n -o ${Output_prefix}_sort2.fastq
	join -1 1 -2 1 ${Output_prefix}_sort1.fastq ${Output_prefix}_sort2.fastq >> ${Output_prefix}.interset.fastq
	join -1 1 -2 1 -v1 ${Output_prefix}_sort1.fastq ${Output_prefix}_sort2.fastq >> ${Output_prefix}.uninterset1.fastq
	join -1 1 -2 1 -v2 ${Output_prefix}_sort1.fastq ${Output_prefix}_sort2.fastq >> ${Output_prefix}.uninterset2.fastq
	

	cut -d " " -f 1-5 ${Output_prefix}.interset.fastq >> ${Output_prefix}.R1.oneline.fastq
	cut -d " " -f 1,6-9 ${Output_prefix}.interset.fastq >> ${Output_prefix}.R2.oneline.fastq
	trans_format ${Output_prefix}.R1.oneline.fastq ${Output_prefix}_interset_R1.fastq
	trans_format ${Output_prefix}.R2.oneline.fastq ${Output_prefix}_interset_R2.fastq
	trans_format ${Output_prefix}.uninterset1.fastq ${Output_prefix}.uninterset_R1.fastq
	trans_format ${Output_prefix}.uninterset2.fastq ${Output_prefix}.uninterset_R2.fastq

	rm ${Output_prefix}_sort1.fastq ${Output_prefix}_sort2.fastq ${Output_prefix}.interset.fastq ${Output_prefix}.uninterset1.fastq ${Output_prefix}.uninterset2.fastq ${Output_prefix}.R1.oneline.fastq ${Output_prefix}.R2.oneline.fastq
