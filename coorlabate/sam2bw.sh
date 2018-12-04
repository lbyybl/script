#!/bin/bash

# Boyuan-Li

# Tue Oct 23 15:06:35 CST 2018

# this is trams sam file to sort bam and bw file

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -s <Sam_file> -o <Output_bwfile> "
	echo ""
	echo " -s STRING          [required] sam file"
	echo ""
	echo " -o STRING          [required] output file"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "s:o:h" optionName
do
	case $optionName in
		s) Sam_file="$OPTARG";;
		o) Output_bwfile="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Sam_file = "" ]]; then
	echo " -s the Sam_file file is needed "
	exit 1
elif [[ ! -f $Sam_file ]]; then
	echo "$Sam_file:   is not found"
	exit 2
fi


if [[ $Output_bwfile = "" ]]; then
	echo " -o the Output_bwfile STRING is needed "
	exit 1
fi

#--- trans sam to sort bam

#	samtools view -bS $Sam_file > ${Sam_file}.bam
	samtools sort ${Sam_file} > ${Sam_file%.*}_sorted.bam
	
#--- index bam file

	samtools index ${Sam_file%.*}_sorted.bam
	
#--- creating bigwig 

	bamCoverage -b ${Sam_file%.*}_sorted.bam -o $Output_bwfile --normalizeUsing RPKM
	
rm $Sam_file ${Sam_file%.*}.bam 
