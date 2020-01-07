#!/bin/bash

# Boyuan_Li

# Wed May  1 15:08:25 2019

# This is used to sample the reads by the number you give.

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -s <Sample> -n <Number> -p <Prefix> "
	echo ""
	echo " -s	string       	[required] The sample you will sampled"
	echo ""
	echo " -n	string       	[required] The number you want sampled"
	echo ""
	echo " -p	string       	[required] The Prefix you will give the output file"
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "s:n:p:h" optionName
do
	case $optionName in
		s) Sample="$OPTARG";;
		n) Number="$OPTARG";;
		p) Prefix="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $Sample = "" ]]; then
	echo " the $Sample string is needed "
	exit 1
fi

if [[ $Number = "" ]]; then
	echo " the $Number string is needed "
	exit 1
fi

if [[ $Prefix = "" ]]; then
	echo " the $Prefix string is needed "
	exit 1
fi

samtools view -H $Sample >> ${Prefix}______.sam
subsample -s 123456 -n $Number <(samtools view -F4 $Sample) >> ${Prefix}______.sam
samtools sort -@ 10 ${Prefix}______.sam >> ${Prefix}.bam
samtools index ${Prefix}.bam
rm -f ${Prefix}______.sam
samtools view -F4 ${Prefix}.bam | wc -l