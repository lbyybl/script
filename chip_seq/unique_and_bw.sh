#!/bin/bash

# Boyuan_Li

# Sat Apr 20 21:27:07 2019

# extract the unique mapping reads and trans it into bigwig and call peak

set -eo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -b <bam_file> -o <output_dir> -p <prefix> -B <blacklist> "
	echo " This script default the species is mouse, if you want apply to human please change the genome size"
	echo ""
	echo " -b	file         	[required] the bam file you want to extract the unique mapping reads"
	echo ""
	echo " -o	dir          	[required] the output dirctory"
	echo ""
	echo " -p	string       	[required] the output prefix"
	echo ""
	echo " -B	file         	[required] the blacklist file"
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "b:o:p:B:h" optionName
do
	case $optionName in
		b) bam_file="$OPTARG";;
		o) output_dir="$OPTARG";;
		p) prefix="$OPTARG";;
		B) blacklist="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $bam_file = "" ]]; then
	echo "the $bam_file file is needed "
	exit 1
elif [[ ! -f $bam_file ]]; then
	echo "$bam_file:   is not found"
	exit 2
fi

if [[ $output_dir = "" ]]; then
	echo "the $output_dir file is needed "
	exit 1
elif [[ ! -d $output_dir ]]; then
	 echo "$output_dir:   is not found"
	exit 2
fi

if [[ $prefix = "" ]]; then
	echo " the $prefix string is needed "
	exit 1
fi

if [[ -z $blacklist ]]; then
	echo "using the default blacklist mm10 "
	blacklist=/DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed
elif [[ ! -f $blacklist ]]; then
	echo "$blacklist:   is not found"
	exit 2
fi

samtools view -hF4 $bam_file | grep -v "XS:" | samtools view -b - >> ${output_dir}/${prefix}_uniqe.bam
samtools index ${output_dir}/${prefix}_uniqe.bam
bamCoverage --bam ${output_dir}/${prefix}_uniqe.bam -o ${output_dir}/${prefix}_uniqe_rpgc.bw --binSize 10 --extendReads 300 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --blackListFileName ${blacklist} --ignoreDuplicates
mkdir -p ${output_dir}/peak
macs2 callpeak -t ${output_dir}/${prefix}_uniqe.bam -n ${prefix} -B -g mm --outdir ${output_dir}/peak --nomodel --shift 100 --extsize 300 -q 0.0001 
