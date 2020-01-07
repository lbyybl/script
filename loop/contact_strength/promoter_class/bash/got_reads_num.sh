#!/bin/bash

# Boyuan_Li

# Fri Jan 25 16:07:05 2019

# This is used to got the reads in a gene region according to the mapping reads region

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -g <gene_bed> -m <gro_seq_mapping_bed> -p <outprefix> -o <outputdir> "
	echo ""
	echo -e " -g	file         	[required] the bed file contian the region you want to counts reads, \n\t\t\tit's format like chr1    3214481 3670498 NM_001011874    0       -"
	echo ""
	echo -e " -m	file         	[required] the bed file contian the reads that mapping to the genome, \n\t\t\tlike chr1    3000040 3000072 33      0       -, \n\t\t\tyou can use R grnomealigament to got it"
	echo ""
	echo " -p	string       	[required] the prefix the output file will be used"
	echo ""
	echo " -o	dir          	[required] the output dir, difult is ./"
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}

outputdir=./
if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "g:m:p:o:h" optionName
do
	case $optionName in
		g) gene_bed="$OPTARG";;
		m) gro_seq_mapping_bed="$OPTARG";;
		p) outprefix="$OPTARG";;
		o) outputdir="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $gene_bed = "" ]]; then
	echo "the $gene_bed file is needed "
	exit 1
elif [[ ! -f $gene_bed ]]; then
	echo "$gene_bed:   is not found"
	exit 2
fi

if [[ $gro_seq_mapping_bed = "" ]]; then
	echo "the $gro_seq_mapping_bed file is needed "
	exit 1
elif [[ ! -f $gro_seq_mapping_bed ]]; then
	echo "$gro_seq_mapping_bed:   is not found"
	exit 2
fi

if [[ $outprefix = "" ]]; then
	echo " the $outprefix string is needed "
	exit 1
fi

if [[ $outputdir = "" ]]; then
	echo "the $outputdir file is needed "
	exit 1
elif [[ ! -d $outputdir ]]; then
	 echo "$outputdir:   is not found"
	mkdir -p $outputdir
fi

# 统计每个区域有多少reads，落在上面；

# gene_bed
# gro_seq_mapping_bed
# outprefix
# outputdir defuld ./

gene_bed=$(readlink -e ${gene_bed})
gro_seq_mapping_bed=$(readlink -e ${gro_seq_mapping_bed})

cd ${outputdir}
bedtools intersect -a ${gene_bed} -b ${gro_seq_mapping_bed} -s -wao | awk '$7!=".",OFS="\t" {print $1,$2,$3,$4,$5,$6}' | \
 sort -k 1,1 -k 2n,2 -k 3n,3 | uniq -c | awk 'OFS="\t" {print $2,$3,$4,$5,$6,$7,$1}' >> ${outprefix}.bed
 
 
