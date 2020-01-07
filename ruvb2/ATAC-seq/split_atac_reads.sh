#!/bin/bash

# Boyuan_Li

# Wed May  8 13:46:46 2019

# This pipeline is used to split the ATAC-READS into <115 bp and 180~247 bp two class.

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -i <Input_fie> -o <Output_dir> -p <prefix> "
	echo ""
	echo " -i	string       	[required] the bam file you want to split"
	echo ""
	echo " -o	dir          	[required] the output dir you want to place the split file"
	echo ""
	echo " -p	string       	[required] the prefix you want to name"
	echo ""
	echo " -h			help"
	echo ""
	echo -e "\tdefault the specise is mm10,if not please modify the scripts"
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "i:o:p:h" optionName
do
	case $optionName in
		i) Input_fie="$OPTARG";;
		o) Output_dir="$OPTARG";;
		p) prefix="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $Input_fie = "" ]]; then
	echo " the $Input_fie string is needed "
	exit 1
fi

if [[ $Output_dir = "" ]]; then
	echo "the $Output_dir file is needed "
	exit 1
elif [[ ! -d $Output_dir ]]; then
	 echo "$Output_dir:   is not found"
	exit 2
fi

if [[ $prefix = "" ]]; then
	echo " the $prefix string is needed "
	exit 1
fi
# blacklist=/DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed
# #--- <=115
# samtools view -H $Input_fie >> ${Output_dir}/${prefix}_free.sam
# samtools view $Input_fie | awk '{if ($7=="=") ch=$8-$4; if (ch>0 && ch < 116) print $0; else if (ch <0 && ch > -116) print $0}' >> ${Output_dir}/${prefix}_free.sam
# samtools sort -@ 3 ${Output_dir}/${prefix}_free.sam >> ${Output_dir}/${prefix}_free.bam 
# samtools index ${Output_dir}/${prefix}_free.bam 
# rm -f ${Output_dir}/${prefix}_free.sam &
# bamCoverage --bam ${Output_dir}/${prefix}_free.bam  -o ${Output_dir}/${prefix}_free_rpgc.bw --binSize 10 --extendReads 300 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --blackListFileName ${blacklist}
#--- 180-247
samtools view -H $Input_fie >> ${Output_dir}/${prefix}_nucleo.sam
samtools view $Input_fie | awk '{if ($7=="=") ch=$8-$4; if (ch>0 && ch < 248 && ch > 179) print $0; else if (ch <0 && ch > -248 && ch < -179) print $0}' >> ${Output_dir}/${prefix}_nucleo.sam
samtools sort -@ 3 ${Output_dir}/${prefix}_nucleo.sam >> ${Output_dir}/${prefix}_nucleo.bam 
rm -f ${Output_dir}/${prefix}_nucleo.sam &
