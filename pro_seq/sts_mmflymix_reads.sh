#!/bin/bash

# Boyuan_Li

# Tue Feb 19 16:05:23 2019

# 写一个能够mapping，转成bigwig文件，还能去接头，call peak的脚本，这个脚本还应该有模块化的功能比如，去接头与质控，mapping与转bigwig，call peak就像hicpro一样；

# 1. 质控 2.去接头 3. mapping 4. 转bigwig 5.call peak 

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -p <Prefix> -o <outputdir>"
	echo ""
	echo " -p	string       	[required] output prefix"
	echo ""
	echo " -o	dir          	[required] output dir"
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}

min_length=25
index=/DATA/software/bcbio/genomes/Mmusculus/mm10/bowtie2/mm10
normalize_method=RPKM
blacklist=/DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed

if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "p:o:h" optionName
do
	case $optionName in
		p) Prefix="$OPTARG";;
		o) outputdir="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $Prefix = "" ]]; then
	echo " the $Prefix string is needed "
	exit 1
fi

if [[ $outputdir = "" ]]; then
	echo "the $outputdir file is needed "
	exit 1
elif [[ ! -d $outputdir ]]; then
	 echo "$outputdir:   is not found"
	exit 2
fi



#--- stastics
#--- total reads
if [[ ! -f ${outputdir}/${Prefix}_mmflystastic.csv ]]; then
	pwd
	#--- mm10
	mm_mapping=$(samtools view -F 4 ${outputdir}/mapping/${Prefix}.bam | grep -v "dm6_chr" | wc -l)
	mm_rmdup=$(samtools view -F 4 ${outputdir}/mapping/${Prefix}_rmdup.bam | grep -v "dm6_chr" | wc -l)
	mm_final_uniq=$(samtools view -F 4 ${outputdir}/unique/${Prefix}_rmdup_uniqe.bam | grep -v "dm6_chr" | wc -l)
	#--- fly
	dm6_mapping=$(samtools view -F 4 ${outputdir}/mapping/${Prefix}.bam | grep "dm6_chr" | wc -l)
	dm6_rmdup=$(samtools view -F 4 ${outputdir}/mapping/${Prefix}_rmdup.bam | grep "dm6_chr" | wc -l)
	dm6_final_uniq=$(samtools view -F 4 ${outputdir}/unique/${Prefix}_rmdup_uniqe.bam | grep "dm6_chr" | wc -l)

	#---- stastic file
	echo -e "sample_name,mm10 mapping reads,mm10 rm duplicate,mm10 final unique,fly mapping reads,fly rm duplicate,fly final unique" >> ${outputdir}/${Prefix}_mmflystastic.csv
	echo -e "${Prefix},${mm_mapping},${mm_rmdup},${mm_final_uniq},${dm6_mapping},${dm6_rmdup},${dm6_final_uniq}" >> ${outputdir}/${Prefix}_mmflystastic.csv
fi
