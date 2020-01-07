#!/bin/bash

# Boyuan-Li

# Fri Jan  4 10:14:08 CST 2019

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -g <Gene_bedfile> -r <Resolution> -o <Output_dir> "
	echo ""
	echo " -g STRING          [required] Ucsc refseq file"
	echo ""
	echo " -r STRING          [Optional] defult 2000"
	echo ""
	echo " -o STRING          [required] Output dir"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
Resolution=2000
while getopts "g:r:o:h" optionName
do
	case $optionName in
		g) Gene_bedfile="$OPTARG";;
		r) Resolution="$OPTARG";;
		o) Output_dir="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Gene_bedfile = "" ]]; then
	echo " -g the Gene_bedfile file is needed "
	exit 1
elif [[ ! -f $Gene_bedfile ]]; then
	echo "$Gene_bedfile:   is not found"
	exit 2
fi


if [[ $Resolution = "" ]]; then
	echo " -r the Resolution STRING is needed "
	exit 1
fi


if [[ $Output_dir = "" ]]; then
	echo " -o the Output_dir directory is needed "
	exit 1
elif [[ ! -d $Output_dir ]]; then
	echo "$Output_dir:   is not found"
	exit 2
fi

#---- distanguish +- strand
awk -v res=${Resolution} '$6=="+", OFS="\t" {st=$2-res;en=$2+res;print $1,st,en,$4,$5,$6}' $Gene_bedfile | sort -k 1,1 -k 2n,2 -k 3n,3 | uniq | awk 'OFS="\t" {cnt++;print $1,$2,$3,$4,$5,$6}'>> ${Output_dir}/Ppromoter_${Resolution}.bed
awk -v res=${Resolution} '$6=="-", OFS="\t" {st=$3-res;en=$3+res;print $1,st,en,$4,$5,$6}' $Gene_bedfile | sort -k 1,1 -k 2n,2 -k 3n,3 | uniq | awk 'OFS="\t" {cnt++;print $1,$2,$3,$4,$5,$6}'>> ${Output_dir}/Npromoter_${Resolution}.bed
cat ${Output_dir}/Ppromoter_${Resolution}.bed ${Output_dir}/Npromoter_${Resolution}.bed >> ${Output_dir}/promoter_${Resolution}.bed
chr=$(echo chr{{1..19},X,Y} | sed 's/ /|/g')
egrep -w "$chr" ${Output_dir}/promoter_${Resolution}.bed >> ${Output_dir}/promoter_final_${Resolution}.bed
rm -f ${Output_dir}/promoter_${Resolution}.bed ${Output_dir}/Ppromoter_${Resolution}.bed ${Output_dir}/Npromoter_${Resolution}.bed
#bedtools merge -i <(sort -k 1,1 -k 2n,2 -k 3n,3 ${Output_dir}/promoter_final_${Resolution}.bed) -c 4,5 -o collapse,collapse >> ${Output_dir}/promoter_merge_${Resolution}.bed
