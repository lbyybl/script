#!/bin/bash

# Boyuan_Li

# Tue Sep 24 16:54:17 2019

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -f <peak1_file> -p <peak2_file> -o <output_mergefile> "
	echo ""
	echo " -f	file         	[required] "
	echo ""
	echo " -p	file         	[required] "
	echo ""
	echo " -o	string       	[required] "
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "f:p:o:h" optionName
do
	case $optionName in
		f) peak1_file="$OPTARG";;
		p) peak2_file="$OPTARG";;
		o) output_mergefile="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $peak1_file = "" ]]; then
	echo "the $peak1_file file is needed "
	exit 1
elif [[ ! -f $peak1_file ]]; then
	echo "$peak1_file:   is not found"
	exit 2
fi

if [[ $peak2_file = "" ]]; then
	echo "the $peak2_file file is needed "
	exit 1
elif [[ ! -f $peak2_file ]]; then
	echo "$peak2_file:   is not found"
	exit 2
fi

if [[ $output_mergefile = "" ]]; then
	echo " the $output_mergefile string is needed "
	exit 1
fi

peak1=$peak1_file
peak2=$peak2_file
out_file=$output_mergefile
	bedtools intersect -a $peak1 -b $peak2 -wa -wb \
	| awk 'OFS="\t"{print $1,$2,$3,$10,$11,$12}' | awk 'OFS="\t"{print $1,$2,$3"\n"$4,$5,$6}' >> $out_file
	sort -k 1,1 -k 2n,2 -k 3n,3 $out_file -o $out_file
	bedtools merge -i $out_file | sort -k 1,1 -k 2n,2 -k 3n,3 -o $out_file
