#!/bin/bash

# Boyuan_Li

# Fri Nov  1 11:37:51 2019

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -f <first> -s <second> -t <third> -f <foth> -F <fifth> -p <prefix> "
	echo ""
	echo " -f	file         	[required] "
	echo ""
	echo " -s	file         	[required] "
	echo ""
	echo " -t	file         	[required] "
	echo ""
	echo " -f	file         	[required] "
	echo ""
	echo " -F	file         	[required] "
	echo ""
	echo " -p	string       	[required] "
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}

first=''
second=''
third=''
foth=''
fifth=''

if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "f:s:t:f:F:p:h" optionName
do
	case $optionName in
		f) first="$OPTARG";;
		s) second="$OPTARG";;
		t) third="$OPTARG";;
		f) foth="$OPTARG";;
		F) fifth="$OPTARG";;
		p) prefix="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done

mkdir ${prefix}_tmp

if [[ $first != "" ]]; then
	echo "the $first file was given "
	#exit 1
	if [[ ! -f $first ]]; then
		echo "$first:   is not found"
		exit 2
	else
		cp $first ${prefix}_tmp/file1.narrowPeak
	fi
fi

if [[ $second != "" ]]; then
	echo "the $second file was given "
	#exit 1
	if [[ ! -f $second ]]; then
		echo "$second:   is not found"
		exit 2
	else
		cp $second ${prefix}_tmp/file2.narrowPeak
	fi
fi

if [[ $third != "" ]]; then
	echo "the $third file was given "
	#exit 1
	if [[ ! -f $third ]]; then
		echo "$third:   is not found"
		exit 2
	else
		cp $third ${prefix}_tmp/file3.narrowPeak
	fi
fi

if [[ $foth != "" ]]; then
	echo "the $foth file was given "
	#exit 1
	if [[ ! -f $foth ]]; then
		echo "$foth:   is not found"
		exit 2
	else
		cp $foth ${prefix}_tmp/file4.narrowPeak
	fi
fi

if [[ $fifth != "" ]]; then
	echo "the $fifth file was given "
	#exit 1
	if [[ ! -f $fifth ]]; then
		echo "$fifth:   is not found"
		exit 2
	else
		cp $fifth ${prefix}_tmp/file5.narrowPeak
	fi
fi

if [[ $prefix = "" ]]; then
	echo " the $prefix string is needed "
	exit 1
fi



	parallel --link "idr --samples ${prefix}_tmp/file{1}.narrowPeak ${prefix}_tmp/file{2}.narrowPeak --input-file-type narrowPeak --output-file ${prefix}_tmp/merge_sub{1}{2}_peak.narrowPeak" ::: 1 1 2 ::: 2 3 3
	ls ${prefix}_tmp/file*.narrowPeak | parallel "bedtools intersect -a {} -b /DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed -v | sort -o {}"
	ls ${prefix}_tmp/file*.narrowPeak | parallel "cut -f 1-10 {} | sort -o {}"
	cat ${prefix}_tmp/file*.narrowPeak > ${prefix}_tmp/merge_all_peak.narrowPeak
	

cp $2 ${prefix}_tmp/merge_ovel1_peak.narrowPeak
#file_num=$(echo $#/2 -1 | bc -l)



for ((i=4;i<=$#-2;i=i+2));
	do 
	n=$(echo "scale=0; ${i}/2" | bc -l)
	p=$(echo "scale=0; ${n} - 1" | bc -l)
	idr --samples ${prefix}_tmp/merge_ovel${p}_peak.narrowPeak ${prefix}_tmp/file${n}.narrowPeak --input-file-type narrowPeak --output-file ${prefix}_tmp/merge_ovel${n}_peak.narrowPeak
	
	#idr --samples merge_ove12_peak.narrowPeak peak/R2doxPolII_overlap_peak.narrowPeak  --input-file-type narrowPeak --output-file merge_ove_peak.narrowPeak
	
done
	
	bedtools intersect -a ${prefix}_tmp/merge_all_peak.narrowPeak -b ${prefix}_tmp/merge_ovel${n}_peak.narrowPeak ${prefix}_tmp/merge_sub*.narrowPeak  -v > ${prefix}_tmp/merge_final_peak.narrowPeak
	echo -e "12\n13\n23" | parallel -j 1 "bedtools intersect -a ${prefix}_tmp/merge_sub{}_peak.narrowPeak -b ${prefix}_tmp/merge_ovel${n}_peak.narrowPeak -v >> ${prefix}_tmp/merge_final_peak.narrowPeak"
	cut -f 1-10 ${prefix}_tmp/merge_ovel${n}_peak.narrowPeak >> ${prefix}_tmp/merge_final_peak.narrowPeak
	awk '{st=$2+$10;en=st+1;print $1"\t"st"\t"en}' ${prefix}_tmp/merge_final_peak.narrowPeak > ${prefix}_tmp/merge_all_peak.summit
	cut -f 1-10 ${prefix}_tmp/merge_final_peak.narrowPeak > ${prefix}_union_merge.narrowPeak
	bedtools merge -i  <(sort -k 1,1 -k 2n,2 -k 3n,3 ${prefix}_union_merge.narrowPeak) -c 4,5,6,7,8,9,10 -o distinct,max,distinct,max,max,max,mean | awk '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%i\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | sort -k 1,1 -k 2n,2 -k 3n,3 -o ${prefix}_union_merge.narrowPeak
	awk '{st=$2+$10;en=st+1;print $1"\t"st"\t"en}' ${prefix}_union_merge.narrowPeak > ${prefix}_union_merge.summit
	#echo $i
	rm -rf ${prefix}_tmp