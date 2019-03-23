#!/bin/bash

# Boyuan-Li

# Wed Jan  2 10:46:47 CST 2019

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -w <Wt> -k <Ko> -o <Out_prefix> "
	echo ""
	echo " -w STRING          [required] WT file"
	echo ""
	echo " -k STRING          [required] KO file"
	echo ""
	echo " -o STRING          [required] Output file prifix"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "w:k:o:h" optionName
do
	case $optionName in
		w) Wt="$OPTARG";;
		k) Ko="$OPTARG";;
		o) Out_prefix="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Wt = "" ]]; then
	echo " -w the Wt file is needed "
	exit 1
elif [[ ! -f $Wt ]]; then
	echo "$Wt:   is not found"
	exit 2
fi


if [[ $Ko = "" ]]; then
	echo " -k the Ko file is needed "
	exit 1
elif [[ ! -f $Ko ]]; then
	echo "$Ko:   is not found"
	exit 2
fi


if [[ $Out_prefix = "" ]]; then
	echo " -o the Out_prefix STRING is needed "
	exit 1
fi

#wt_file=$(readlink -e $Wt)
#ko_file=$(readlink -e $Ko)

wt_file=$(basename $Wt)
ko_file=$(basename $Ko)

#--- select reads long than 20k
	awk '{if (($2!=$5)) print $0}' ${wt_file} >> ${wt_file%.*}_select.pairs
	awk '{if (($2!=$5)) print $0}' ${ko_file} >> ${ko_file%.*}_select.pairs

#---got reads num
	wt_num=$(wc -l ${wt_file%.*}_select.pairs | awk '{print $1}')
	ko_num=$(wc -l ${ko_file%.*}_select.pairs | awk '{print $1}')
#---sample large equal to samll
	if [[ ${wt_num} > ${ko_num} ]]; then
		subsample -s 123456 -n ${ko_num} ${wt_file%.*}_select.pairs >> ${wt_file%.*}_select_sam.pairs
		mv ${wt_file%.*}_select_sam.pairs ${Out_prefix}_wt_sub.allValidPairs
		mv ${ko_file%.*}_select.pairs ${Out_prefix}_ko_sub.allValidPairs
	else
		subsample -s 123456 -n ${wt_num} ${ko_file%.*}_select.pairs >> ${ko_file%.*}_select_sam.pairs
		mv ${ko_file%.*}_select_sam.pairs ${Out_prefix}_ko_sub.allValidPairs
		mv ${wt_file%.*}_select.pairs ${Out_prefix}_wt_sub.allValidPairs
	fi
	
echo -e "${ko_num}\t${wt_num}"
wc -l *_sub.allValidPairs
