#!/bin/bash

# Boyuan-Li

# Wed Nov  7 17:12:41 CST 2018

# It's used to got the diff interaction between wt and ko with .hic file

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -w <Wt_file> -k <Ko_file> -r <Resolution> "
	echo ""
	echo " -w STRING          [required] WT .hic file"
	echo ""
	echo " -k STRING          [required] KO .hic file"
	echo ""
	echo " -r STRING          [required] The resolution to call difference"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "w:k:r:h" optionName
do
	case $optionName in
		w) Wt_file="$OPTARG";;
		k) Ko_file="$OPTARG";;
		r) Resolution="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Wt_file = "" ]]; then
	echo " -w the Wt_file file is needed "
	exit 1
elif [[ ! -f $Wt_file ]]; then
	echo "$Wt_file:   is not found"
	exit 2
fi


if [[ $Ko_file = "" ]]; then
	echo " -k the Ko_file file is needed "
	exit 1
elif [[ ! -f $Ko_file ]]; then
	echo "$Ko_file:   is not found"
	exit 2
fi


if [[ $Resolution = "" ]]; then
	echo " -r the Resolution STRING is needed "
	exit 1
fi

#--- got the prefix 
	filename=${Ko_file##*/}
	prefix=$(echo $filename | sed -e "s/_ko_\|_KO_\|_Ko_\|.hic/_/g")
	

#--- 利用HiCcompare call difference
	echo {{1..19},X,Y} | xargs -n1 | parallel "java -jar /DATA/work/lbyybl/tools/juicer_fml/juicer_tools.1.8.9_jcuda.0.8.jar dump observed NONE $Ko_file {} {} BP $Resolution ${prefix}Ko_{.}.observed"
	echo {{1..19},X,Y} | xargs -n1 | parallel "java -jar /DATA/work/lbyybl/tools/juicer_fml/juicer_tools.1.8.9_jcuda.0.8.jar dump observed NONE $Wt_file {} {} BP $Resolution ${prefix}Wt_{.}.observed"
	echo {{1..19},X,Y} | xargs -n1 | parallel -j10 "/usr/bin/Rscript /home/boyuanli/bashscript/bin/HiCcompare/hicdifference.R -w ${prefix}Wt_{}.observed -k ${prefix}Ko_{}.observed -o ${prefix}chr{}.diff -c chr{}"


