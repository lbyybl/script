#!/bin/bash

# Boyuan-Li

# Tue Jun 19 21:56:43 CST 2018

# This script used trimlinker to trimlinker, but sometimes deal many files in the same directory, trimlinkdr is inconvenient

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -i <Input> -o <Output> -a <AdaptorA> -b <AdaptorB> "
	echo ""
	echo " -i STRING          [required]  the input directoryy"
	echo ""
	echo " -o STRING          [required]  the output directory"
	echo ""
	echo " -a STRING          [required]  the adaptor A in trimlinker"
	echo ""
	echo " -b STRING          [required]  the adaptor B in trimlinker"
	echo " If you want to adjust other parameter of Trimlinker, please modify this script"
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "i:o:a:b:h" optionName
do
	case $optionName in
		i) Input="$OPTARG";;
		o) Output="$OPTARG";;
		a) AdaptorA="$OPTARG";;
		b) AdaptorB="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Input = "" ]]; then
	echo " -i the Input directory is needed "
	exit 1
elif [[ ! -d $Input ]]; then
	echo "$Input:   is not found"
	exit 2
fi


if [[ $Output = "" ]]; then
	echo " -o the Output directory is needed "
	exit 1
elif [[ ! -d $Output ]]; then
	echo "$Output :   is not found"
	echo "$Output : will be mkdir ......."
	mkdir -p $Output
fi


if [[ $AdaptorA = "" ]]; then
	echo " -a the AdaptorA STRING is needed "
	exit 1
fi


if [[ $AdaptorB = "" ]]; then
	echo " -b the AdaptorB STRING is needed "
	exit 1
fi

Abs_output=$(readlink -e $Output)
readlink -e $Input/* | parallel basename {} | sed -e 's/_R[1|2].*//g' | uniq >> $Input/sample.list
cat $Input/sample.list | parallel trimLinker -t 2 -m 1 -k 1 -l 16 -o $Abs_output -n {} -A $AdaptorA -B $AdaptorB $Input/{}_R1* $Input/{}_R2*

