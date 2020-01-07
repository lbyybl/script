#!/bin/bash

# Boyuan-Li

# Fri Sep  7 22:22:52 CST 2018

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -i <Input> -s <String> -o <Output> "
	echo ""
	echo " -i STRING          [required] "
	echo ""
	echo " -s STRING          [required] "
	echo ""
	echo " -o STRING          [required] "
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "i:s:o:h" optionName
do
	case $optionName in
		i) Input="$OPTARG";;
		s) String="$OPTARG";;
		o) Output="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Input = "" ]]; then
	echo " -i the Input file is needed "
	exit 1
elif [[ ! -f $Input ]]; then
	echo "$Input:   is not found"
	exit 2
fi


if [[ $String = "" ]]; then
	echo " -s the String STRING is needed "
	exit 1
fi


if [[ $Output = "" ]]; then
	echo " -o the Output STRING is needed "
	exit 1
fi
Old_string1="N"${String}
Old_string2="P"${String}
New_string="T"${String}
#cat $Input | sed -e "s/${Old_string1}/${New_string}/g" -e "s/${Old_string2}/${New_string}/g" | sort | uniq -c -w 9 | awk '{printf "%s\t%s\t%s\n",$2,$3,$4*$1}' | sort -o $Output
cat $Input | sed -e "s/${Old_string1}/${New_string}/g" -e "s/${Old_string2}/${New_string}/g" | sort -o ___$Output
Rscript /home/boyuanli/bashscript/bin/geneinteraction/mergesum.r -i ___$Output -o $Output
rm ___$Output