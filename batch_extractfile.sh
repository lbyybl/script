#!/bin/bash

# Boyuan-Li

# Tue Jun 19 22:26:05 CST 2018

# This script is used to batch extract file in one direcotry;

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -i <Input> -o <Output> -n <Number> "
	echo ""
	echo " -i STRING          [required]  input direcotry"
	echo ""
	echo " -o STRING          [required]  output direcotry"
	echo ""
	echo " -n STRING          [required]  how many reads you want to extract"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "i:o:n:h" optionName
do
	case $optionName in
		i) Input="$OPTARG";;
		o) Output="$OPTARG";;
		n) Number="$OPTARG";;
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
	echo "$Output:   is not found"
	exit 2
fi


if [[ $Number = "" ]]; then
	echo " -n the Number STRING is needed "
	exit 1
fi

Abs_output=$(readlink -e $Output)
readlink -e $Input/* | parallel basename {} | sed -e 's/.gz//g' | uniq >> $Input/sample.list
cat $Input/sample.list | parallel head -$Number '<''('less $Input/{}*')' '>>' $Abs_output/{}
rm $Input/sample.list
