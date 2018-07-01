#!/bin/bash

# Boyuan-Li

# Tue Jun 19 21:43:14 CST 2018

# This scipt is used to produce link file, for the directory

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -i <Input> -o <Output> "
	echo ""
	echo " -i STRING          [required] the input should be directory"
	echo ""
	echo " -o STRING          [required] the output direcoty can unexit"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "i:o:h" optionName
do
	case $optionName in
		i) Input="$OPTARG";;
		o) Output="$OPTARG";;
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
	echo "whill mkdir $Output......"
	mkdir -p $Output
fi
Abs_output=$(readlink -e $Output)
readlink -e $Input/* | parallel ln {} $Abs_output
