#!/bin/bash

# Boyuan-Li

# Sun Jun 17 22:35:53 CST 2018

# This script is used to trans the row into col

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -f <Input file> -s <Delimiter> "
	echo ""
	echo " -f STRING          [required] Input file that you want to trans"
	echo ""
	echo " -s STRING          [required] Delimiter you can assign, defult '\t'"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}

Delimiter="\t"

if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "f:s:h" optionName
do
	case $optionName in
		f) Input="$OPTARG";;
		s) Delimiter="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Input = "" ]]; then
	echo " -f the Input file is needed "
	exit 1
elif [[ ! -f $Input ]]; then
	echo "$Input:   is not found"
	exit 2
fi


if [[ $Delimiter = "" ]]; then
	echo " -s the Delimiter STRING is needed "
	exit 1
fi

awk -F "$Delimiter" '{ for(i=1;i<=length($0);i++)
        col[i]=col[i]"\t"$i
    }
    END { i=1
        while (i in col)
        {
            print col[i]"\t"
            i++
        }
    }'  $Input | column -t 
