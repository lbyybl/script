#!/bin/bash

# Boyuan-Li

# Mon Jul 30 22:29:06 CST 2018

# This stript is used to got the tad in suitable range for the metatad

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -s <St_bin> -e <En_bin> "
	echo ""
	echo " -s STRING          [required] start bin"
	echo ""
	echo " -e STRING          [required] end bin"
	echo " -e STRING          [required] input file"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "s:e:i:h" optionName
do
	case $optionName in
		s) St_bin="$OPTARG";;
		e) En_bin="$OPTARG";;
		i) Inputfile="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $St_bin = "" ]]; then
	echo " -s the St_bin STRING is needed "
	exit 1
fi


if [[ $En_bin = "" ]]; then
	echo " -e the En_bin STRING is needed "
	exit 1
fi

awk -v st_bin="${St_bin}" -v en_bin="${En_bin}" '{if ($2 >= 2*st_bin-en_bin && $2 <= 2*en_bin-st_bin && $1>= 2*st_bin-en_bin && $1 <= 2*en_bin-st_bin) print $0}' ${Inputfile}
