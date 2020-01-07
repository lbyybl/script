#!/bin/bash

# Boyuan_Li

# Wed Nov 13 17:30:30 2019

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -d <work_dir> "
	echo ""
	echo " -d	dir          	[required] "
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "d:h" optionName
do
	case $optionName in
		d) work_dir="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $work_dir = "" ]]; then
	echo "the $work_dir file is needed "
	exit 1
elif [[ ! -d $work_dir ]]; then
	 echo "$work_dir:   is not found"
	exit 2
fi

for i in `ls ${work_dir}/*.bed`; do
	n=$(wc -l ${i} | awk '{print $1}');
	if [[ $n = 0 ]];then
		#region=$region $i;
		echo -e "$i will be remove, because it's null!";
		rm -f $i
	fi
done