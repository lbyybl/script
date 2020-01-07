#!/bin/bash

# Boyuan_Li

# Fri Oct 11 00:06:25 2019

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -i <input> "
	echo ""
	echo " -i	file         	[required] "
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "i:h" optionName
do
	case $optionName in
		i) input="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $input = "" ]]; then
	echo "the $input file is needed "
	exit 1
elif [[ ! -f $input ]]; then
	echo "$input:   is not found"
	exit 2
fi

input_dir=${input%%/*}
file_name=${input##*/}
file_prefix=${file_name%.*}
out_dir=${input_dir}/unique

sambamba view -h --nthreads 10 -f bam -F "[XS] == null and not unmapped and paired and proper_pair and not duplicate" $input >> ${out_dir}/${file_name}
samtools index ${out_dir}/${file_name}
samtools flagstat ${out_dir}/${file_name} >> ${out_dir}/${file_prefix}.flag
