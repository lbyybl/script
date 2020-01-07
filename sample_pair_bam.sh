#!/bin/bash

# Boyuan_Li

# Tue Oct 15 23:05:28 2019

# This file is used to sample down Paired-end BAM file according to the number you given

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -i <Input> -o <ouput_file> -n <number> -t <genome is mix or mm10>"
	echo ""
	echo " -i	file         	[required] Input file you want to sample down"
	echo ""
	echo " -o	string       	[required] Output file name"
	echo ""
	echo " -n	string       	[required] The final reads number you want to down sample; If it's a mix genome, this is the reads of fly"
	echo ""
	echo " -t	string       	[optional] if 'mix'; The genome is mix of mm10 and fly and normalized by fly reads"
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}
type=""

if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "i:o:n:t:h" optionName
do
	case $optionName in
		i) Input="$OPTARG";;
		o) ouput_file="$OPTARG";;
		n) number="$OPTARG";;
		t) type="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $Input = "" ]]; then
	echo "the $Input file is needed "
	exit 1
elif [[ ! -f $Input ]]; then
	echo "$Input:   is not found"
	exit 2
fi

if [[ $ouput_file = "" ]]; then
	echo " the $ouput_file string is needed "
	exit 1
fi

if [[ $number = "" ]]; then
	echo " the $number string is needed "
	exit 1
fi
if [[ $type = "" ]]; then
	echo " it is not a mix genome "
	reads_n=$(samtools view $Input | wc -l)
	probability=$(echo ${number}/${reads_n} | bc -l)
	/usr/bin/java -jar /DATA/software/bcbio/anaconda/share/picard-2.17.0-0/picard.jar DownsampleSam I=${Input} O=${ouput_file} RANDOM_SEED=123 PROBABILITY=${probability}
elif [[ $type = "mix" ]]; then
	echo " Note!!! it is a mix genome "
	reads_n=$(samtools view $Input | grep "dm6_chr" | wc -l)
	probability=$(echo ${number}/${reads_n} | bc -l)
	/usr/bin/java -jar /DATA/software/bcbio/anaconda/share/picard-2.17.0-0/picard.jar DownsampleSam I=${Input} O=${ouput_file} RANDOM_SEED=123 PROBABILITY=${probability}
else
	echo "You gave the wrong type!!!"
fi

