#!/bin/bash

# Boyuan_Li

# Wed Dec 18 21:39:29 2019

# this script is used to normalize the pro-seq

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -N <fly_num> -p <prefix> -i <bam_file> -I <spike_in> "
	echo ""
	# echo " -n	string       	[required] the mm10 reads num"
	# echo ""
	echo " -N	string       	[required] the fly reads num"
	echo ""
	echo " -p	string       	[required] the output prefix"
	echo ""
	echo " -i	file         	[required] the input bam file"
	echo ""
	echo " -I	string       	[required] Is there a spike in? YES or NO"
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "N:p:i:I:h" optionName
do
	case $optionName in
		# n) mm10_num="$OPTARG";;
		N) fly_num="$OPTARG";;
		p) prefix="$OPTARG";;
		i) bam_file="$OPTARG";;
		I) spike_in="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


# if [[ $mm10_num = "" ]]; then
	# echo " the $mm10_num string is needed "
	# exit 1
# fi

if [[ $fly_num = "" ]]; then
	echo " the $fly_num string is needed "
	exit 1
fi

if [[ $prefix = "" ]]; then
	echo " the $prefix string is needed "
	exit 1
fi

if [[ $bam_file = "" ]]; then
	echo "the $bam_file file is needed "
	exit 1
elif [[ ! -f $bam_file ]]; then
	echo "$bam_file:   is not found"
	exit 2
fi

if [[ $spike_in = "YES" ]]; then
	echo " You used spike in!!! "
	scale_factor=$(echo 1000000/${fly_num} | bc -l)
	bamCoverage --bam $bam_file -o ${prefix}_forward.bw --binSize 1 -p 10 --extendReads --normalizeUsing None --scaleFactor ${scale_factor} --filterRNAstrand forward --outFileFormat bigwig
	bamCoverage --bam $bam_file -o ${prefix}_reverse.bw --binSize 1 -p 10 --extendReads --normalizeUsing None --scaleFactor -${scale_factor} --filterRNAstrand reverse --outFileFormat bigwig
	
elif [[ $spike_in = "NO" ]]; then
	echo " You didn't use spike in!!! "
	bamCoverage --bam $bam_file -o ${prefix}_forward.bw --binSize 1 -p 10 --extendReads --normalizeUsing RPKM --filterRNAstrand forward --outFileFormat bigwig
	bamCoverage --bam $bam_file -o ${prefix}_reverse.bw --binSize 1 -p 10 --extendReads --normalizeUsing RPKM --scaleFactor -1 --filterRNAstrand reverse --outFileFormat bigwig

else 
	echo " You given a wrong TYPE for -I, please use YES or NO !!! "
fi

