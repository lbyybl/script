#!/bin/bash

# Boyuan-Li

# Sat Sep 22 10:16:43 CST 2018

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -n <Namelist2_file> -o <Output_file> -i <Input_bedfile> "
	echo ""
	echo " -n STRING          [required] "
	echo ""
	echo " -o STRING          [required] "
	echo ""
	echo " -i STRING          [required] "
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "n:o:i:h" optionName
do
	case $optionName in
		n) Namelist2_file="$OPTARG";;
		o) Output_file="$OPTARG";;
		i) Input_bedfile="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Namelist2_file = "" ]]; then
	echo " -n the Namelist2_file file is needed "
	exit 1
elif [[ ! -f $Namelist2_file ]]; then
	echo "$Namelist2_file:   is not found"
	exit 2
fi


if [[ $Output_file = "" ]]; then
	echo " -o the Output_file STRING is needed "
	exit 1
fi


if [[ $Input_bedfile = "" ]]; then
	echo " -i the Input_bedfile file is needed "
	exit 1
elif [[ ! -f $Input_bedfile ]]; then
	echo "$Input_bedfile:   is not found"
	exit 2
fi
Prefix=${Output_file%.*}
cat $Namelist2_file | parallel -C '\t' " egrep \"{1}\"  $Input_bedfile >> _${Prefix}_{2}_-_{1}.bed "
ls _${Prefix}*_-_*.bed | parallel -k -v "Rscript /home/boyuanli/bashscript/bin/geneinteraction/intervalsum.r -i {}" >> _${Prefix}intersum.list
sed -e "s/Rscript.*_${Prefix}_//g" -e 's/.bed//g' -e 's/\[1\]//g' -e 's/_-_.*$//g' _${Prefix}intersum.list | xargs -l2 | awk '{printf "%s\t%s\n",$1,$2}' | sort -o $Output_file
#--- 将 EN IN的区间分别除二
awk '{if ($1=="SE" || $1=="TE" || $1=="IN") printf "%s\t%s\n",$1,$2/2; else print $0}' $Output_file | sort -o $Output_file
rm _${Prefix}*_-_*.bed _${Prefix}intersum.list
