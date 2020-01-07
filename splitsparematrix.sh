#!/bin/bash

# Boyuan-Li

# Mon Jul 16 23:48:14 CST 2018

# Used to extract submatrix and sub bed

set -eu
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -f <Input_matrix> -c <Chr_you_want_extract> -o <Output_file> -b <Corresponding_bed_file> "
	echo ""
	echo " -f STRING          [required] Sparse matrix by hic-pro"
	echo ""
	echo " -c STRING          [required] chr number"
	echo ""
	echo " -o STRING          [required] Output prefix"
	echo ""
	echo " -b STRING          [required] Corresponding_bed_file"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "f:c:o:b:h" optionName
do
	case $optionName in
		f) Input_matrix="$OPTARG";;
		c) Chr_you_want_extract="$OPTARG";;
		o) Output_file="$OPTARG";;
		b) Corresponding_bed_file="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Input_matrix = "" ]]; then
	echo " -f the Input_matrix file is needed "
	exit 1
elif [[ ! -f $Input_matrix ]]; then
	echo "$Input_matrix:   is not found"
	exit 2
fi


if [[ $Chr_you_want_extract = "" ]]; then
	echo " -c the Chr_you_want_extract STRING is needed "
	exit 1
fi


if [[ $Output_file = "" ]]; then
	echo " -o the Output_file STRING is needed "
	exit 1
fi


if [[ $Corresponding_bed_file = "" ]]; then
	echo " -b the Corresponding_bed_file file is needed "
	exit 1
elif [[ ! -f $Corresponding_bed_file ]]; then
	echo "$Corresponding_bed_file:   is not found"
	exit 2
fi
	egrep -w "$Chr_you_want_extract" $Corresponding_bed_file >> ${Output_file}.bed
	First_bin_number=$(egrep -w "$Chr_you_want_extract" $Corresponding_bed_file | head -1  | cut -f 4)
	Last_bin_number=$(egrep -w "$Chr_you_want_extract" $Corresponding_bed_file | tail -1  | cut -f 4)
	awk -v First="${First_bin_number}" -v Last="${Last_bin_number}" '{if ($1 < Last && $1 > First && $2 < Last && $1> First) print $0}' $Input_matrix >> ${Output_file}.matrix
