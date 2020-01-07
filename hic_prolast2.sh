#!/bin/bash

# Boyuan-Li

# Sat Nov 24 13:44:30 CST 2018

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -i <Input_dir_include_allvalidparis> -o <Output_dir> -c <Config_file> "
	echo ""
	echo " -i STRING          [required] two depth"
	echo ""
	echo " -o STRING          [required] "
	echo ""
	echo " -c STRING          [required] "
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "i:o:c:h" optionName
do
	case $optionName in
		i) Input_dir_include_allvalidparis="$OPTARG";;
		o) Output_dir="$OPTARG";;
		c) Config_file="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Input_dir_include_allvalidparis = "" ]]; then
	echo " -i the Input_dir_include_allvalidparis directory is needed "
	exit 1
elif [[ ! -d $Input_dir_include_allvalidparis ]]; then
	echo "$Input_dir_include_allvalidparis:   is not found"
	exit 2
fi


if [[ $Output_dir = "" ]]; then
	echo " -o the Output_dir STRING is needed "
	exit 1
fi


if [[ $Config_file = "" ]]; then
	echo " -c the Config_file file is needed "
	exit 1
elif [[ ! -f $Config_file ]]; then
	echo "$Config_file:   is not found"
	exit 2
fi


HiC-Pro -i $Input_dir_include_allvalidparis -o $Output_dir -c $Config_file -s build_contact_maps
HiC-Pro -i $Output_dir/hic_results/matrix/ -o $Output_dir -c $Config_file -s ice_norm

