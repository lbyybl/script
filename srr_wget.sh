#!/bin/bash

# Boyuan_Li

# Mon Apr  8 11:23:21 2019

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -s <Srr_entriy> "
	echo ""
	echo " -s	string       	[required] "
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "s:h" optionName
do
	case $optionName in
		s) Srr_entriy="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $Srr_entriy = "" ]]; then
	echo " the $Srr_entriy string is needed "
	exit 1
fi

wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${Srr_entriy:0:6}/${Srr_entriy}/${Srr_entriy}.sra
