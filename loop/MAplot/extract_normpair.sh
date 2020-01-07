#!/bin/bash

# Boyuan_Li

# Wed Jan  9 12:47:07 2019

# This is used to extract correspording bin pair interaction from .hic file

loop_bin_file=loop_pair_bin.bedpe
reso=5000
set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -j <hic_file> -r <reso> -l <loop_bin_file> -p <output_prefix> "
	echo ""
	echo " -j	file         	[required] hic file"
	echo ""
	echo " -r	string       	[required] resolutin default 5000"
	echo ""
	echo " -l	file         	[required] loop_bin_file;default loop_pair_bin.bedpe"
	echo ""
	echo " -p	string       	[required] output_prefix"
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "j:r:l:p:h" optionName
do
	case $optionName in
		j) hic_file="$OPTARG";;
		r) reso="$OPTARG";;
		l) loop_bin_file="$OPTARG";;
		p) output_prefix="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $hic_file = "" ]]; then
	echo "the $hic_file file is needed "
	exit 1
elif [[ ! -f $hic_file ]]; then
	echo "$hic_file:   is not found"
	exit 2
fi

if [[ $reso = "" ]]; then
	echo " the $reso string is needed "
	exit 1
fi

if [[ $loop_bin_file = "" ]]; then
	echo "the $loop_bin_file file is needed "
	exit 1
elif [[ ! -f $loop_bin_file ]]; then
	echo "$loop_bin_file:   is not found"
	exit 2
fi

if [[ $output_prefix = "" ]]; then
	echo " the $output_prefix string is needed "
	exit 1
fi

out_put_dir=tmp
mkdir -p ${output_prefix}_${out_put_dir}
echo chr{{1..19},X,Y} | xargs -n1 | parallel "java -jar /DATA/work/lbyybl/tools/juicer_fml/juicer_tools.1.8.9_jcuda.0.8.jar dump oe \
VC_SQRT ${hic_file} {} {} BP ${reso}  | awk -v chr={} 'OFS=\"\t\"{print chr\":\"\$1\":\"\$2,\$3}'>> ${output_prefix}_${out_put_dir}/${output_prefix}_{}_juicer_dump.bed"

awk 'OFS="\t" {print $1":"$2":"$5,$1}' ${loop_bin_file} | parallel -C "\t" -k "egrep -w '{1}' ${output_prefix}_${out_put_dir}/${output_prefix}_{2}_juicer_dump.bed" >> ${output_prefix}_loop.bed

rm -rf ${output_prefix}_${out_put_dir}
