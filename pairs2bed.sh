#!/bin/bash

# Boyuan_Li

# Tue Feb 19 19:29:55 2019

# it's used to convert the output of hicpro into bigwig

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -i <input_dir> -o <output_dir> -p <prefix> "
	echo ""
	echo " -i	dir          	[required] "
	echo ""
	echo " -o	dir          	[required] "
	echo ""
	echo " -p	string       	[required] "
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "i:o:p:h" optionName
do
	case $optionName in
		i) input_dir="$OPTARG";;
		o) output_dir="$OPTARG";;
		p) prefix="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $input_dir = "" ]]; then
	echo "the $input_dir file is needed "
	exit 1
elif [[ ! -d $input_dir ]]; then
	 echo "$input_dir:   is not found"
	exit 2
fi

if [[ $output_dir = "" ]]; then
	echo "the $output_dir file is needed "
	exit 1
elif [[ ! -d $output_dir ]]; then
	 echo "$output_dir:   is not found"
	exit 2
fi

if [[ $prefix = "" ]]; then
	echo " the $prefix string is needed "
	exit 1
fi

/DATA/work/lbyybl/tools/bedItemOverlapCount/pairs2bed.py -i ${input_dir} -o ${output_dir}
rm -f ${output_dir}/all.Pairs.tmp
cat ${output_dir}/allpairs.bed.tmp | egrep -w "chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chrX|chrY"  | sort -k 1,1 -k 2n,2| /DATA/work/lbyybl/tools/bedItemOverlapCount/bedItemOverlapCount mm10 -chromSize=/DATA/work/lbyybl/tools/bedItemOverlapCount/mm10.chrom.sizes stdin | sort -k1,1 -k2,2n >> ${output_dir}/allpairs.bedgraph
rm -f ${output_dir}/allpairs.bed.tmp ${output_dir}/all_select.Pairs.tmp ${output_dir}/all_select.Pairs.tmp
cd ${output_dir}/
/DATA/work/lbyybl/tools/bedItemOverlapCount/bedGraphToBigWig allpairs.bedgraph /DATA/work/lbyybl/tools/bedItemOverlapCount/mm10.chrom.sizes ${prefix}.bw
#rm -f ${output_dir}/allpairs.bed.tmp ${output_dir}/allpairs.bedgraph ${output_dir}/all.Pairs.tmp ${output_dir}/allpairs.bed
