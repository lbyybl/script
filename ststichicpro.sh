#!/bin/bash

# Boyuan-Li

# Sat Jun 23 23:06:25 CST 2018

# this script is used to extract the stastic info from the hic-pro result
# can't parallel you only to input the dir of hic-pro result and this script
# will produce a dir _stastic_ and in the sumsta, one file named final2_formal is the 
# stastic file; you can trans it into excel in the excel;

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -i <Input_hicpro_dir> -o <Output_file> "
	echo ""
	echo " -i STRING          [required] the result fo hicpro"
	echo ""
	echo " -o STRING          [required]  Output_file"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "i:o:h" optionName
do
	case $optionName in
		i) Input_hicpro_dir="$OPTARG";;
		o) Output_file="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Input_hicpro_dir = "" ]]; then
	echo " -i the Input_hicpro_dir directory is needed "
	exit 1
elif [[ ! -d $Input_hicpro_dir ]]; then
	echo "$Input_hicpro_dir:   is not found"
	exit 2
fi


if [[ $Output_file = "" ]]; then
	echo " -o the Output_file directory is needed "
	exit 1
#elif [[ ! -d $Output_file ]]; then
#	echo "$Output_file:   is not found"
#	exit 2
fi

Input_file_name=${Input_hicpro_dir##*/}
Output_file_name=${Output_file##*/}

Output_file_abs=$(pwd)/$Output_file
mkdir -p _${Input_file_name}_${Output_file_name}_
ls -d $Input_hicpro_dir/hic_results/data/* | parallel cp {}/*_allValidPairs.mergestat ./_${Input_file_name}_${Output_file_name}_
ls -d $Input_hicpro_dir/bowtie_results/bwt2/* | parallel -v cp {}/*.mpairstat ./_${Input_file_name}_${Output_file_name}_
ls -d $Input_hicpro_dir/hic_results/data/* | parallel -v cp {}/*.mRSstat ./_${Input_file_name}_${Output_file_name}_

cd ./_${Input_file_name}_${Output_file_name}_
mkdir sumsta
sed -i -e 's/#.*//g' -e '/^$/d' *.mRSstat
sed -i -e 's/#.*//g' -e '/^$/d' *.mpairstat

# exchangerow-col.sh -f *.mRSstat >> ./sumsta/filter.stat
# exchangerow-col.sh -f *_allValidPairs.mergestat >> ./sumsta/interaction.stat
# exchangerow-col.sh -f *.mpairstat >> ./sumsta/pair.stat
# awk 'NR<=2{print $0}' ./sumsta/pair.stat >> ./sumsta/pair.state
 
ls *.mRSstat | parallel -j1 -v cat {} >> mRSstat.stats
sed -i -e 's/cat\s/sample_name\t/g' -e 's/.mRSstat.*$//g' mRSstat.stats
split -l 11 mRSstat.stats mRSstatsplit
ls mRSstatsplit* | parallel -j1 exchangerow-col.sh -f {} >> ./sumsta/filter.stat

ls *_allValidPairs.mergestat | parallel -j1 -v cat {} >> allValidPairs.stat
sed -i -e 's/cat\s/sample_name\t/g' -e 's/_allValidPairs.mergestat.*$//g' allValidPairs.stat
split -l 7 allValidPairs.stat allValidPairssplit
ls allValidPairssplit* | parallel -j1 exchangerow-col.sh -f {} >> ./sumsta/interaction.stat

ls *.mpairstat | parallel -j1 -v cat {} >> mpairstat.stat
sed -i -e 's/cat\s/sample_name\t/g' -e 's/.mpairstat.*$//g' mpairstat.stat
cut -f 1-2 mpairstat.stat >> mpairstat.stat2
split -l 11 mpairstat.stat2 mpairstatsplit
ls mpairstatsplit* | parallel -j1 exchangerow-col.sh -f {} >> ./sumsta/pair.state 

cd sumsta

paste pair.state filter.stat >> pair_filter
paste pair_filter interaction.stat | column -t >> final2.stastic 

Rscript /home/boyuanli/bashscript/bin/makecsv.R -i final2.stastic -o $Output_file_abs 

rm -r ${Output_file_abs%/*}/_${Input_file_name}_${Output_file_name}_
