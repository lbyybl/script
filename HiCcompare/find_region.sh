#!/bin/bash

# Boyuan-Li

# Thu Nov  8 16:18:59 CST 2018

# this to find the gene or gene and regulation elements interaction with the diff loop callel by hiccompare
# and don't distanguish the strand

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -g <Genome_Bed_file> -d <Diff_loop_fill> -o <Output_file> "
	echo ""
	echo " -g STRING          [required] Gene_bed_file can be file merged by promoter enhancer gff file"
	echo ""
	echo " -d STRING          [required] Diff_loop_file called by HiCcompare"
	echo ""
	echo " -o STRING          [required] Output_file name and the prefix also will be used and the dir will be used as output dir"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "g:d:o:h" optionName
do
	case $optionName in
		g) Genome_Bed_file="$OPTARG";;
		d) Diff_loop_fill="$OPTARG";;
		o) Output_file="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Genome_Bed_file = "" ]]; then
	echo " -g the Genome_Bed_file file is needed "
	exit 1
elif [[ ! -f $Genome_Bed_file ]]; then
	echo "$Genome_Bed_file:   is not found"
	exit 2
fi


if [[ $Diff_loop_fill = "" ]]; then
	echo " -d the Diff_loop_fill file is needed "
	exit 1
elif [[ ! -f $Diff_loop_fill ]]; then
	echo "$Diff_loop_fill:   is not found"
	exit 2
fi


if [[ $Output_file = "" ]]; then
	echo " -o the Output_file STRING is needed "
	exit 1
fi

Gene_file=$(readlink -e $Genome_Bed_file)
Diff_file=$(readlink -e $Diff_loop_fill)
Output_file_name=${Output_file##*/}
Output_dir="./"
Prefix=${Output_file_name%.*}
# dir=""
# dir=$(echo -e "${Output_file}" | grep "/")
# if [[ $dir != "" && $dir != "./" ]]; then
	# Output_dir=${Output_file%/*}
	# mkdir -p ${Output_dir}
# fi	
Output_dir=$(readlink -e ${Output_dir})

#--- sort input
	cp ${Gene_file} ${Output_dir}/${Prefix}_${Gene_file##*/}
	sort -k 1,1 -k 2n,2 -k 3n,3 ${Output_dir}/${Prefix}_${Gene_file##*/} -o ${Output_dir}/${Prefix}_${Gene_file##*/}
	sort -k 1,1 -k 2n,2 -k 3n,3 ${Diff_file} -o ${Diff_file}

	#--- got the reads corresponding gene
		bedtools intersect -a ${Diff_file} -b ${Output_dir}/${Prefix}_${Gene_file##*/} -sorted -wao -sorted -wao | \
		awk '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}' \
		| sort -k 1,1 -k 2,2n -k 3,3n | bedtools intersect -a stdin -b ${Output_dir}/${Prefix}_${Gene_file##*/} -sorted -wao | \
		awk '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29}' \
		>> ${Output_dir}/$Output_file_name &
		PS5=$!
		
		wait $PS5 
	
		rm ${Output_dir}/${Prefix}_${Gene_file##*/}
	
	


