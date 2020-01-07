#!/bin/bash

# Boyuan-Li

# Thu Nov 15 21:47:15 CST 2018

# this to find the gene or gene and regulation elements interaction to draw interaction heatmap
# and don't distanguish the strand 
# and used the loop file called by hichipper

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -g <Gene_bed_file> -l <Loop_file> -o <Output_file> "
	echo ""
	echo " -g STRING          [required] Gene_bed_file can be file merged by promoter enhancer gff file."
	echo ""
	echo " -l STRING          [required] Loop files called by hichipper"
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
while getopts "g:l:o:h" optionName
do
	case $optionName in
		g) Gene_bed_file="$OPTARG";;
		l) Loop_file="$OPTARG";;
		o) Output_file="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Gene_bed_file = "" ]]; then
	echo " -g the Gene_bed_file file is needed "
	exit 1
elif [[ ! -f $Gene_bed_file ]]; then
	echo "$Gene_bed_file:   is not found"
	exit 2
fi


if [[ $Loop_file = "" ]]; then
	echo " -l the Loop_file file is needed "
	exit 1
elif [[ ! -f $Loop_file ]]; then
	echo "$Loop_file:   is not found"
	exit 2
fi


if [[ $Output_file = "" ]]; then
	echo " -o the Output_file STRING is needed "
	exit 1
fi

#--- get abs path

		Gene_file=$(readlink -e $Gene_bed_file)
		looppe_file=$(readlink -e $Loop_file)
		Output_file_name=${Output_file##*/}
		Output_dir=${Output_file%%/*}
		Prefix=${Output_file_name%.*}
        dir=$(echo -e "${Output_file}" | grep "/")
        if [[ $dir != "" && $dir != "./" ]]; then
                        Output_dir=${Output_file%/*}
                        if [[ ! -d $Output_dir ]]; then
                                mkdir -p ${Output_dir}
                        fi
        fi

		Output_dir=$(readlink -e ${Output_dir})
#--- get the pairs interaction file according to the input Allavalidpairs_file

#--- stastic the interaction intre gene and different elements
	#--- split the different strand
		#--- for gene.bed
		sort -k 1,1 -k 2,2n -k 3,3n ${Gene_file} -o ${Output_dir}/${Prefix}_${Gene_file##*/}

	#--- got the reads corresponding gene
		bedtools intersect -a <(awk '{OFS="\t";print $1,$2,$3,$4,$5,$6}' ${looppe_file}) -b ${Output_dir}/${Prefix}_${Gene_file##*/} -wao | awk '{OFS="\t";print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}' | sort -k 1,1 -k 2,2n -k 3,3n | bedtools intersect -a stdin -b ${Output_dir}/${Prefix}_${Gene_file##*/} -sorted -wao | awk '{OFS="\t";print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' >> ${Output_dir}/_${Prefix}final.bed &
		#bedtools intersect -a <(awk '{OFS="\t";print $1,$2,$3,$4,$5,$6}' ${looppe_file}) -b ${Output_dir}/${Prefix}_${Gene_file##*/} -wao | sort -k 1,1 -k 2,2n -k 3,3n | bedtools intersect -a stdin -b ${Output_dir}/${Prefix}_${Gene_file##*/} -sorted -wao | cut -f 1-11,16-17 >> ${Output_dir}/_${Prefix}nnnfinal.bed &
		PS5=$!
		
		wait $PS5
		
		ls ${Output_dir}/_${Prefix}final.bed | parallel gzip {}	
		
	# #--- extract the interaction reads
	ls ${Output_dir}/_${Prefix}final.bed.gz | parallel -k "awk '{if (\$12!=0 && \$18!=0) printf \"%s\\t%s\\t%s\\t%s\\n\",\$11,\$17,\$7,\$13}'"  '<(zcat {})'  | sort -k 1 -k 2 | uniq -c | awk '{printf "%s\t%s\t%s\t%s\t%s\n",$2,$3,$1,$4,$5}' | sort -k 3n >> ${Output_dir}/${Output_file_name}
	ls ${Output_dir}/_${Prefix}* | parallel mv {} ${Output_dir}/${Prefix}{/}
		