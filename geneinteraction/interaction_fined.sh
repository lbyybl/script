#!/bin/bash

# Boyuan-Li

# Thu Aug 16 20:41:37 CST 2018

# this to find the gene or gene and regulation elements interaction to draw interaction heatmap
# and don't distanguish the strand

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -g <Gene_bed_file> -p <Allavalidpairs_file> -o <Output_file> "
	echo ""
	echo " -g STRING          [required] Gene_bed_file can be file merged by promoter enhancer gff file."
	echo ""
	echo " -p STRING          [required] Allavalidpairs_file produced by hic-pro"
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
while getopts "g:p:o:h" optionName
do
	case $optionName in
		g) Gene_bed_file="$OPTARG";;
		p) Allavalidpairs_file="$OPTARG";;
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


if [[ $Allavalidpairs_file = "" ]]; then
	echo " -p the Allavalidpairs_file file is needed "
	exit 1
elif [[ ! -f $Allavalidpairs_file ]]; then
	echo "$Allavalidpairs_file:   is not found"
	exit 2
fi


if [[ $Output_file = "" ]]; then
	echo " -o the Output_file STRING is needed "
	exit 1
fi


#--- get abs path

Gene_file=$(readlink -e $Gene_bed_file)
Pairs_file=$(readlink -e $Allavalidpairs_file)
Output_file_name=${Output_file##*/}
Output_dir=${Output_file%%/*}
Prefix=${Output_file_name%.*}
mkdir -p ${Output_dir}
Output_dir=$(readlink -e ${Output_dir})
#--- get the pairs interaction file according to the input Allavalidpairs_file

#awk '{pos2=$3+50;pos4=$6+50;printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$2,$3,pos2,$4,$5,$6,pos4,$7,$1}' ${Pairs_file} | sort -k 1,1 -k 2,2n -k 3,3n >> _${Prefix}pairs.bed

#--- stastic the interaction intre gene and different elements
	#--- split the different strand
		#--- for gene.bed
		sort -k 1,1 -k 2,2n -k 3,3n ${Gene_file} -o ${Output_dir}/${Prefix}_${Gene_file##*/}
#		awk '{if ($4=="-") print $0}' ${Gene_file} >> _${Prefix}Ngene.bed
#		awk '{if ($4=="+") print $0}' ${Gene_file} >> _${Prefix}Pgene.bed
		#--- for pairs
		awk '{pos2=$3+150;pos4=$6+150;if ($4=="+" && $7=="+") printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$2,$3,pos2,$4,$5,$6,pos4,$7,$1}' ${Pairs_file} | sort -k 1,1 -k 2,2n -k 3,3n >> ${Output_dir}/_${Prefix}pairspp.bed &
		PS1=$!
		awk '{pos2=$3+150;pos4=$6-150;if ($4=="+" && $7=="-") printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$2,$3,pos2,$4,$5,pos4,$6,$7,$1}' ${Pairs_file} | sort -k 1,1 -k 2,2n -k 3,3n >> ${Output_dir}/_${Prefix}pairspn.bed &
		PS2=$!
		awk '{pos2=$3-150;pos4=$6-150;if ($4=="-" && $7=="-") printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$2,pos2,$3,$4,$5,pos4,$6,$7,$1}' ${Pairs_file} | sort -k 1,1 -k 2,2n -k 3,3n >> ${Output_dir}/_${Prefix}pairsnn.bed &
		PS3=$!
		awk '{pos2=$3-150;pos4=$6+150;if ($4=="-" && $7=="+") printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$2,pos2,$3,$4,$5,$6,pos4,$7,$1}' ${Pairs_file} | sort -k 1,1 -k 2,2n -k 3,3n >> ${Output_dir}/_${Prefix}pairsnp.bed &
		PS4=$!
		wait $PS1 $PS2 $PS3 $PS4
		
	#--- merge four files and sort
	#	cat ${Output_dir}/_${Prefix}pairspp.bed ${Output_dir}/_${Prefix}pairspn.bed ${Output_dir}/_${Prefix}pairsnn.bed ${Output_dir}/_${Prefix}pairsnp.bed | sort -k 1,1 -k 2,2n -k 3,3n >> ${Output_dir}/_${Prefix}pairsall.bed
	

	#--- got the reads corresponding gene
		bedtools intersect -a ${Output_dir}/_${Prefix}pairsnn.bed -b ${Output_dir}/${Prefix}_${Gene_file##*/} -sorted -wao | cut -f 1-9,14-15 | awk '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$5,$6,$7,$8,$1,$2,$3,$4,$9,$10,$11}' | sort -k 1,1 -k 2,2n -k 3,3n | bedtools intersect -a stdin -b ${Output_dir}/${Prefix}_${Gene_file##*/} -sorted -wao | cut -f 1-11,16-17 >> ${Output_dir}/_${Prefix}nnnfinal.bed &
		PS5=$!
		bedtools intersect -a ${Output_dir}/_${Prefix}pairspp.bed -b ${Output_dir}/${Prefix}_${Gene_file##*/} -sorted -wao | cut -f 1-9,14-15 | awk '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$5,$6,$7,$8,$1,$2,$3,$4,$9,$10,$11}' | sort -k 1,1 -k 2,2n -k 3,3n | bedtools intersect -a stdin -b ${Output_dir}/${Prefix}_${Gene_file##*/} -sorted -wao | cut -f 1-11,16-17 >> ${Output_dir}/_${Prefix}pppfinal.bed &
		PS6=$!
		bedtools intersect -a ${Output_dir}/_${Prefix}pairsnp.bed -b ${Output_dir}/${Prefix}_${Gene_file##*/} -sorted -wao | cut -f 1-9,14-15 | awk '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$5,$6,$7,$8,$1,$2,$3,$4,$9,$10,$11}' | sort -k 1,1 -k 2,2n -k 3,3n | bedtools intersect -a stdin -b ${Output_dir}/${Prefix}_${Gene_file##*/} -sorted -wao | cut -f 1-11,16-17 >> ${Output_dir}/_${Prefix}npfinal.bed &
		PS7=$!
		bedtools intersect -a ${Output_dir}/_${Prefix}pairspn.bed -b ${Output_dir}/${Prefix}_${Gene_file##*/} -sorted -wao | cut -f 1-9,14-15 | awk '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$5,$6,$7,$8,$1,$2,$3,$4,$9,$10,$11}' | sort -k 1,1 -k 2,2n -k 3,3n | bedtools intersect -a stdin -b ${Output_dir}/${Prefix}_${Gene_file##*/} -sorted -wao | cut -f 1-11,16-17 >> ${Output_dir}/_${Prefix}pnfinal.bed &
		PS8=$!
		wait $PS5 $PS6 $PS7 $PS8
		
		ls ${Output_dir}/_${Prefix}*final.bed | parallel gzip {}	
		ls ${Output_dir}/_${Prefix}pair*.bed | parallel gzip {}
	# #--- extract the interaction reads
	ls ${Output_dir}/_${Prefix}*final.bed.gz | parallel -k "awk '{if (\$11!=0 && \$13!=0) printf \"%s\\t%s\\t%s\\t%s\\n\",\$10,\$12,\$5,\$1}'"  '<(zcat {})'  | sort -k 1 -k 2 | uniq -c | awk '{printf "%s\t%s\t%s\t%s\t%s\n",$2,$3,$1,$4,$5}' | sort -k 3n >> ${Output_dir}/${Output_file_name}
	ls ${Output_dir}/_${Prefix}* | parallel mv {} ${Output_dir}/${Prefix}{/}
	rm -f 	${Output_dir}/${Prefix}*.gz

	
	# #--- got the reads corresponding gene
		# bedtools intersect -a ${Output_dir}/_${Prefix}pairsall.bed -b ${Output_dir}/${Gene_file##*/} -sorted -wao | cut -f 1-9,14-15 | awk '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$5,$6,$7,$8,$1,$2,$3,$4,$9,$10,$11}' | sort -k 1,1 -k 2,2n -k 3,3n | bedtools intersect -a stdin -b ${Output_dir}/${Gene_file##*/} -sorted -wao | cut -f 1-11,16-17 >> ${Output_dir}/_${Prefix}final.bed &
		# PS5=$!
		# wait $PS5 
		
		# cat ${Output_dir}/_${Prefix}final.bed | parallel --pipe --block 20M -k "gzip >> ${Output_dir}/_${Prefix}final.bed.gz"
		# rm ${Output_dir}/_${Prefix}final.bed
		# ls _${Prefix}pair*.bed | parallel gzip {}
	# # #--- extract the interaction reads
		# # ls _${Prefix}*final.bed | parallel -k "awk '{if (\$11!=0 && \$13!=0) print \$0}'"  {} | awk '{printf "%s\t%s\t%s\t%s\n",$10,$12,$5,$1}' >> _${Prefix}gene.interaction
	# # #--- stastic interaction count
		# # sort -k 1 -k 2 _${Prefix}gene.interaction | uniq -c | awk '{printf "%s\t%s\t%s\t%s\t%s\n",$2,$3,$1,$4,$5}' | sort -k 3n >> $Output_file
	# #--- modify above to one line
	# zcat ${Output_dir}/_${Prefix}final.bed.gz | parallel --pipe -k "awk '{if (\$11!=0 && \$13!=0) printf \"%s\\t%s\\t%s\\t%s\\n\",\$10,\$12,\$5,\$1}'"  | sort -k 1 -k 2 | uniq -c | awk '{printf "%s\t%s\t%s\t%s\t%s\n",$2,$3,$1,$4,$5}' | sort -k 3n >> ${Output_dir}/Output_file_name
	# ls ${Output_dir}/_${Prefix}* | parallel mv {} ${Output_dir}/${Prefix}{/}
	
