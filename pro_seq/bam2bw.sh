#!/bin/bash

# Boyuan_Li

# Mon Sep 23 23:13:15 2019

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -b <bam_file> -o <output_dir> -c <chr_file> "
	echo ""
	echo " -b	file         	[required] "
	echo ""
	echo " -o	dir          	[required] "
	echo ""
	echo " -c	file         	[required] "
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "b:o:c:h" optionName
do
	case $optionName in
		b) bam_file="$OPTARG";;
		o) output_dir="$OPTARG";;
		c) chr_file="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $bam_file = "" ]]; then
	echo "the $bam_file file is needed "
	exit 1
elif [[ ! -f $bam_file ]]; then
	echo "$bam_file:   is not found"
	exit 2
fi

if [[ $output_dir = "" ]]; then
	echo "the $output_dir file is needed "
	exit 1
elif [[ ! -d $output_dir ]]; then
	 echo "$output_dir:   is not found"
	exit 2
fi

if [[ $chr_file = "" ]]; then
	echo "the $chr_file file is needed "
	exit 1
elif [[ ! -f $chr_file ]]; then
	echo "$chr_file:   is not found"
	exit 2
fi

TMPDIR=$output_dir
OUTPUT=$output_dir
f=$bam_file
CHINFO=$chr_file
echo " "
echo "Writing bigWigs:"

   j=`echo $f | awk -F"/" '{print $NF}' | cut -d \. -f 1 `
   echo $j > ${OUTPUT}/${j}.align.log

# in SE, MAP5 alwasys TRUE


   #if [[ "${RNA5}" == "R1_5prime" && "${OPP}" == "FALSE" ]] ; then ## report The 5 prime end of the RNA.   #like GRO-seq
   # if [[ "$SE_OUTPUT" == "G" ]] ; then
      bedtools bamtobed -i $f 2> ${TMPDIR}/kill.warnings| awk 'BEGIN{OFS="\t"} ($5 > 0){print $0}' | \
      awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,$6}; ($6 == "-") {print $1,$3-1,$3,$4,$5,$6}' | gzip > ${TMPDIR}/$j.bed.gz
   #elif [[ "${RNA3}" == "R1_5prime" && "${OPP}" == "TRUE" ]] ; then  #like PRO-seq
    # elif [[ "$SE_OUTPUT" == "P" ]] ; then
      # bedtools bamtobed -i $f 2> ${TMPDIR}/kill.warnings| awk 'BEGIN{OFS="\t"} ($5 > 0){print $0}' | \
      # awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,"-"}; ($6 == "-") {print $1,$3-1,$3,$4,$5,"+"}' | gzip > ${TMPDIR}/$j.bed.gz
   # fi


   echo 'Number of mappable reads:' >> ${OUTPUT}/${j}.align.log
   readCount=`zcat ${TMPDIR}/$j.bed.gz | grep "" -c`
   echo ${readCount} >> ${OUTPUT}/${j}.align.log
   
   ## Remove rRNA and reverse the strand (PRO-seq).
   zcat ${TMPDIR}/$j.bed.gz | grep "rRNA\|chrM" -v | grep "_" -v | sort-bed - | gzip > ${TMPDIR}/$j.nr.rs.bed.gz
   echo 'Number of mappable reads (excluding rRNA):' >> ${OUTPUT}/${j}.align.log
   echo `zcat ${TMPDIR}/$j.nr.rs.bed.gz | grep "" -c` >> ${OUTPUT}/${j}.align.log
   
   ## Convert to bedGraph ... Cannot gzip these, unfortunately.
   bedtools genomecov -bg -i ${TMPDIR}/$j.nr.rs.bed.gz -g ${CHINFO} -strand + > ${TMPDIR}/$j\_plus.bedGraph
   bedtools genomecov -bg -i ${TMPDIR}/$j.nr.rs.bed.gz -g ${CHINFO} -strand - > ${TMPDIR}/$j\_minus.noinv.bedGraph
   
   ## Invert minus strand.
   cat ${TMPDIR}/$j\_minus.noinv.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' > ${TMPDIR}/$j\_minus.bedGraph ## Invert read counts on the minus strand.
   
   ## normalized by RPM
   cat ${TMPDIR}/$j\_plus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*1000*1000/'$readCount'/1}' > ${TMPDIR}/$j\_plus.rpm.bedGraph  &
   cat ${TMPDIR}/$j\_minus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*1000*1000/'$readCount'/1}' > ${TMPDIR}/$j\_minus.rpm.bedGraph  &
   wait
   ## Then to bigWig (nomalized and non-nomrmalized ones)
   bedGraphToBigWig ${TMPDIR}/$j\_plus.rpm.bedGraph ${CHINFO} ${OUTPUT}/$j\_plus.rpm.bw 
   bedGraphToBigWig ${TMPDIR}/$j\_minus.rpm.bedGraph ${CHINFO} ${OUTPUT}/$j\_minus.rpm.bw 
   wait
   bedGraphToBigWig ${TMPDIR}/$j\_plus.bedGraph ${CHINFO} ${OUTPUT}/$j\_plus.bw 
   bedGraphToBigWig ${TMPDIR}/$j\_minus.bedGraph ${CHINFO} ${OUTPUT}/$j\_minus.bw 
#   rm ${TMPDIR}/$j.nr.rs.bed.gz ${TMPDIR}/$j.bed.gz ${TMPDIR}/$j*.bedGraph
