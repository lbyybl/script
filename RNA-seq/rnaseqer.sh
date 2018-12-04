#!/bin/bash

# Boyuan-Li

# Tue Dec  4 15:32:47 CST 2018

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -i <Input_dir> -d <Index> -g <GFF_file> -o <Output_dir>-p <CUP_num> "
	echo ""
	echo " -i STRING          [required] Input dir contain all sample"
	echo ""
	echo " -d STRING          [required] Index path and the prefix of index"
	echo ""
	echo " -g STRING          [required] GFF/GTF file"
	echo ""
	echo " -o STRING          [required] Output dir"
	echo ""
	echo " -p STRING          [optional] CPU num of per sample, defuld is 1"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}

CUP_num=1
if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "i:d:g:o:p:h" optionName
do
	case $optionName in
		i) Input_dir="$OPTARG";;
		d) Index="$OPTARG";;
		g) GFF_file="$OPTARG";;
		o) Output_dir="$OPTARG";;
		p) CUP_num="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Input_dir = "" ]]; then
	echo " -i the Input_dir directory is needed "
	exit 1
elif [[ ! -d $Input_dir ]]; then
	echo "$Input_dir:   is not found"
	exit 2
fi


if [[ $Index = "" ]]; then
	echo " -d the Index STRING is needed "
	exit 1
fi


if [[ $GFF_file = "" ]]; then
	echo " -g the GFF_file file is needed "
	exit 1
elif [[ ! -f $GFF_file ]]; then
	echo "$GFF_file:   is not found"
	exit 2
fi

if [[ $Output_dir = "" ]]; then
	echo " -i the Output_dir directory is needed "
	exit 1
elif [[ ! -d $Output_dir ]]; then
	echo "$Output_dir:   is not found"
	echo "mkdir $Output_dir"
	mkdir -p $Output_dir
fi

if [[ $CUP_num = "" ]]; then
	echo " -p the CUP_num STRING is needed "
	exit 1
fi

index_dir=$(readlink -e ${Index%%/*})
index=${index_dir}/${Index##*/}
gff_file=$(readlink -e $GFF_file)
outdir=$(readlink -e $Output_dir)
echo "enter input dir....."

	cd $Input_dir
	sample_num=$(ls | wc -l)
	cup2=$(echo $CUP_num*sample_num/2 | bc -l)
	mkdir -p ${outdir}/mapping
	ls | sed -e 's/_1.*$//g' -e 's/_2.*$//g' -e 's/_R1.*$//g' -e 's/_R2.*$//g' | sort | uniq >> ${outdir}/mapping/namelist.txt
	ls | xargs -n 2 | parallel -C " " -j ${cup2} "hisat2 -x ${index} -1 {1} -2 {2} -S ${outdir}/mapping/{=1 s/_1\..*$//g;s/_R1\..*$//g;s/_1_.*$//g;s/_1_.*$//g;s/_R1.*$//g=}.sam --dta --dta-cufflinks"
	cd ${outdir}/mapping 
	ls *.sam | parallel -j ${cup2} "samtools view -bS {} | samtools sort - -o {.}_sort.bam; samtools index {.}_sort.bam;rm {}"
	mkdir -p ${outdir}/stringtie 
	ls *.bam | parallel -j ${cup2} "stringtie -p $CUP_num -G ${gff_file} -o ${outdir}/stringtie/{.}.gtf -l {.} {}"
	cd ${outdir}/stringtie
	ls *.gtf >> mergelist.txt
	stringtie --merge -p ${cup2} -G ${gff_file} -o stringtie_merged.gtf mergelist.txt
	gffcompare -r ${gff_file} -G -o merge ./stringtie_merged.gtf
	mkdir -p ${outdir}/ballgown
	cat ${outdir}/mapping/namelist.txt | parallel "stringtie -e -B -p $CUP_num -G ${outdir}/stringtie/stringtie_merged.gtf -o ${outdir}/ballgown/{}/{}.gtf ${outdir}/mapping/{}_sort.bam"
	
