#!/bin/bash

# Boyuan-Li

# Wed Dec 12 22:15:14 CST 2018

# This script is used to ypj's 4c data local mapping
# 1. select enzyme site 2. select 3' primer 4. cut 5' adaptor 5. local mapping 6. rm duplicate 7. rm religation and self ligation 8. rm dumped 9. get bigwig 10. call peak
# 11. stastic
set -eo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -f <Reads1_file> -r <Reads2_file> -o <Output_dir> -m <Primer3_file> -p <Prefix> "
	echo ""
	echo " -f STRING          [required] forward fastq file"
	echo ""
	echo " -r STRING          [required] reverse fastq file"
	echo ""
	echo " -o STRING          [required] Output dir"
	echo ""
	echo " -m STRING          [required] forward reads1 5' primer file used to select reads"
	echo ""
	echo " -p STRING          [required] output file Prefix "
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "f:r:o:m:p:h" optionName
do
	case $optionName in
		f) Reads1_file="$OPTARG";;
		r) Reads2_file="$OPTARG";;
		o) Output_dir="$OPTARG";;
		m) Primer3_file="$OPTARG";;
		p) Prefix="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Reads1_file = "" ]]; then
	echo " -f the Reads1_file file is needed "
	exit 1
elif [[ ! -f $Reads1_file ]]; then
	echo "$Reads1_file:   is not found"
	exit 2
fi


if [[ $Reads2_file = "" ]]; then
	echo " -r the Reads2_file file is needed "
	exit 1
elif [[ ! -f $Reads2_file ]]; then
	echo "$Reads2_file:   is not found"
	exit 2
fi

if [[ $Output_dir = "" ]]; then
	echo " -o the Output_dir directory is needed "
	exit 1
elif [[ ! -d $Output_dir ]]; then
	echo "$Output_dir:   is not found"
	exit 2
fi


if [[ $Primer3_file = "" ]]; then
	echo " -m the Primer3_file file is needed "
	exit 1
elif [[ ! -f $Primer3_file ]]; then
	echo "$Primer3_file:   is not found"
	exit 2
fi


if [[ $Prefix = "" ]]; then
	echo " -p the Prefix STRING is needed "
	exit 1
fi

linkerA=CATG
linkerB=CATG
primer3_file=$(readlink -e ${Primer3_file})
#primer5_file=/gpfs1/xiongji_pkuhpc/lbyybl/collaberate/YB/rawdata/DE19-DOX-2293_TKR180700929_HNM5HCCXY_L6/cut_adaptor/split_barcode/test/AD38A2.5primer
Output_dir=$(readlink -e $Output_dir)
reads1=$(readlink -e ${Reads1_file})
reads2=$(readlink -e ${Reads2_file})
prefix=${Prefix}
blacklist_file=/DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed
bowtie2_index=/DATA/genomes/mm10/mm10
fragment_bed_file=/home/boyuanli/tools/hic-pro/HiC-Pro_2.9.0/annotation/mm10_nlaIII.bed
#self_re_lig=/DATA/work/lbyybl/ypjiang/4c/4c20181216/self_re.bed
self_re_lig=$(readlink -e ${selfrelig})
mkdir -p ${Output_dir}
ln -s ${reads1} ${Output_dir}/read1.fq
ln -s ${reads2} ${Output_dir}/read2.fq
cd ${Output_dir}
mkdir -p trimlinker select_primter mapping

#--- select reads with enzyme site no enzyme site, so use the link to avoid more modify;
	# trimLinker -t 12 -m 2 -k 1 -l 16 -o trimlinker/ -n ${prefix} -A ${linkerA} -B ${linkerB} read1.fq read2.fq
	# mv trimlinker/${prefix}_1.valid.fastq trimlinker/${prefix}.1.fastq
	# mv trimlinker/${prefix}_2.valid.fastq trimlinker/${prefix}.2.fastq
	ln -s ${reads1} trimlinker/${prefix}.1.fastq
	ln -s ${reads2} trimlinker/${prefix}.2.fastq
#--- select 3'primer
	fastq-multx -B ${primer3_file} -b trimlinker/${prefix}.1.fastq trimlinker/${prefix}.2.fastq -o select_primter/%.1.fastq  -o select_primter/%.2.fastq -m 3 -x

#--- mapping 	
	# bowtie2 --very-sensitive-local -L 30 --score-min G,20,8 --reorder --rg-id BML --rg SM:${prefix} -p 30 -x ${bowtie2_index}  -U select_primter/${prefix}.2.fastq -S mapping/${prefix}_local.sam 2> mapping/${prefix}_local.log
	bowtie2 -q -x ${bowtie2_index} -X 2000 -1 select_primter/${prefix}.1.fastq -2 select_primter/${prefix}.2.fastq --rg-id ${prefix} --rg SM:${prefix} --rg PL:illumina -p 3 -S mapping/${prefix}_local.sam 2> mapping/${prefix}_local.log
#---rm self re, duplicate, dumped
	pwd
	cd mapping/
	pwd
	mkdir -p rm_dup unique # rm_self_re rm_dumped bigwig
	#--- rm duplicate		
		samtools sort -@ 10 ${prefix}_local.sam >> ${prefix}_local.bam
		sambamba markdup -r -t 10 ${prefix}_local.bam rm_dup/${prefix}_local.rmdup.bam
	#--- rm religation and self ligation 
		# bedtools intersect -a rm_dup/${prefix}_local.rmdup.bam -b ${self_re_lig} -wa -v >> rm_self_re/${prefix}_local.rmselfre.bam
	# #--- rm dumped
		# awk '{st=$2-900;en=$2+900;print $1,st,en}' ${fragment_bed_file} >> dump.bed
		# awk '{if ($2 < 0) $2=0; printf "%s\t%s\t%s\n",$1,$2,$3}' dump.bed >> dump.txt
		# sed -i '/HBV/ d' dump.txt
		# bedtools intersect -a rm_self_re/${prefix}_local.rmselfre.bam -b dump.txt -wa >> rm_dumped/${prefix}_local.rmdump.bam
		# samtools index rm_dumped/${prefix}_local.rmdump.bam
	#--- rm multiple
		samtools view -hF4 rm_dup/${prefix}_local.rmdup.bam | grep -v "XS:" | samtools view -b - >> unique/${prefix}_uniqe.bam
		samtools index unique/${prefix}_uniqe.bam
		# /home/boyuanli/bashscript/bin/pro_seq/unique_and_bw.sh -b rm_dup/${prefix}_local.rmdup.bam -o unique -p ${Prefix}_uniqe -B /DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed
	#--- get bigwig
		bedtools merge -i ${blacklist_file} >> bigwig/blacklist
		bamCoverage -b unique/${prefix}_uniqe.bam -o bigwig/${prefix}_rpkm.bw --extendReads --binSize 10 --normalizeUsing RPKM --blackListFileName bigwig/blacklist
		bamCoverage --bam unique/${prefix}_uniqe.bam -o bigwig/${prefix}_rpgc.bw --extendReads --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --blackListFileName bigwig/blacklist
#--- call peak
	mkdir -p peak
	macs2 callpeak -t rm_dumped/${prefix}_local.rmdump.bam -f BAM -g 2.7e9 --nomodel --outdir peak/ -n ${prefix} --broad -q 0.05
#--- stastics
	#--- total reads
	pwd
	total_resds=$(grep "^@" <(zcat ../read1.fq) | wc -l)
	# sel_enz=$(grep "^@" ../trimlinker/${prefix}.1.fastq | wc -l)
	sel_3p=$(grep "^@" ../select_primter/${prefix}.1.fastq | wc -l)
	#sel_5p=$(grep "^@" ../select_5_primer/${prefix}.2.fastq | wc -l)
	mapping=$(samtools view -F 4 ${prefix}_local.sam | wc -l)
	rmdup=$(samtools view -F 4 rm_dup/${prefix}_local.rmdup.bam | wc -l)
	unique=$(samtools view -F 4 unique/${prefix}_uniqe.bam | wc -l)
	# rm_re_self=$(samtools view -F 4 rm_self_re/${prefix}_local.rmselfre.bam | wc -l)
	# rmdump=$(samtools view -F 4 rm_dumped/${prefix}_local.rmdump.bam | wc -l)
	#---- stastic file
	echo -e "sample_name,total_resds,select 3' primer, mapping reads,rm duplicate,unique" >> stastic.csv
	echo -e "${prefix},${total_resds},${sel_3p},${mapping},${rmdup},${unique}" >> stastic.csv
