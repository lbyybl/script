#!/bin/bash

# Boyuan_Li

# Tue Feb 19 16:05:23 2019

# 写一个能够mapping，转成bigwig文件，还能去接头，call peak的脚本，这个脚本还应该有模块化的功能比如，去接头与质控，mapping与转bigwig，call peak就像hicpro一样；

# 1. 质控 2.去接头 3. mapping 4. 转bigwig 5.call peak 

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -f <fastq1> -p <Prefix> -o <outputdir> -b <blacklist> -m <min_length> -x <index> -n <normalize_method> -s <step>"
	echo ""
	echo " -f	file         	[required] fastq1"
	echo ""
	echo " -p	string       	[required] output prefix"
	echo ""
	echo " -o	dir          	[required] output dir"
	echo ""
	echo " -b	file         	[optional] blacklist difult:/DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed"
	echo ""
	echo " -m	number         	[optional] min_length, defult:25"
	echo ""
	echo " -x	path and Prefix [optional] index using for mapping defult:/DATA/software/bcbio/genomes/Mmusculus/mm10/bowtie2/mm10"
	echo ""
	echo " -n	string         	[optional] normalize_method, defult:RPKM;RPKM, CPM, BPM, RPGC, None"
	echo ""
	echo " -s	string         	[optional] which step you want to run: QC,cutadapt,mapping,bigwig,callpeak"
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}

min_length=25
index=/DATA/software/bcbio/genomes/Mmusculus/mm10/bowtie2/mm10
normalize_method=RPKM
blacklist=/DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed

if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "f:r:a:A:p:o:b:m:x:n:s:h" optionName
do
	case $optionName in
		f) fastq1="$OPTARG";;
		r) fastq2="$OPTARG";;
		a) forward_adapt="$OPTARG";;
		A) reverse_adapt="$OPTARG";;
		p) Prefix="$OPTARG";;
		o) outputdir="$OPTARG";;
		b) blacklist="$OPTARG";;
		m) min_length="$OPTARG";;
		x) index="$OPTARG";;
		n) normalize_method="$OPTARG";;
		s) step="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $fastq1 = "" ]]; then
	echo "the $fastq1 file is needed "
	exit 1
elif [[ ! -f $fastq1 ]]; then
	echo "$fastq1:   is not found"
	exit 2
fi

if [[ $Prefix = "" ]]; then
	echo " the $Prefix string is needed "
	exit 1
fi

if [[ $outputdir = "" ]]; then
	echo "the $outputdir file is needed "
	exit 1
elif [[ ! -d $outputdir ]]; then
	 echo "$outputdir:   is not found"
	exit 2
fi

if [[ $blacklist = "" ]]; then
	echo "the $blacklist file is needed "
	exit 1
elif [[ ! -f $blacklist ]]; then
	echo "$blacklist:   is not found"
	exit 2
fi

# 写一个能够mapping，转成bigwig文件，还能去接头，call peak的脚本，这个脚本还应该有模块化的功能比如，去接头与质控，mapping与转bigwig，call peak就像hicpro一样；

# 1. 质控 2.去接头 3. mapping 4. 转bigwig 5.call peak 

# fastq1
# fastq2
# forward_adapt
# reverse_adapt
# min_length # defult 25
# Prefix
# index # defult is /DATA/software/bcbio/genomes/Mmusculus/mm10/bowtie2/mm10
# normalize_method # defult is RPKM
# blacklist # defult is /DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed
# outputdir

#mkdir -p ${outputdir}/{bigwig,cutadapt,mapping,QC,rawdata}
mkdir -p ${outputdir}/QC
#--- 1. 质控
fastqc ${fastq1} -o ${outputdir}/QC

#--- 2. 去低质量 
mkdir -p ${outputdir}/cutadapt
cutadapt -q 15,15 -m ${min_length} -o ${outputdir}/cutadapt/${Prefix}_1.fq ${fastq1}
fastqc ${outputdir}/cutadapt/${Prefix}_1.fq -o ${outputdir}/cutadapt
#--- 3. mapping
mkdir -p ${outputdir}/mapping

bowtie2 -q -x ${index} --very-sensitive-local -L 30 --score-min G,20,8 --reorder --rg SM:${Prefix} -U ${outputdir}/cutadapt/${Prefix}_1.fq --rg-id ${Prefix} -p 10 -S ${outputdir}/mapping/${Prefix}.sam 2> ${outputdir}/mapping/${Prefix}_mapping.stat

#bowtie2 -q -x ${index} -X 2000 -U ${outputdir}/cutadapt/${Prefix}_1.fq --rg-id ${Prefix} --rg SM:${Prefix} --rg PL:illumina -p 10 -S ${outputdir}/mapping/${Prefix}.sam 2> ${outputdir}/mapping/${Prefix}_mapping.stat
samtools sort -@ 10 ${outputdir}/mapping/${Prefix}.sam >> ${outputdir}/mapping/${Prefix}.bam
#rm -f ${outputdir}/mapping/${Prefix}.sam
sambamba markdup -r -t 10 ${outputdir}/mapping/${Prefix}.bam ${outputdir}/mapping/${Prefix}_rmdup.bam



#--- 4. 转bigwig
mkdir -p ${outputdir}/bigwig
bedtools merge -i ${blacklist} >> ${outputdir}/bigwig/${Prefix}_blacklist
bamCoverage -b ${outputdir}/mapping/${Prefix}_rmdup.bam -o ${outputdir}/bigwig/${Prefix}_rmdup.bw --normalizeUsing ${normalize_method} --blackListFileName ${outputdir}/bigwig/${Prefix}_blacklist

#--- 5. call peak
mkdir -p ${outputdir}/peak

macs2 callpeak -t ${outputdir}/mapping/${Prefix}_rmdup.bam -f BAM -n ${Prefix} -B --nolambda -g mm --outdir ${outputdir}/peak --nomodel --shift 100 --extsize 200 -q 0.00001 -m 10 50

#--- stastics
#--- total reads
pwd
suffix=$(echo ${fastq1} | awk -F "." '{print $NF}')
if [[ $suffix == 'gz' ]]; then
	total_resds=$(grep "^@" <(zcat ${fastq1}) | wc -l)
elif [[ $suffix != 'gz' ]]; then
	total_resds=$(grep "^@" ${fastq1} | wc -l)
fi

cutadapt=$(grep "^@" ${outputdir}/cutadapt/${Prefix}_1.fq | wc -l)
mapping=$(samtools view -F 4 ${outputdir}/mapping/${Prefix}.bam | wc -l)
rmdup=$(samtools view -F 4 ${outputdir}/mapping/${Prefix}_rmdup.bam | wc -l)
#---- stastic file
echo -e "sample_name,total_resds,cutadapt,mapping reads,rm duplicate" >> ${outputdir}/${Prefix}_stastic.csv
echo -e "${Prefix},${total_resds},${cutadapt},${mapping},${rmdup}" >> ${outputdir}/${Prefix}_stastic.csv
sed -e "s/^[ \t]*//g" -e 's/ /,/g' ${outputdir}/mapping/${Prefix}_mapping.stat >> ${outputdir}/mapping/${Prefix}_mapping_format.stat
Rscript /home/boyuanli/bashscript/bin/chip_ATAC/stastic_reads_info.r -b ${outputdir}/mapping/${Prefix}_mapping_format.stat -s ${outputdir}/${Prefix}_stastic.csv -o ${outputdir}/${Prefix}_stastic_all.csv
rm -f ${outputdir}/mapping/${Prefix}.bam