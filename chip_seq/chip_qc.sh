#!/bin/bash

# Boyuan_Li

# Fri Nov 22 16:39:53 2019

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -b <bam_file> -n <output_file_name> "
	echo ""
	echo " -b	file         	[required] The bam file that with duplicate"
	echo ""
	echo " -n	string       	[required] output file name"
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "b:n:h" optionName
do
	case $optionName in
		b) bam_file="$OPTARG";;
		n) output_file_name="$OPTARG";;
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

if [[ $output_file_name = "" ]]; then
	echo " the $output_file_name string is needed "
	exit 1
fi
filename=$(/usr/bin/basename $bam_file)
prefix=${filename%%.*}

#--- make bigwig
if [[ ! -d unique/bw ]]; then
	mkdir -p unique/bw
	bamCoverage --bam unique/${prefix}_rmdup_uniqe.bam -o unique/bw/${prefix}_rmdup_uniqe_rpgc.bw --extendReads --binSize 1 -p 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500

else
	if [[ ! -f unique/bw/${prefix}_rpgc.bw ]]; then
		bamCoverage --bam unique/${prefix}_rmdup_uniqe.bam -o unique/bw/${prefix}_rmdup_uniqe_rpgc.bw --extendReads --binSize 1 -p 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500
	fi
fi

#-- call peak
if [[ ! -d unique/peak ]]; then
	mkdir -p unique/peak
	if [[ ! -f unique/peak/${prefix}_peaks.narrowPeak ]]; then		
		macs2 callpeak -t unique/${prefix}_rmdup_uniqe.bam -f BAMPE --call-summits -n ${prefix} -B -g mm --outdir unique/peak
	fi
else
	if [[ ! -f unique/peak/${prefix}_peaks.narrowPeak ]]; then
		macs2 callpeak -t unique/${prefix}_rmdup_uniqe.bam -f BAMPE --call-summits -n ${prefix} -B -g mm --outdir unique/peak
	fi
fi

	#--- 为了后续统计PBC1，PBC2和NRF，先筛选出coordinate的reads
if [[ ! -d CHIP_QC ]]; then
	
	mkdir -p CHIP_QC
	samtools view -H $bam_file > CHIP_QC/${prefix}-cor.sam
	samtools view $bam_file | grep "YT:Z:CP" >> CHIP_QC/${prefix}-cor.sam
	samtools sort -n CHIP_QC/${prefix}-cor.sam > CHIP_QC/${prefix}-cor.bam
	rm -f CHIP_QC/${prefix}-cor.sam

else
	if [[ ! -f CHIP_QC/${prefix}-cor.bam ]]; then
		samtools view -H $bam_file > CHIP_QC/${prefix}-cor.sam
		samtools view $bam_file | grep "YT:Z:CP" >> CHIP_QC/${prefix}-cor.sam
		samtools sort -n CHIP_QC/${prefix}-cor.sam > CHIP_QC/${prefix}-cor.bam
		rm -f CHIP_QC/${prefix}-cor.sam
	fi
fi


	
#--- fragment length 
if [[ ! -f CHIP_QC/${prefix}_fragmentsize.pdf ]]; then
	
	# samtools view unique/${prefix}_rmdup_uniqe.bam | awk '{cnt++;if ($9>0) {sum+=$9;print $9} else {sum+=(-$9);print -$9}} END{print sum/cnt}' > CHIP_QC/${prefix}_FL_count.txt
	bamPEFragmentSize -b unique/${prefix}_rmdup_uniqe.bam -hist CHIP_QC/${prefix}_fragmentsize.pdf -T "Fragment size" --samplesLabel ${prefix} --table CHIP_QC/${prefix}_fragmentsize.txt -p 10
	FL=$(awk 'NR==2{print $5}' CHIP_QC/${prefix}_fragmentsize.txt)
	RL=$(awk 'NR==2{print $23}' CHIP_QC/${prefix}_fragmentsize.txt)
else 
	FL=$(awk 'NR==2{print $5}' CHIP_QC/${prefix}_fragmentsize.txt)
	RL=$(awk 'NR==2{print $23}' CHIP_QC/${prefix}_fragmentsize.txt)
fi
	
if [[ ! -f $output_file_name ]]; then
	
	#--- PBC matrix
	if [[ ! -f CHIP_QC/${prefix}_pbc.qc.txt ]]; then
		bedtools bamtobed -bedpe -i CHIP_QC/${prefix}-cor.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0}($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n", mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > CHIP_QC/${prefix}_pbc.qc.txt
		PBC1=$(cut -f 6 CHIP_QC/${prefix}_pbc.qc.txt)
		PBC2=$(cut -f 7 CHIP_QC/${prefix}_pbc.qc.txt)
		NRF=$(cut -f 5 CHIP_QC/${prefix}_pbc.qc.txt)
	else
		PBC1=$(cut -f 6 CHIP_QC/${prefix}_pbc.qc.txt)
		PBC2=$(cut -f 7 CHIP_QC/${prefix}_pbc.qc.txt)
		NRF=$(cut -f 5 CHIP_QC/${prefix}_pbc.qc.txt)	
	fi
	#结果为mt(总coordinate),m0(去除dup之后的reads),m1(没有dup的reads),m2(2次dup的reads), NRF,PBC1,PBC2,
	#--- FRIP
	if [[ ! -f CHIP_QC/${prefix}_enrichment.pdf ]]; then
		# awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"} {print $4, $1, $2, $3, "."}' unique/peak/${prefix}_peaks.narrowPeak > CHIP_QC/${prefix}.saf
		# featureCounts -p -a CHIP_QC/${prefix}.saf -F SAF -o CHIP_QC/${prefix}_readCountInPeaks.txt unique/${prefix}_rmdup_uniqe.bam 2> CHIP_QC/${prefix}_featurecounts.txt
		# total_reads=$(grep "Total fragments" CHIP_QC/${prefix}_featurecounts.txt | awk '{print $5}')
		# reads_in_peaks=$(grep "Successfully assigned fragments" CHIP_QC/${prefix}_featurecounts.txt | awk '{print $6}')
		# FRiP=$(awk "BEGIN {print "${reads_in_peaks}"/"${total_reads}"}")
		awk 'BEGIN{FS=OFS="\t"} ($6=="+"){print $1, $2-1000, $2+1000} ($6=="-") {print $1, $3-1000, $3+1000}' /DATA/work/lbyybl/genomes/mm10/DNA_elements_made_by_Boyuan/mm10_gene.bed > CHIP_QC/${prefix}_tss.bed
		bedtools merge -i <(sort -k 1,1 -k 2n,2 -k 3n,3 CHIP_QC/${prefix}_tss.bed) > CHIP_QC/${prefix}_tss2.bed
		
		plotEnrichment -b unique/${prefix}_rmdup_uniqe.bam --BED unique/peak/${prefix}_peaks.narrowPeak CHIP_QC/${prefix}_tss2.bed /DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed --regionLabels "peak" "TSS" "blacklist" \
		-o CHIP_QC/${prefix}_enrichment.pdf --outRawCounts CHIP_QC/${prefix}_enrichment.txt -p 10 -e --perSample --labels ${prefix} --variableScales
		FRIP=$(awk 'NR==3{print $3}' CHIP_QC/${prefix}_enrichment.txt)
		TSSE=$(awk 'NR==4{print $3}' CHIP_QC/${prefix}_enrichment.txt)
		Blacklist=$(awk 'NR==2{print $3}' CHIP_QC/${prefix}_enrichment.txt)
		
	else 
		FRIP=$(awk 'NR==3{print $3}' CHIP_QC/${prefix}_enrichment.txt)
		TSSE=$(awk 'NR==4{print $3}' CHIP_QC/${prefix}_enrichment.txt)
		Blacklist=$(awk 'NR==2{print $3}' CHIP_QC/${prefix}_enrichment.txt)
		
	fi
	# #---TSS enrichment score
	if [[ ! -f CHIP_QC/${prefix}_GCfreq.pdf ]]; then
		computeGCBias -b unique/${prefix}_rmdup_uniqe.bam --effectiveGenomeSize 2652783500 -g /DATA/work/lbyybl/genomes/mm10/mm10.2bit --GCbiasFrequenciesFile CHIP_QC/${prefix}_GCfreq.txt -p 10 --biasPlot CHIP_QC/${prefix}_GCfreq.pdf 1> CHIP_QC/${prefix}_GCfreq.log 2>&1

	fi

	echo "sample_ID,frag_leng,read_length,PBC1,PBC2,NRF,FRiP,TSSE,blacklist" > $output_file_name
	echo -e "${prefix},${FL},${RL},${PBC1},${PBC2},${NRF},${FRIP},${TSSE},${Blacklist}" >> $output_file_name
	

fi

if [[ -f CHIP_QC/${prefix}_tss.txt ]]; then
	echo '' > CHIP_QC/${prefix}-cor.bam
	echo '' > CHIP_QC/${prefix}_readCountInPeaks.txt 
	echo '' > CHIP_QC/${prefix}.saf 
	echo '' > CHIP_QC/${prefix}_tss_readCountInPeaks.txt 
	echo '' > CHIP_QC/${prefix}_tss.saf
	
fi