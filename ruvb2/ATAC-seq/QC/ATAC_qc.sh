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


#--- 将bam 文件shift
if [[ ! -d shift ]]; then
	
	mkdir -p shift
	alignmentSieve -b unique/${prefix}_rmdup_uniqe.bam --ATACshift -p 10 -o shift/${prefix}_rmdup_uniqe_shift.bam

else
	if [[ ! -f shift/${prefix}_rmdup_uniqe_shift.bam ]]; then
		alignmentSieve -b unique/${prefix}_rmdup_uniqe.bam --ATACshift -p 10 -o shift/${prefix}_rmdup_uniqe_shift.bam
	fi
fi

#-- 去除chrM
if [[ ! -d rmchrM ]]; then
	
	mkdir -p rmchrM
	samtools view -h shift/${prefix}_rmdup_uniqe_shift.bam | grep -v 'chrM' | samtools view -b - > rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM.bam
	samtools sort -@5 rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM.bam > rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM_sort.bam
	samtools index rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM_sort.bam

else
	if [[ ! -f rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM_sort.bam ]]; then
		if [[ ! -f rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM.bam ]]; then
			samtools view -h shift/${prefix}_rmdup_uniqe_shift.bam | grep -v 'chrM' | samtools view -b - > rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM.bam
			samtools sort -@5 rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM.bam > rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM_sort.bam
			samtools index rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM_sort.bam
		else
			samtools sort -@5 rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM.bam > rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM_sort.bam
			samtools index rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM_sort.bam
		fi
	fi
fi

#--- make bigwig
if [[ ! -d rmchrM/bw ]]; then
	mkdir -p rmchrM/bw
	bamCoverage --bam rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM_sort.bam -o rmchrM/bw/${prefix}_rpgc.bw --extendReads --binSize 1 -p 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500

else
	if [[ ! -f rmchrM/bw/${prefix}_rpgc.bw ]]; then
		bamCoverage --bam rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM_sort.bam -o rmchrM/bw/${prefix}_rpgc.bw --extendReads --binSize 1 -p 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500
	fi
fi

#-- call peak
if [[ ! -d rmchrM/peak ]]; then
	mkdir -p rmchrM/peak
	if [[ ! -f rmchrM/peak/${prefix}_peaks.narrowPeak ]]; then		
		macs2 callpeak -t rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM_sort.bam -f BAMPE --call-summits -n ${prefix} -B -g mm --outdir rmchrM/peak
	fi
else
	if [[ ! -f rmchrM/peak/${prefix}_peaks.narrowPeak ]]; then
		macs2 callpeak -t rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM_sort.bam -f BAMPE --call-summits -n ${prefix} -B -g mm --outdir rmchrM/peak
	fi
fi

	#--- 为了后续统计PBC1，PBC2和NRF，先筛选出coordinate的reads
if [[ ! -d ATAC_QC ]]; then
	
	mkdir -p ATAC_QC
	samtools view -H $bam_file > ATAC_QC/${prefix}-cor.sam
	samtools view $bam_file | grep "YT:Z:CP" >> ATAC_QC/${prefix}-cor.sam
	samtools sort -n ATAC_QC/${prefix}-cor.sam > ATAC_QC/${prefix}-cor.bam
	rm -f ATAC_QC/${prefix}-cor.sam

else
	if [[ ! -f ATAC_QC/${prefix}-cor.bam ]]; then
		samtools view -H $bam_file > ATAC_QC/${prefix}-cor.sam
		samtools view $bam_file | grep "YT:Z:CP" >> ATAC_QC/${prefix}-cor.sam
		samtools sort -n ATAC_QC/${prefix}-cor.sam > ATAC_QC/${prefix}-cor.bam
		rm -f ATAC_QC/${prefix}-cor.sam
	fi
fi
	
#--- fragment length 
if [[ ! -f ATAC_QC/${prefix}_FL_count.txt ]]; then
	
	samtools view rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM_sort.bam | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > ATAC_QC/${prefix}_FL_count.txt

fi
	
if [[ ! -f $output_file_name ]]; then
	
	#--- chrM的reads
	file=unique/${prefix}_rmdup_uniqe.bam
	mtReads=$(samtools idxstats $file | grep 'chrM' | cut -f 3)
	totalReads=$(samtools idxstats $file | awk '{SUM += $3} END {print SUM}')
	chrM_ratio=$(echo $mtReads/$totalReads | bc -l)
	#--- PBC matrix
	if [[ ! -f ATAC_QC/${prefix}_pbc.qc.txt ]]; then
		bedtools bamtobed -bedpe -i ATAC_QC/${prefix}-cor.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0}($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n", mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ATAC_QC/${prefix}_pbc.qc.txt
		PBC1=$(cut -f 6 ATAC_QC/${prefix}_pbc.qc.txt)
		PBC2=$(cut -f 7 ATAC_QC/${prefix}_pbc.qc.txt)
		NRF=$(cut -f 5 ATAC_QC/${prefix}_pbc.qc.txt)
	else
		PBC1=$(cut -f 6 ATAC_QC/${prefix}_pbc.qc.txt)
		PBC2=$(cut -f 7 ATAC_QC/${prefix}_pbc.qc.txt)
		NRF=$(cut -f 5 ATAC_QC/${prefix}_pbc.qc.txt)	
	fi
	#结果为mt(总coordinate),m0(去除dup之后的reads),m1(没有dup的reads),m2(2次dup的reads), NRF,PBC1,PBC2,
	#--- FRIP
	if [[ ! -f ATAC_QC/${prefix}_featurecounts.txt ]]; then
		awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"} {print $4, $1, $2, $3, "."}' rmchrM/peak/${prefix}_peaks.narrowPeak > ATAC_QC/${prefix}.saf
		featureCounts -p -a ATAC_QC/${prefix}.saf -F SAF -o ATAC_QC/${prefix}_readCountInPeaks.txt rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM_sort.bam 2> ATAC_QC/${prefix}_featurecounts.txt
		total_reads=$(grep "Total fragments" ATAC_QC/${prefix}_featurecounts.txt | awk '{print $5}')
		reads_in_peaks=$(grep "Successfully assigned fragments" ATAC_QC/${prefix}_featurecounts.txt | awk '{print $6}')
		FRiP=$(awk "BEGIN {print "${reads_in_peaks}"/"${total_reads}"}")
	else 
		total_reads=$(grep "Total fragments" ATAC_QC/${prefix}_featurecounts.txt | awk '{print $5}')
		reads_in_peaks=$(grep "Successfully assigned fragments" ATAC_QC/${prefix}_featurecounts.txt | awk '{print $6}')
		FRiP=$(awk "BEGIN {print "${reads_in_peaks}"/"${total_reads}"}")
		
	fi
	#---TSS enrichment score
	if [[ ! -f ATAC_QC/${prefix}_tss.txt ]]; then
		awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"} ($6=="+"){print $4, $1, $2-1000, $2+1000, "+"} ($6=="-") {print $4, $1, $3-1000, $3+1000, "-"}' /DATA/work/lbyybl/genomes/mm10/DNA_elements_made_by_Boyuan/mm10_gene.bed > ATAC_QC/${prefix}_tss.saf
		featureCounts -p -a ATAC_QC/${prefix}_tss.saf -F SAF -o ATAC_QC/${prefix}_tss_readCountInPeaks.txt rmchrM/${prefix}_rmdup_uniqe_shift_rmchrM_sort.bam 2> ATAC_QC/${prefix}_tss.txt
		total_reads=$(grep "Total fragments" ATAC_QC/${prefix}_tss.txt | awk '{print $5}')
		reads_in_tss=$(grep "Successfully assigned fragments" ATAC_QC/${prefix}_tss.txt | awk '{print $6}')
		TSSE=$(awk "BEGIN {print "${reads_in_tss}"/"${total_reads}"}")
	else
		total_reads=$(grep "Total fragments" ATAC_QC/${prefix}_tss.txt | awk '{print $5}')
		reads_in_tss=$(grep "Successfully assigned fragments" ATAC_QC/${prefix}_tss.txt | awk '{print $6}')
		TSSE=$(awk "BEGIN {print "${reads_in_tss}"/"${total_reads}"}")	
	fi

	echo "sample_ID,chrM,PBC1,PBC2,NRF,FRiP,TSSE" > $output_file_name
	echo -e "${prefix},${chrM_ratio},${PBC1},${PBC2},${NRF},${FRiP},${TSSE}" >> $output_file_name
	

fi

if [[ -f ATAC_QC/${prefix}_tss.txt ]]; then
	echo '' > ATAC_QC/${prefix}-cor.bam
	echo '' > ATAC_QC/${prefix}_readCountInPeaks.txt 
	echo '' > ATAC_QC/${prefix}.saf 
	echo '' > ATAC_QC/${prefix}_tss_readCountInPeaks.txt 
	echo '' > ATAC_QC/${prefix}_tss.saf
	echo '' > shift/${prefix}_rmdup_uniqe_shift.bam
	
fi