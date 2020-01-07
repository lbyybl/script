#!/bin/bash

# Boyuan_Li

# Fri Nov 29 15:10:59 2019

# This script is used to evalute the chip-seq quality;

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -d <dir> -p <prefix>"
	echo ""
	echo " -d	dir          	[required] the dir inclucde the mapping,unique/{peak,bw}"
	echo ""
	echo " -p	string          [required] the prefix"
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "d:p:h" optionName
do
	case $optionName in
		d) dir="$OPTARG";;
		p) prefix="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $dir = "" ]]; then
	echo "the $dir file is needed "
	exit 1
elif [[ ! -d $dir ]]; then
	 echo "$dir:   is not found"
	exit 2
fi

if [[ $prefix = "" ]]; then
        echo " the $prefix string is needed "
        exit 1
fi

#--- merge peak; calculate coorlation; PCA

#cd ${dir}/unique/peak/
cd ${dir}
	# ls *.narrowPeak | xargs -n2 | parallel -C " " "idr --samples {1} {2} --input-file-type narrowPeak --output-file {=1 s/-.*//g=}_overlap_peak.narrowPeak"
	cat ${dir}/unique/peak/*.narrowPeak > ${dir}/${prefix}_union_merge.narrowPeak
	
	bedtools merge -i  <(sort -k 1,1 -k 2n,2 -k 3n,3 ${dir}/${prefix}_union_merge.narrowPeak) -c 4,5,6,7,8,9,10 -o distinct,max,distinct,max,max,max,mean | awk '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%i\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | sort -k 1,1 -k 2n,2 -k 3n,3 -o ${dir}/${prefix}_union_merge.narrowPeak
	bedtools intersect -a ${dir}/${prefix}_union_merge.narrowPeak -b /DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed -v | cut -f 1-3 | sort -o ${dir}/${prefix}_union_merge.narrowPeak
	#awk '{st=$2+$10;en=st+1;print $1"\t"st"\t"en}' ${prefix}_union_merge.narrowPeak > ${prefix}_union_merge.summit
	#--- calculate coorlation and PCA
	
	multiBigwigSummary BED-file -b ${dir}/unique/bw/*.bw --BED ${dir}/${prefix}_union_merge.narrowPeak --smartLabels -out ${dir}/${prefix}_summary_peak.npz --outRawCounts ${dir}/${prefix}_summary_peak.tab -p 10 -bl /DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed
	plotCorrelation -in ${dir}/${prefix}_summary_peak.npz --corMethod spearman --skipZeros --plotTitle "spearman Correlation" --whatToPlot scatterplot -o ${dir}/${prefix}_summary_spearman_scatter.pdf --outFileCorMatrix ${dir}/${prefix}_summary_spearman_scatter.tab &
	plotCorrelation -in ${dir}/${prefix}_summary_peak.npz --corMethod spearman --skipZeros --plotTitle "spearman Correlation" --whatToPlot heatmap --colorMap YlOrRd --plotNumbers -o ${dir}/${prefix}_summary_spearman_heatmap.pdf --removeOutliers & 
	
	plotPCA -in ${dir}/${prefix}_summary_peak.npz -o ${dir}/${prefix}_summary_peak_PCA.pdf -T "PCA" --outFileNameData ${dir}/${prefix}_summary_peak_PCA.txt

	multiBigwigSummary bins -b ${dir}/unique/bw/*.bw --binSize 10000 --smartLabels -out ${dir}/${prefix}_summary_bin.npz --outRawCounts ${dir}/${prefix}_summary_bin.tab -p 10 -bl /DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed
	plotCorrelation -in ${dir}/${prefix}_summary_bin.npz --corMethod spearman --skipZeros --plotTitle "spearman Correlation" --whatToPlot scatterplot -o ${dir}/${prefix}_summary_bin_spearman_scatter.pdf --outFileCorMatrix ${dir}/${prefix}_summary_bin_spearman_scatter.tab &
	plotCorrelation -in ${dir}/${prefix}_summary_bin.npz --corMethod spearman --skipZeros --plotTitle "spearman Correlation" --whatToPlot heatmap --colorMap YlOrRd --plotNumbers -o ${dir}/${prefix}_summary_bin_spearman_heatmap.pdf --removeOutliers & 
	
	plotPCA -in ${dir}/${prefix}_summary_bin.npz -o ${dir}/${prefix}_summary_bin_PCA.pdf -T "PCA" --outFileNameData ${dir}/${prefix}_summary_bin_PCA.txt
		
	#--- enrich sufficiention
	plotFingerprint --skipZeros -bl /DATA/work/lbyybl/genomes/mm10/mm10.blacklist.bed -b ${dir}/unique/*.bam --smartLabels -T "Fingerprints of different samples" --plotFile ${dir}/${prefix}_fingerprints.pdf --outRawCounts ${dir}/${prefix}_fingerprints.txt \
	--outQualityMetrics ${dir}/${prefix}_fingerprints.metrics -p 10 1> ${dir}/${prefix}_fingerprints.log 2>&1
	
	#--- sequence deepth
	plotCoverage -b ${dir}/unique/*.bam --plotFile ${dir}/${prefix}_coverage.pdf -n 1000000 --plotTitle "coverage" --outRawCounts ${dir}/${prefix}_coverage.tab -e 1> ${dir}/${prefix}_coverage.log 2>&1
	
	