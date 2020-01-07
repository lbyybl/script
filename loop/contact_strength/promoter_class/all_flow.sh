#!/bin/bash

# Boyuan_Li

# Sun Jan 27 20:22:11 2019

# it's used to calculate four class promoter; all the bed file should the standard format;
# but you need debug first

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -g <gene_bed> -m <mapping_bed> -n <non_gene_bed> -d <output_dir> -p <prefix> "
	echo ""
	echo " -g	file         	[required] the gene bed file download form ucsc; and required bed fomat"
	echo ""
	echo " -m	file         	[required] the file convert from bam file; required bed format "
	echo ""
	echo " -n	file         	[required] the non gene bed file; chr st en . id"
	echo ""
	echo " -d	dir          	[required] the output dir"
	echo ""
	echo " -p	string       	[required] the output file prefix"
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "g:m:n:d:p:h" optionName
do
	case $optionName in
		g) gene_bed="$OPTARG";;
		m) mapping_bed="$OPTARG";;
		n) non_gene_bed="$OPTARG";;
		d) output_dir="$OPTARG";;
		p) prefix="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $gene_bed = "" ]]; then
	echo "the $gene_bed file is needed "
	exit 1
elif [[ ! -f $gene_bed ]]; then
	echo "$gene_bed:   is not found"
	exit 2
fi

if [[ $mapping_bed = "" ]]; then
	echo "the $mapping_bed file is needed "
	exit 1
elif [[ ! -f $mapping_bed ]]; then
	echo "$mapping_bed:   is not found"
	exit 2
fi

if [[ $non_gene_bed = "" ]]; then
	echo "the $non_gene_bed file is needed "
	exit 1
elif [[ ! -f $non_gene_bed ]]; then
	echo "$non_gene_bed:   is not found"
	exit 2
fi

if [[ $output_dir = "" ]]; then
	echo "the $output_dir file is needed "
	exit 1
elif [[ ! -d $output_dir ]]; then
	 echo "$output_dir:   is not found"
	exit 2
fi

if [[ $prefix = "" ]]; then
	echo " the $prefix string is needed "
	exit 1
fi

gene_bed=$(readlink -e ${gene_bed})
mapping_bed=$(readlink -e ${mapping_bed})
non_gene_bed=$(readlink -e ${non_gene_bed})
#--- 1. 根据gene bed file得到unique 的；2. 根据得到的unique的，得到缩短的延长的gene body和promoter,以及non gene bed中大于500k的区段 2.1 去除non gene bed中为N的区段；
# 3. 将mapping的bed文件unique一下；4. 对两类gene body和non gene算reads；5. 对promoter算mapping上的reads；6. 对prmote reads bed算最大的50bp；
# 7. 使用R计算四类gene 8. 根据gene 得到promoter
bash_dir=/home/boyuanli/bashscript/bin/loop/contact_strength/promoter_class
cd ${output_dir}
#--- 1. 这个脚本会根据输入的六列的标准的bed格式分正负链将gene uniqe；unique时取第一个，长的和短的共起点时使用短的（这一步真的有必要吗？因为看到是4类promoter，
# 故我怕短的是active，长的是unactive，最后看4类promoter的时候同一个promoter会在不同的类别中出现）；
awk '{ch=$3-$2; if(ch > 3000) print $0}' ${gene_bed} >> ${prefix}_gene.bed
Rscript /home/boyuanli/bashscript/bin/loop/contact_strength/promoter_class/unique_gene.r -i ${prefix}_gene.bed -o unique_short_gene.bed -d ./

#--- 2. 得到延长的gene body和缩短的gene body 和promoter的bed文件
${bash_dir}/bash/gotgenebodyformucsc.sh -g unique_short_gene.bed -r -1000 -o ./
${bash_dir}/bash/gotgenebodyformucsc.sh -g unique_short_gene.bed -r 1000 -o ./
${bash_dir}/bash/gotpromoterformucsc.sh -g unique_short_gene.bed -r 1000 -o ./

awk '$3-$2 > 500000 {print $0}' ${non_gene_bed} >> non_gene_region.bed
#---2.1 去除non gene bed中为N的区段
Rscript ${bash_dir}/clean_nongene.r -b non_gene_region.bed -o non_gene_region.bed -d ./
awk 'OFS="\t"{print $1,$2,$3,$5,$4,"."}' non_gene_region.bed | sort -k 1,1 -k 2n,2 -k 3n,3 -o non_gene_region.bed
#--- 3 将mapping的bed文件unique一下
sort -k 1,1 -k 2n,2 -k 3n,3 ${mapping_bed} | uniq >> gro_seq_mapping.bed

#--- 4 对两类gene body和non gene算reads；
${bash_dir}/bash/got_reads_num.sh -g genebody_final_-1000.bed -m gro_seq_mapping.bed -p genebody_final_-1000_reads_num ./
${bash_dir}/bash/got_reads_num.sh -g genebody_final_1000.bed -m gro_seq_mapping.bed -p genebody_final_1000_reads_num ./
bedtools intersect -a non_gene_region.bed -b gro_seq_mapping.bed -wao | awk '$7!=".",OFS="\t" {print $1,$2,$3,$4,$5,$6}' | sort -k 1,1 -k 2n,2 -k 3n,3 | uniq -c | awk 'OFS="\t" {print $2,$3,$4,$5,$6,$7,$1}' >> non_gene.bed
#--- 5 对promoter算mapping上的reads；
${bash_dir}/bash/got_reads.sh -g promoter_final_1000.bed -m gro_seq_mapping.bed -p promoter_final_1000_reads_num ./

#--- 6 对prmote reads bed算最大的50bp；
python ${bash_dir}/got_max50windw.py -i promoter_final_1000_reads_num.bed -p promoter_final_1000_stis.bed
sed -i '1 d' promoter_final_1000_stis.bed
awk 'OFS="\t" {print $2,$3,$4,$5,$6,$7,$1}' promoter_final_1000_stis.bed | sort -k 1,1 -k 2n,2 -k 3n,3 -o promoter_final_1000_stis.bed

#--- 7 使用R计算四类gene
mkdir -p four_class_gene
Rscript ${bash_dir}/four_class_promoter.r -b promoter_final_1000_stis.bed -g genebody_final_1000_reads_num.bed -d four_class_gene -p ${prefix} -n non_gene.bed -s genebody_final_-1000_reads_num.bed \
-u unique_short_gene.bed 

#--- 8 根据四类gene得到四类promoter
bash2_dir=/DATA/work/lbyybl/genomes/mm10/DNA_elements_made_by_Boyuan/got_elements_sh/
mkdir -p four_promoter
${bash2_dir}/gotpromoterformucsc.sh -g four_class_gene/${prefix}paused_active.bed -r 1000 -o four_promoter
mv four_promoter/promoter_merge_1000.bed four_promoter/paused_active_promoter_final_1000.bed
rm -f four_promoter/promoter_final_1000.bed
${bash2_dir}/gotpromoterformucsc.sh -g four_class_gene/${prefix}paused_unactive.bed -r 1000 -o four_promoter
mv four_promoter/promoter_merge_1000.bed four_promoter/paused_unactive_promoter_final_1000.bed
rm -f four_promoter/promoter_final_1000.bed
${bash2_dir}/gotpromoterformucsc.sh -g four_class_gene/${prefix}unpaused_active.bed -r 1000 -o four_promoter
mv four_promoter/promoter_merge_1000.bed four_promoter/unpaused_active_promoter_final_1000.bed
rm -f four_promoter/promoter_final_1000.bed
${bash2_dir}/gotpromoterformucsc.sh -g four_class_gene/${prefix}unpaused_unactive.bed -r 1000 -o four_promoter
mv four_promoter/promoter_merge_1000.bed four_promoter/unpaused_unactive_promoter_final_1000.bed
rm -f four_promoter/promoter_final_1000.bed
#--- 为了适应自己算loop算apa的程序，将上述文件的格式改一下
cd four_promoter
awk -v name=PAAC 'OFS="\t" {print $1,$2,$3,$4,name"-"$4}' paused_active_promoter_final_1000.bed >> paused_active_promoter_format.bed
awk -v name=PAUA 'OFS="\t" {print $1,$2,$3,$4,name"-"$4}' paused_unactive_promoter_final_1000.bed >> paused_unactive_promoter_format.bed
awk -v name=UPAC 'OFS="\t" {print $1,$2,$3,$4,name"-"$4}' unpaused_active_promoter_final_1000.bed >> unpaused_active_promoter_format.bed
awk -v name=UPUC 'OFS="\t" {print $1,$2,$3,$4,name"-"$4}' unpaused_unactive_promoter_final_1000.bed >> unpaused_unactive_promoter_format.bed
cat  *_promoter_format.bed >> promoter_format.bed

