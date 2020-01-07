#!/bin/bash

# Boyuan_Li

# Mon Jan 21 20:53:30 2019

set -euo pipefail
function helps
{
	echo ""
	echo -e "\tIt's used to calculate the APA used loop_file called by hichipper; Gene_bed file which interaction \n \
	interested by you and a prefix will be used by the output file "
	echo ""
	echo -e "Usage: $0 [options] -l <loop_file> -g <gene_bed_file> -p <Prefix> "
	echo ""
	echo " -l	file         	[required] "
	echo ""
	echo " -g	file         	[required] "
	echo ""
	echo " -p	string       	[required] "
	echo ""
	echo " -h			help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "l:g:p:h" optionName
do
	case $optionName in
		l) loop_file="$OPTARG";;
		g) gene_bed_file="$OPTARG";;
		p) Prefix="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $loop_file = "" ]]; then
	echo "the $loop_file file is needed "
	exit 1
elif [[ ! -f $loop_file ]]; then
	echo "$loop_file:   is not found"
	exit 2
fi

if [[ $gene_bed_file = "" ]]; then
	echo "the $gene_bed_file file is needed "
	exit 1
elif [[ ! -f $gene_bed_file ]]; then
	echo "$gene_bed_file:   is not found"
	exit 2
fi

if [[ $Prefix = "" ]]; then
	echo " the $Prefix string is needed "
	exit 1
fi

elements=$(awk '{s=substr($5,1,4);print s}' ${gene_bed_file} | sort | uniq | xargs)
#--- 7. 对于contact strength的图
	#---将loop bedpe文件修改成/home/boyuanli/bashscript/bin/HiCcompare/find_region.sh可以使用的格式；
	awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,$8,"1","2","3","4","5","6","7","8","9","10"}' ${loop_file}  >> ${Prefix}_format.bed
	#---使用/home/boyuanli/bashscript/bin/HiCcompare/find_region.sh找出loop相连接的原件
	/home/boyuanli/bashscript/bin/HiCcompare/find_region.sh -g ${gene_bed_file} -d ${Prefix}_format.bed -o ${Prefix}_all_interaction.bed
	#--- got interaction with elements
	awk '{if ($23!="." && $28 != ".") print $0}' ${Prefix}_all_interaction.bed >> ${Prefix}_focus_interaction.bed
	#--- 2. 分开各种互做的种类
	parallel -k "awk -v loc1=\"{1}\" -v loc2=\"{2}\" '{e1=substr(\$23,1,4);e2=substr(\$28,1,4);if ((e1==loc1 && e2==loc2) || (e1==loc2 && e2==loc1)) print \$0}' ${Prefix}_focus_interaction.bed >> ${Prefix}_{1}_{2}.bed" ::: ${elements}  :::  ${elements}
	