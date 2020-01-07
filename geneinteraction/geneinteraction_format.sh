#!/bin/bash

# Boyuan-Li

# Tue Sep 11 09:20:07 CST 2018

# it's used to got the final interaction format

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -i <Input> -l <Namelist> -o <Output_file> "
	echo ""
	echo " -i STRING          [required] "
	echo ""
	echo " -l STRING          [required] namelist file"
	echo ""
	echo " -o STRING          [required] "
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}

#--- namelist format like fellow 第一列必须是四个字母
# Insu    Insulator
# INT_    Typical_enhancer
# NPro    Promoter
# PPro    Promoter
# Stan    Gene_body
# Supe    Super_enhancer

if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "i:l:o:h" optionName
do
	case $optionName in
		i) Input="$OPTARG";;
		l) Namelist="$OPTARG";;
		o) Output_file="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Input = "" ]]; then
	echo " -i the Input file is needed "
	exit 1
elif [[ ! -f $Input ]]; then
	echo "$Input:   is not found"
	exit 2
fi


if [[ $Namelist = "" ]]; then
	echo " -l the Namelist file is needed "
	exit 1
elif [[ ! -f $Namelist ]]; then
	echo "$Namelist:   is not found"
	exit 2
fi


if [[ $Output_file = "" ]]; then
	echo " -o the Output_file STRING is needed "
	exit 1
fi

Prefix=${Output_file%.*}
cut -f 1 $Namelist | awk '{st=substr($1,1,4); print st}' | sort | uniq >> namelist_${Prefix}
parallel echo ::: `cat namelist_${Prefix}` ::: `cat namelist_${Prefix}` >> interaction_${Prefix}.list
cat interaction_${Prefix}.list | parallel -C " " -k "awk -v str1={1} -v str2={2} '{st=substr(\$1,1,4); en=substr(\$2,1,4);if ((st==str1 && en==str2) || (en==str1 && st==str2)) print \$0}' ${Input} >> ${Prefix}_{1}_-_{2}.bed"
ls ${Prefix}* | parallel -k -v "Rscript /home/boyuanli/bashscript/bin/geneinteraction/sum.r -i {}" >> interaction_${Prefix}.stastic
sed -e 's/_-_/\t/g' -e "s/Rscript.*${Prefix}_//g" -e 's/.bed//g' -e 's/\[1\] //g'  interaction_${Prefix}.stastic | xargs -l2 >> interaction_${Prefix}_format.stastic
cp $Namelist namelist_${Prefix}2
for i in `cut -f 1 $Namelist | awk '{st=substr($1,2,4); print st}' | sort | uniq -c | awk '{if ($1 > 1) print $2}' `; do
	/home/boyuanli/bashscript/bin/geneinteraction/merge.sh -i interaction_${Prefix}_format.stastic -s ${i} -o interaction_${Prefix}_format.stastic ;
	var=$(egrep -w "^.${i}" namelist_${Prefix}2 | cut -f 2 | uniq)
	sed -i "/^.${i}\b/d" namelist_${Prefix}2
	echo -e "T${i}\t${var}" >> namelist_${Prefix}2
done
	
cat namelist_${Prefix}2 | parallel -j 1 -C '\t' sed -i "s/{1}/{2}/g" interaction_${Prefix}_format.stastic
mv interaction_${Prefix}_format.stastic $Output_file
#rm namelist_${Prefix} interaction_${Prefix}.list namelist_${Prefix}2

