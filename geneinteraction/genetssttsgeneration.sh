#!/bin/bash

# Boyuan-Li

# Fri Sep 21 11:50:49 CST 2018

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -l <Gene_length> -g <Gene_file> -i <Insulator_file> -e <Enhancer_file> "
	echo ""
	echo " -l STRING          [required] "
	echo ""
	echo " -g STRING          [required] "
	echo ""
	echo " -i STRING          [required] "
	echo ""
	echo " -e STRING          [required] "
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "l:g:i:e:h" optionName
do
	case $optionName in
		l) Gene_length="$OPTARG";;
		g) Gene_file="$OPTARG";;
		i) Insulator_file="$OPTARG";;
		e) Enhancer_file="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Gene_length = "" ]]; then
	echo " -l the Gene_length STRING is needed "
	exit 1
fi


if [[ $Gene_file = "" ]]; then
	echo " -g the Gene_file file is needed "
	exit 1
elif [[ ! -f $Gene_file ]]; then
	echo "$Gene_file:   is not found"
	exit 2
fi


if [[ $Insulator_file = "" ]]; then
	echo " -i the Insulator_file file is needed "
	exit 1
elif [[ ! -f $Insulator_file ]]; then
	echo "$Insulator_file:   is not found"
	exit 2
fi


if [[ $Enhancer_file = "" ]]; then
	echo " -e the Enhancer_file file is needed "
	exit 1
elif [[ ! -f $Enhancer_file ]]; then
	echo "$Enhancer_file:   is not found"
	exit 2
fi

		awk -v len="${Gene_length}000" '{lent=$5-$4;if (lent > len) print $0}' $Gene_file >> genelen${Gene_length}k.bed
		awk '{if ($3=="-") print $0}' genelen${Gene_length}k.bed >> Ngenelen${Gene_length}k.bed
		awk '{if ($3=="+") print $0}' genelen${Gene_length}k.bed >> Pgenelen${Gene_length}k.bed
		awk '{st=$4+1000;en=$5-1000;printf "%s\t%s\t%s\t%s\t%s\n",$2,st,en,$3,$1}' Ngenelen${Gene_length}k.bed | awk '{printf "%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,"Ge"$5}' >> gene_body/Ngene_body${Gene_length}K.bed
		awk '{st=$4+1000;en=$5-1000;printf "%s\t%s\t%s\t%s\t%s\n",$2,st,en,$3,$1}' Pgenelen${Gene_length}k.bed | awk '{printf "%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,"Ge"$5}' >> gene_body/Pgene_body${Gene_length}K.bed
		#--- tss
		awk '{st=$5-1000;en=$5+1000;printf "%s\t%s\t%s\t%s\t%s\n",$2,st,en,$3,$1}' Ngenelen${Gene_length}k.bed | awk '{cnt++;printf "%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,"NTSS"cnt}' >> TSS/NTSS${Gene_length}K.bed
		awk '{st=$4-1000;en=$4+1000;printf "%s\t%s\t%s\t%s\t%s\n",$2,st,en,$3,$1}' Pgenelen${Gene_length}k.bed | awk '{cnt++;printf "%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,"PTSS"cnt}' >> TSS/PTSS${Gene_length}K.bed
		#--- tts
		awk '{st=$4-1000;en=$4+1000;printf "%s\t%s\t%s\t%s\t%s\n",$2,st,en,$3,$1}' Ngenelen${Gene_length}k.bed | awk '{cnt++;printf "%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,"NTTS"cnt}' >> TTS/NTTS${Gene_length}K.bed
		awk '{st=$5-1000;en=$5+1000;printf "%s\t%s\t%s\t%s\t%s\n",$2,st,en,$3,$1}' Pgenelen${Gene_length}k.bed | awk '{cnt++;printf "%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,"PTTS"cnt}' >> TTS/PTTS${Gene_length}K.bed
		cat $Enhancer_file $Insulator_file gene_body/*${Gene_length}K.bed TSS/*${Gene_length}K.bed TTS/*${Gene_length}K.bed >> final_gene_element${Gene_length}k.bed
		