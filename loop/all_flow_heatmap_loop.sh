#!/bin/bash

# Boyuan-Li

# Fri Nov 23 18:29:35 CST 2018

# this file interagate all flow to draw heatmap whith loop bedpe file produced by hichipper 
# include the region_interaction.sh geneinteraction_format.sh stasinterval.sh gotnormalize.sh

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -g <Genome_bed_file> -w <Wt_loop_bedpe> -k <Ko_loop_bedpe> -1 <Namelist1> -2 <Namelist2> -o <Output_prefilx> "
	echo ""
	echo " -g STRING          [required] Gene_bed_file can be file merged by promoter enhancer gff file"
	echo ""
	echo " -w STRING          [required] wt Loop files called by hichipper"
	echo ""
	echo " -k STRING          [required] ko Loop files called by hichipper"
	echo ""
	echo " -1 STRING          [required] namelist1 file"
	echo ""
	echo " -2 STRING          [required] namelist2 file"
	echo ""
	echo " -o STRING          [required] Output_file name and the prefix also will be used and the dir will be used as output dir"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "g:w:k:1:2:o:h" optionName
do
	case $optionName in
		g) Genome_bed_file="$OPTARG";;
		w) Wt_loop_bedpe="$OPTARG";;
		k) Ko_loop_bedpe="$OPTARG";;
		1) Namelist1="$OPTARG";;
		2) Namelist2="$OPTARG";;
		o) Output_prefilx="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Genome_bed_file = "" ]]; then
	echo " -g the Genome_bed_file file is needed "
	exit 1
elif [[ ! -f $Genome_bed_file ]]; then
	echo "$Genome_bed_file:   is not found"
	exit 2
fi

if [[ $Wt_loop_bedpe = "" ]]; then
	echo " -k the Wt_loop_bedpe file is needed "
	exit 1
elif [[ ! -f $Wt_loop_bedpe ]]; then
	echo "$Wt_loop_bedpe:   is not found"
	exit 2
fi

if [[ $Ko_loop_bedpe = "" ]]; then
	echo " -k the Ko_loop_bedpe file is needed "
	exit 1
elif [[ ! -f $Ko_loop_bedpe ]]; then
	echo "$Ko_loop_bedpe:   is not found"
	exit 2
fi


if [[ $Namelist1 = "" ]]; then
	echo " -1 the Namelist1 file is needed "
	exit 1
elif [[ ! -f $Namelist1 ]]; then
	echo "$Namelist1:   is not found"
	exit 2
fi


if [[ $Namelist2 = "" ]]; then
	echo " -2 the Namelist2 file is needed "
	exit 1
elif [[ ! -f $Namelist2 ]]; then
	echo "$Namelist2:   is not found"
	exit 2
fi


if [[ $Output_prefilx = "" ]]; then
	echo " -o the Output_prefilx STRING is needed "
	exit 1
fi


# namelist1 like below 

# Enha    TE
# GeEN    GB
# Insu    IN
# NTSS    TSS
# NTTS    TTS
# PTSS    TSS
# PTTS    TTS
# Supe    SE
# NONG    NG

# namelist2 linke bewlow

# Enha    TE
# GeEN    GB
# Insu    IN
# TTS     TTS
# TSS     TSS
# Supe    SE
# NONG    NG




#--- get abs path

		Gene_file=$(readlink -e $Genome_bed_file)
		wt_loop_file=$(readlink -e $Wt_loop_bedpe)
		ko_loop_file=$(readlink -e $Ko_loop_bedpe)
		namelist1=$(readlink -e $Namelist1)
		namelist2=$(readlink -e $Namelist2)
		Output_dir="./"
		# Output_dir=${Output_prefilx%%/*}
		# dir="./"
        # dir=$(echo -e "${Output_prefilx}" | grep "/")
        # if [[ $dir != "" && $dir != "./" ]]; then
                        # Output_dir=${Output_prefilx%/*}
                        # if [[ ! -d $Output_dir ]]; then
                                # mkdir -p ${Output_dir}
                        # fi
        # fi

		Output_dir=$(readlink -e ${Output_dir})
		medfix=$(date | sed -e 's/ /./g' -e 's/:/./g')
#--- pipeline
	#--- generate interaction file
	mkdir -p tmp.${medfix}
		#--- wt
		/home/boyuanli/bashscript/bin/loop/region_interaction.sh -g ${Gene_file} -l ${wt_loop_file} -o tmp.${medfix}/wt.interaction.${medfix}.bed
		#--- ko
		/home/boyuanli/bashscript/bin/loop/region_interaction.sh -g ${Gene_file} -l ${ko_loop_file} -o tmp.${medfix}/ko.interaction.${medfix}.bed
		#rm -f 
	#--- generate format
	mkdir -p tmp.${medfix}/format
		cd tmp.${medfix}/format
		#--- wt
		/home/boyuanli/bashscript/bin/geneinteraction/geneinteraction_format.sh -i ../wt.interaction.${medfix}.bed -l ${namelist1} -o wt.format.bed
		#--- ko
		/home/boyuanli/bashscript/bin/geneinteraction/geneinteraction_format.sh -i ../ko.interaction.${medfix}.bed -l ${namelist1} -o ko.format.bed
		
		cd -
		
	#--- calculate region	
	mkdir -p tmp.${medfix}/format/region
	cd tmp.${medfix}/format/region
		/home/boyuanli/bashscript/bin/geneinteraction/stasinterval.sh -n ${namelist2} -i ${Gene_file} -o region.bed
	cd -	
	#--- generate file with region
		mkdir -p tmp.${medfix}/format/final
		#--- wt
		join -1 1 -2 1 <(sort -k 1,1 -k 2,2 tmp.${medfix}/format/wt.format.bed) <(sort -k1 tmp.${medfix}/format/region/region.bed) | sort -k 2,2 -k 1,1 -o tmp.${medfix}/format/final/wt.final.bed
		join -1 2 -2 1 tmp.${medfix}/format/final/wt.final.bed <(sort -k1 tmp.${medfix}/format/region/region.bed) | sort -o tmp.${medfix}/format/final/wt.final.bed
		#--- ko
		join -1 1 -2 1 <(sort -k 1,1 -k 2,2 tmp.${medfix}/format/ko.format.bed) <(sort -k1 tmp.${medfix}/format/region/region.bed) | sort -k 2,2 -k 1,1 -o tmp.${medfix}/format/final/ko.final.bed
		join -1 2 -2 1 tmp.${medfix}/format/final/ko.final.bed <(sort -k1 tmp.${medfix}/format/region/region.bed) | sort -o tmp.${medfix}/format/final/ko.final.bed
		
	#--- finally
		cp tmp.${medfix}/format/final/ko.final.bed ${Output_dir}/${Output_prefilx}_ko.bed
		cp tmp.${medfix}/format/final/wt.final.bed ${Output_dir}/${Output_prefilx}_wt.bed
		rm -rf tmp.${medfix}
		