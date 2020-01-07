#!/bin/bash

# Boyuan-Li

# Thu Aug  9 18:51:56 CST 2018

# call compartment, resort the matrix and draw saddleplot
# need the script juicer2compartment.sh sort_matrix_compartment.R and draw_saddleplot.r

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -b <Bed_file_by_hicpro_corresponding_to_resoluton> -j <Hic_file_by_juicer> -r <Resolution> -g <Genome_version> -t <Juicer_type> -o <Output_filename> "
	echo ""
	echo " -b STRING          [required] Bed_file_by_hicpro_corresponding_to_resoluton"
	echo ""
	echo " -j STRING          [required] .hic file produced by juicer pre"
	echo ""
	echo " -r STRING          [required] bin size"
	echo ""
	echo " -g STRING          [required] The genome version you used, such as mm10,mm9,hg19 et al."
	echo ""
	echo " -t STRING          [required] the juicer observed or observed/expected"
	echo ""
	echo " -o STRING          [required] the output pdf file name of saddleplot"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "b:j:r:g:t:o:h" optionName
do
	case $optionName in
		b) Bed_file_by_hicpro_corresponding_to_resoluton="$OPTARG";;
		j) Hic_file_by_juicer="$OPTARG";;
		r) Resolution="$OPTARG";;
		g) Genome_version="$OPTARG";;
		t) Juicer_type="$OPTARG";;
		o) Output_filename="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Bed_file_by_hicpro_corresponding_to_resoluton = "" ]]; then
	echo " -b the Bed_file_by_hicpro_corresponding_to_resoluton file is needed "
	exit 1
elif [[ ! -f $Bed_file_by_hicpro_corresponding_to_resoluton ]]; then
	echo "$Bed_file_by_hicpro_corresponding_to_resoluton:   is not found"
	exit 2
fi


if [[ $Hic_file_by_juicer = "" ]]; then
	echo " -j the Hic_file_by_juicer file is needed "
	exit 1
elif [[ ! -f $Hic_file_by_juicer ]]; then
	echo "$Hic_file_by_juicer:   is not found"
	exit 2
fi


if [[ $Resolution = "" ]]; then
	echo " -r the Resolution STRING is needed "
	exit 1
fi


if [[ $Genome_version = "" ]]; then
	echo " -g the Genome_version STRING is needed "
	exit 1
fi


if [[ $Juicer_type = "" ]]; then
	echo " -t the Juicer_type STRING is needed "
	exit 1
fi


if [[ $Output_filename = "" ]]; then
	echo " -o the Output_filename STRING is needed "
	exit 1
fi

Hicfile=${Hic_file_by_juicer##*/}
Hicfilename=${Hicfile%.*}
Hic_file_by_juicer=$(readlink -e ${Hic_file_by_juicer})
Bed_file_by_hicpro_corresponding_to_resoluton=$(readlink -e ${Bed_file_by_hicpro_corresponding_to_resoluton})
Pwd_path=$(pwd)
mkdir -p ${Pwd_path}/_zscore/${Resolution}/  ${Pwd_path}/_ploter/${Resolution}/
#--- call compartment
	juicer2compartment.sh -b $Bed_file_by_hicpro_corresponding_to_resoluton -j $Hic_file_by_juicer -r $Resolution -o ${Pwd_path}/_zscore/${Resolution}/ -g $Genome_version -t $Juicer_type
#--- got observed/expected matrix
	echo {{1..19},x,y} | xargs -n1 | parallel java -jar /DATA/work/lbyybl/tools/juicer_fml/juicer_tools.1.8.9_jcuda.0.8.jar dump oe KR $Hic_file_by_juicer {} {} BP $Resolution ${Pwd_path}/_ploter/${Resolution}/chr{}.oe
#--- resort the matrix
	cd ${Pwd_path}/_ploter/${Resolution}/
	echo chr{1..19} | xargs -n1 | parallel Rscript /home/boyuanli/bashscript/bin/sort_matrix_compartment.R -c ${Pwd_path}/_zscore/${Resolution}/{}_*_dense.addedHeaders.zScore.compartments -m {}.oe -o {}.oe.resort -r $Resolution
#--- draw saddleplot
	Rscript /home/boyuanli/bashscript/bin/saddleplot/draw_saddleplot.r -p ./ -o ${Pwd_path}/$Output_filename
#--- rm dir
	#rm -rf ${Pwd_path}/_zscore/${Resolution}/ ${Pwd_path}/_ploter/${Resolution}/	
