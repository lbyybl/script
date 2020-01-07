#!/bin/bash

# Boyuan_Li

# Thu Jan 24 00:11:54 2019

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -w <wt_loop> -k <ko_loop> -g <gene_bed_file> -c <wt_hic> -e <ko_hic> -p <prefix> "
	echo ""
	echo " -w	file         	[required] "
	echo ""
	echo " -k	file         	[required] "
	echo ""
	echo " -g	file         	[required] "
	echo ""
	echo " -c	file         	[required] "
	echo ""
	echo " -e	file         	[required] "
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
while getopts "w:k:g:c:e:p:h" optionName
do
	case $optionName in
		w) wt_loop="$OPTARG";;
		k) ko_loop="$OPTARG";;
		g) gene_bed_file="$OPTARG";;
		c) wt_hic="$OPTARG";;
		e) ko_hic="$OPTARG";;
		p) prefix="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done


if [[ $wt_loop = "" ]]; then
	echo "the $wt_loop file is needed "
	exit 1
elif [[ ! -f $wt_loop ]]; then
	echo "$wt_loop:   is not found"
	exit 2
fi

if [[ $ko_loop = "" ]]; then
	echo "the $ko_loop file is needed "
	exit 1
elif [[ ! -f $ko_loop ]]; then
	echo "$ko_loop:   is not found"
	exit 2
fi

if [[ $gene_bed_file = "" ]]; then
	echo "the $gene_bed_file file is needed "
	exit 1
elif [[ ! -f $gene_bed_file ]]; then
	echo "$gene_bed_file:   is not found"
	exit 2
fi

if [[ $wt_hic = "" ]]; then
	echo "the $wt_hic file is needed "
	exit 1
elif [[ ! -f $wt_hic ]]; then
	echo "$wt_hic:   is not found"
	exit 2
fi

if [[ $ko_hic = "" ]]; then
	echo "the $ko_hic file is needed "
	exit 1
elif [[ ! -f $ko_hic ]]; then
	echo "$ko_hic:   is not found"
	exit 2
fi

if [[ $prefix = "" ]]; then
	echo " the $prefix string is needed "
	exit 1
fi
wt_hic=$(readlink -e $wt_hic)
ko_hic=$(readlink -e $ko_hic)
wt_loop=$(readlink -e $wt_loop)
ko_loop=$(readlink -e $ko_loop)
 mkdir ko 
	cd ko
	/home/boyuanli/bashscript/bin/loop/contact_strength/classify_loop.sh -l ${ko_loop} -g ${gene_bed_file} -p ${prefix}_ko
	cd ../
	mkdir wt
	cd wt
	/home/boyuanli/bashscript/bin/loop/contact_strength/classify_loop.sh -l ${wt_loop} -g ${gene_bed_file} -p ${prefix}_wt

 cd ../	
#---merge loop
function merge_file
{
	wt_file=$1
	ko_file=$2
	out_file=$3
	cut -f 1-6 $wt_file >> $out_file
	cut -f 1-6 $ko_file >> $out_file
	sort $out_file | uniq | awk 'OFS="\t"{print $0,".","."}' | sort -o $out_file
}

# function merge_file
# {
	# wt_file=$1
	# ko_file=$2
	# out_file=$3
	# cut -f 1-8 $wt_file | sort | uniq >> $out_file
# }
	elements=$(awk '{s=substr($5,1,4);print s}' ${gene_bed_file} | sort | uniq | xargs)
	#merge_file wt/hichip_wt_Supe_Supe.bed ko/hichip_ko_Supe_Supe.bed hichip_all_Supe_Supe.bed
	for i in $(parallel echo {1}_{2}.bed ::: ${elements}  :::  ${elements})
	do 
		merge_file wt/${prefix}_wt_${i} ko/${prefix}_ko_${i} all_${i}
	done	
		
	#parallel merge_file wt/${prefix}_wt_{1}_{2}.bed ko/${prefix}_ko_{1}_{2}.bed all_{1}_{2}.bed ::: ${elements}  :::  ${elements}
	
	cd wt 
	#--- 3. APA
	parallel -j 5 -k "java -jar /DATA/work/lbyybl/tools/juicer_fml/juicer_tools.1.8.9_jcuda.0.8.jar apa -u -w 10 -k VC_SQRT -r 10000 $wt_hic ../all_{1}_{2}.bed ${prefix}_wt_{1}_{2}" ::: ${elements}  :::  ${elements} &
	
	cd ../
	cd ko
	#--- 3. APA
	parallel -j 5 -k "java -jar /DATA/work/lbyybl/tools/juicer_fml/juicer_tools.1.8.9_jcuda.0.8.jar apa -u -w 10 -k VC_SQRT -r 10000 $ko_hic ../all_{1}_{2}.bed ${prefix}_ko_{1}_{2}" ::: ${elements}  :::  ${elements} 
	