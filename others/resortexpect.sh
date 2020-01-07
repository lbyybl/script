#!/bin/bash

# Boyuan-Li

# Fri Jul 13 21:35:06 CST 2018

# this script is used to resort the matrix according to the eg1 produced by the cworld (x2compartment.pl);
# means change the bin number according to the eg1

set -eu 
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -z <Zscore_file_produced_by_cworld> -m <Correspording_spare_matrix> -o <Output_spare_file> -c <chr name>"
	echo ""
	echo " -z STRING          [required] the Zscore_file_produced_by_cworld"
	echo ""
	echo " -m STRING          [required] the Correspording_spare_matrix produced by HiC-Pro"
	echo ""
	echo " -o STRING          [required] the output file name"
	echo ""
	echo " -c STRING          [optional] the chr name, if you want to see the whole genome, please NULL"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
Chr_name=""
while getopts "z:m:o:c:h" optionName
do
	case $optionName in
		z) Zscore_file_produced_by_cworld="$OPTARG";;
		m) Correspording_spare_matrix="$OPTARG";;
		o) Output_spare_file="$OPTARG";;
		c) Chr_name="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Zscore_file_produced_by_cworld = "" ]]; then
	echo " -z the Zscore_file_produced_by_cworld file is needed "
	exit 1
elif [[ ! -f $Zscore_file_produced_by_cworld ]]; then
	echo "$Zscore_file_produced_by_cworld:   is not found"
	exit 2
fi


if [[ $Correspording_spare_matrix = "" ]]; then
	echo " -m the Correspording_spare_matrix file is needed "
	exit 1
elif [[ ! -f $Correspording_spare_matrix ]]; then
	echo "$Correspording_spare_matrix:   is not found"
	exit 2
fi


if [[ $Output_spare_file = "" ]]; then
	echo " -o the Output_spare_file STRING is needed "
	exit 1
fi

if [[ $Chr_name = "" ]]; then
	echo " The resort is used to the whole genome "
	Chr_name="wholegenome"
fi


#########################################################################################
# this is used to change the file produce by the  cworld into soted file;$1 is the file 
# produced by the cworld, the $2 is the bed file in the raw dir produced by hic-pro
#########################################################################################

cat $Zscore_file_produced_by_cworld | cut -f 1-3,6 >> ${Chr_name}_chrbineg1_
sed -i '/^#chr.*$/d' ${Chr_name}_chrbineg1_
#sed -i -e 's/\tnan\t/\t0.00\t/g' -e 's/\tnan\t/\t0.00\t/g' -e 's/\tnan$/\t0.00/g' ${Chr_name}_chrbineg1_
awk '{++cnt;printf "%s \t %i \t %i \t %f \t %i \n",$1,$2,$3,$4,cnt}' ${Chr_name}_chrbineg1_ >> ${Chr_name}_chrbineg1plus_
#--- find the EG1==NAN bin
number=$(awk '{if ($4=="nan") print $5}' ${Chr_name}_chrbineg1plus_ | xargs | sed 's/\s/|/g' )
egrep -wv "${number}" $Correspording_spare_matrix >> ${Chr_name}_matrix_rm_nancorbin
sort -rg -k4 ${Chr_name}_chrbineg1plus_ | cut -f 1-3,5 >> ${Chr_name}_sorted-eg1_

##########################################################################################
# below is change the sorted file prduct by the cworld according to the eg1; _sorted-eg1_ is the 
# sorted file, $3 is the spare matrix, $4 is the output-spare file;
##########################################################################################
#bedname=${Correspording_bed_file##*/}
awk '{++cn;printf "%s \t %i \t %i \t %i \t %i \n",$1,$2,$3,$4,cn}' ${Chr_name}_sorted-eg1_  >> ${Chr_name}_addnewcount_
cat ${Chr_name}_addnewcount_ | cut -f 4-5 >> ${Chr_name}_correspondingnumber_
echo "#!/bin/bash" >> ${Chr_name}_numberarray_
awk '{printf "arr[%i]=%i\n",$1,$2}' ${Chr_name}_correspondingnumber_ >> ${Chr_name}_numberarray_
chmod +x ${Chr_name}_numberarray_


echo "#!/bin/bash" >> ${Chr_name}_change_bin_number_
echo -e "source ${Chr_name}_numberarray_" >> ${Chr_name}_change_bin_number_
awk '{printf "echo -e \"${arr[%i]}\t${arr[%i]}\t%f\"\n",$1,$2,$3}' ${Chr_name}_matrix_rm_nancorbin >> ${Chr_name}_change_bin_number_
chmod +x ${Chr_name}_change_bin_number_
./${Chr_name}_change_bin_number_ >> ${Output_spare_file}
#cat ${Chr_name}_addnewcount_ | cut -f 1-3,5 >> ${bedname%.*}_resort.bed # maybe this file is no use;
rm ${Chr_name}_addnewcount_ ${Chr_name}_chrbineg1_  ${Chr_name}_chrbineg1plus_  ${Chr_name}_correspondingnumber_ ${Chr_name}_numberarray_ ${Chr_name}_change_bin_number_ ${Chr_name}_sorted-eg1_ 
rm ${Chr_name}_matrix_rm_nancorbin


