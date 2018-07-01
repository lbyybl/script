#!/bin/bash

# Boyuan-Li

# Wed Jun 20 14:14:45 CST 2018

# This script is used to modified the name of the files in one specified directory;
# Though there mv to modified the name, but when we want to modified part of the name batchly,
# maybe mv will not very convenient, that times you can use this script

#set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -i <Input> -b <Modify_before_name> -a <Modify_after_name> "
	echo ""
	echo " -i STRING          [required]  the directory storing the file you want to change name"
	echo ""
	echo " -b STRING          [required]  the part of string in the name you want to change batchly, notice: only can batch!!"
	echo ""
	echo " -a STRING          [required]  the part of string you want to transfer into"
	echo ""
	echo " -h                 help"
	echo ""
	echo " This script is used to modified the name of the files in one specified directory;\
Though there mv to modified the name, but when we want to modified part of the name batchly\
 maybe mv will not very convenient, that times you can use this script"
#	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "i:b:a:h" optionName
do
	case $optionName in
		i) Input="$OPTARG";;
		b) Modify_before_name="$OPTARG";;
		a) Modify_after_name="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Input = "" ]]; then
	echo " -i the Input directory is needed "
	exit 1
elif [[ ! -d $Input ]]; then
	echo "$Input:   is not found"
	exit 2
fi


if [[ $Modify_before_name = "" ]]; then
	echo " -b the Modify_before_name STRING is needed "
	exit 1
fi


if [[ $Modify_after_name = "" ]]; then
	echo " -a the Modify_after_name STRING is needed "
	exit 1
fi

	ls $Input/* | parallel -j1 echo {/} '>>' $Input/__sample_md_name.list 
	sed -e 's/'"${Modify_before_name}"'/'"${Modify_after_name}"'/g' $Input/__sample_md_name.list >> $Input/__sample_md_af_name.list
	paste $Input/__sample_md_name.list $Input/__sample_md_af_name.list >> $Input/_trans_name.list
	cat $Input/_trans_name.list | parallel -C "\s" mv $Input/{1} $Input/{2}
	rm $Input/__sample_md_name.list $Input/_trans_name.list $Input/__sample_md_af_name.list
