#!/bin/bash

# Boyuan-Li
# 2018/6/16
set -euo pipefail
# This script is used to add the judge for the option of getopt

##########
# common 3 type options: file, str or dir; so there are function to create the common
# conditonal statement
##########
mk_dir=""
function helps
{
        echo ""
        echo -e "Usage: $0 [required] -f <script you are operating> -t <variable type> [required] -n <variable name> -o <the option letter> -c <if not exit then create the dir>"
        echo ""
	echo " -f STRING          [required] the script you are operating"
        echo " -t STRING          [required] the variable type you want to create"
	echo " -n STRING          [required] the variable name you want"
	echo " -o STRING	  [required] the option litter is needed, and must be a letter "
	echo " -c STRING          [optional] if the dir you input don't exit, does need to create it?, your answer must be yes or no"
        echo ""

}
if [ $# -eq 0 ]; then
        helps
        exit 0
fi
while getopts "f:t:n:o:ch" optionName
do
        case "$optionName" in
		f) script="$OPTARG";;
                t) var_type="$OPTARG";;
                n) var_name="$OPTARG";;
		o) op_tion="$OPTARG";;
		c) mk_dir="$OPTARG";;
                h)
                        helps
                        exit 0
                        ;;
        esac
done

if [[ $script = "" ]]; then
        echo " -f the bash script file you are operation is required "
        exit 1
fi

if [[ $var_type = "" ]]; then
        echo " -t the var type must be file, string or dir "
        exit 1
fi

if [[ $var_name = "" ]]; then
        echo " -n the var name is needed "
        exit 1
fi

if [[ $op_tion = "" ]]; then
        echo " -o the option litter is needed, and must be a letter "
        exit 1
fi


if [[ $mk_dir = "" ]]; then
#        echo " -c mk_dir is need, must be yes or no!\
#	default no "
	mk_dir="no"
#        exit 1
fi


function filestatment
{
	op_tion=$1
	file_name=$2
	script=$3
	echo "" >> $script
	echo -e "if [[ \$$file_name = \"\" ]]; then" >> $script
	echo -e "\techo \" -$op_tion the $file_name file is needed \"" >> $script
	echo -e "\texit 1" >> $script
	echo -e "elif [[ ! -f \$$file_name ]]; then" >> $script
	echo -e "\techo \"\$$file_name:   is not found\"" >> $script
	echo -e "\texit 2" >> $script
	echo -e "fi" >> $script
	echo "" >> $script
}

function make_dir
{
	op_tion=$1
	dir_name=$2
	script=$3
#	ty_pe=$4
	echo "" >> $script
	echo -e "if [[ \$$dir_name = \"\" ]]; then" >> $script
	echo -e "\techo \" -$op_tion the $dir_name directory is needed \"" >> $script
	echo -e "\texit 1" >> $script
	echo -e "elif [[ ! -d \$$dir_name ]]; then" >> $script
	echo -e "\techo \"\$$dir_name:   is not found\"" >> $script
	echo -e "\tcreat dir" >> $script
	echo -e "mkdir \$$dir_name" >> $script
	#echo -e "\texit 2" 
	echo -e "fi" >> $script
	echo "" >> $script
}

function unmake_dir
{
	op_tion=$1
	dir_name=$2
	script=$3
	echo "" >> $script
	echo -e "if [[ \$$dir_name = \"\" ]]; then" >> $script
	echo -e "\techo \" -$op_tion the $dir_name directory is needed \"" >> $script
	echo -e "\texit 1" >> $script
	echo -e "elif [[ ! -d \$$dir_name ]]; then" >> $script
	echo -e "\techo \"\$$dir_name:   is not found\"" >> $script
	#echo -e "\tcreat dir" 
	#echo -e "mkdir \$$dir_name"
	echo -e "\texit 2" >> $script
	echo -e "fi" >> $script
	echo "" >> $script
}


function stringstatement
{
	op_tion=$1
	str_name=$2
	script=$3
	echo "" >> $script
	echo -e "if [[ \$$str_name = \"\" ]]; then" >> $script
	echo -e "\techo \" -$op_tion the $str_name STRING is needed \"" >> $script
	echo -e "\texit 1" >> $script
	#echo -e "elif [[ ! -f \$$file_name ]]; then"
	#echo -e "\techo \"\$$file_name:   is not found\""
	#echo -e "\texit 2"
	echo -e "fi" >> $script
	echo "" >> $script
}



if [[ $var_type = "file" ]]; then
	filestatment $op_tion $var_name $script
fi


if [[ $var_type = "string" ]]; then
	stringstatement $op_tion $var_name $script
fi


if [[ $var_type = "dir" && $mk_dir = "no" ]]; then
	unmake_dir $op_tion $var_name $script
fi

if [[ $var_type = "dir" && $mk_dir = "yes" ]]; then
        make_dir $op_tion $var_name $script
	
fi


