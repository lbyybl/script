#!/bin/bash

# Boyuan-Li

# 2018/6/16

# This script is used to add header for the bash script

set -euo pipefail

auther='Boyuan-Li'
function helps
{
	echo ""
	echo -e "Usage: $0 [required] -f <shell script you will add header> [required] -u <auther> [required]"
	echo "-f STRING	[required] the shell script you will add header"
	echo "-u STRING	[optional] the auther of the script; defult Boyuan-Li"
	echo ""
}

if [ $# -eq 0 ]; then
        helps
        exit 0
fi

while getopts "f:uh" optionName
do
        case "$optionName" in
                f) script="$OPTARG";;
                u) auther="$OPTARG";;
                h)
                        helps
                        exit 0
                        ;;
        esac
done

if [[ $auther = "" ]]; then
	echo " -u the auther file is needed "
	exit 1
fi


if [[ $script = "" ]]; then
	echo " -f the script file is needed "
	exit 1
#elif [[ ! -f $script ]]; then
#	echo "$script:   is not found"
#	exit 2
fi
	
#	echo "set -euo pipefail" >> $script
	echo "#!/bin/bash" >> $script
	echo "" >> $script
	echo -e "# $auther" >> $script
	echo "" >> $script
	echo -e "# `date`" >> $script
	echo "" >> $script
	echo "set -euo pipefail" >> $script
	echo ""
