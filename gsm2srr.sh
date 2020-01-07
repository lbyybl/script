#!/bin/bash

# Boyuan-Li

# Wed Aug  8 21:24:19 CST 2018

# It's used to convert the GSM TO SRR ID according to the SRA_Accessions.tab downloading from SRA ftp server(ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata)

set -euo pipefail
function helps
{
	echo ""
	echo -e "Usage: $0 [options] -g <Gsm_id> -f <SRA2GSM_file> "
	echo ""
	echo " -g STRING          [required] GSM ID"
	echo ""
	echo " -f STRING          [required] SRA_Accessions.tab file"
	echo ""
	echo " -h                 help"
	echo ""
	echo ""
}


if [ $# -eq 0 ]; then
	helps
	exit 0
fi
while getopts "g:f:h" optionName
do
	case $optionName in
		g) Gsm_id="$OPTARG";;
		f) SRA2GSM_file="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done



if [[ $Gsm_id = "" ]]; then
	echo " -g the Gsm_id STRING is needed "
	exit 1
fi


if [[ $SRA2GSM_file = "" ]]; then
	echo " -f the SRA2GSM_file file is needed "
	exit 1
elif [[ ! -f $SRA2GSM_file ]]; then
	echo "$SRA2GSM_file:   is not found"
	exit 2
fi


mediavar=$(egrep -w "${Gsm_id}" ${SRA2GSM_file} | egrep "SRX|SRR" | cut -f 1)
judgevar=$(echo -e "${mediavar}" | cut -b 1-3)
if [[ ${judgevar} = "SRR" ]]; then
	result=${mediavar}
else 
	result=$(egrep -w "${mediavar}" ${SRA2GSM_file} | cut -f 1 | grep "SRR")
fi

echo -e "${Gsm_id} was covert to ......"
echo -e "${result}"
