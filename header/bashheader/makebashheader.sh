#!/bin/bash

# Boyuan-Li

# Sat Jun 16 22:26:50 CST 2018

# this scipt is used to cread the getoption accordiong to the parameter you give;
# every parameter format should option,option_name,type; such as i,Input_dir,dir;
# option is the letter will be regard option, option_name is name you give the option;
# type is dir,string,file
# this script need to in the same dir with the getopthelper.sh and headhelper.sh

k=""
for k in {A..G}
do
	eval ${k}_par=""
done

z=" "
for z in {A..G}
do
#        if [[ ${z}_par != "" ]]; then
#               do
        eval ${z}_option=""  ## very amazing, yestaday if no eval it' doesn't work, but tonday for some var if with eval will not work
        eval ${z}_var_name=""
        eval ${z}_op_type=""
#        fi
done


function echo_help
{
	option=$1
	echo -e " -$option	[optional] the $2 optional"
}

function helps
{
	echo ""
	echo -e "Usage: $0 [options] -s <script> -a <option,option_name,type> "
	echo " -s STRING          [required] the script you want to operating"
	echo "all the parameter format should be option,option_name,type"
	echo_help a first
	echo_help b second
	echo_help c third
	echo_help d forth
	echo_help e fifth
	echo_help f sixth
	echo_help g seventh
	echo " -h                 help"
	echo ""
	echo ""
}

if [ $# -eq 0 ]; then
	helps
	exit 0
fi

while getopts "s:a:b:c:d:e:f:g:h" optionName
do
	case "$optionName" in 
		s) script="$OPTARG";;
		a) A_par="$OPTARG";;
		b) B_par="$OPTARG";;
		c) C_par="$OPTARG";;
		d) D_par="$OPTARG";;
		e) E_par="$OPTARG";;
		f) F_par="$OPTARG";;
		g) G_par="$OPTARG";;
		h)
			helps
			exit 0
			;;
	esac
done

if [[ $script = "" ]]; then
	echo " -s <script you want operate> is required "
	exit 1
fi

i=" "
for i in {A..G}
do
#	eval parame_ter=${i}_par
	if [[ $(eval echo -e "\$${i}_par") != "" ]]; then
#		do
#		eval parame_ter=${i}_par 
		eval ${i}_option=$(eval echo -e "\$${i}_par" | cut -d "," -f 1)
		eval ${i}_var_name=$(eval echo -e "\$${i}_par" | cut -d "," -f 2)
		eval ${i}_op_type=$(eval echo -e "\$${i}_par" | cut -d "," -f 3)
	fi
done

#########
# print help function
#########
p=" "
usage=""
getopt_par=""
for p in {A..G}
do
	eval _opt=$(eval echo -e "\$${p}_option")
	eval _varname=$(eval echo -e "\$${p}_var_name")
	if [[ $_opt != "" ]]; then
		usage=${usage}"-"${_opt}" ""<"${_varname}">"" "     ###very amazing
		eval getopt_par=${getopt_par}${_opt}:
	fi
done

	headhelper.sh -f $script
	echo -e "function helps" >> $script
	echo -e "{" >> $script
	echo -e "\techo \"\"" >> $script
	echo -e "\techo -e \"Usage: \$0 [options] ${usage}\"" >> $script
	echo -e "\techo \"\"" >> $script
	
q=" "
for q in {A..G}
do 
	eval _opt=$(eval echo -e "\$${q}_option")
	eval _varname=$(eval echo -e "\$${q}_var_name")
	if [[ $_opt != "" ]]; then
		echo -e "\techo \" -${_opt} STRING          [required] \"" >> $script
		echo -e "\techo \"\"" >> $script
	fi
done

	echo -e "\techo \" -h                 help\"" >> $script
	echo -e "\techo \"\"" >> $script
	echo -e "\techo \"\"" >> $script
	echo -e "}" >> $script
	echo -e "" >> $script
	echo -e "" >> $script


#########
# print the getopt
#########
	

	echo -e "if [ \$# -eq 0 ]; then" >> $script
	echo -e "\thelps" >> $script
	echo -e "\texit 0" >> $script
	echo -e "fi" >> $script

	echo -e ""
	echo -e ""

## print getopt

	echo -e "while getopts \"${getopt_par}h\" optionName" >> $script
	echo -e "do" >> $script
	echo -e "\tcase "\$optionName" in" >> $script

r=""
for r in {A..G}
do
	eval _opt=$(eval echo -e "\$${r}_option")
	eval _varname=$(eval echo -e "\$${r}_var_name")
	if [[ $_opt != "" ]]; then
		echo -e "\t\t${_opt}) ${_varname}=\"\$OPTARG\";;" >> $script
#	echo -e "\t\t${B_option}) ${B_var_name}=\"\$OPTARG\";;" >> $script
#	echo -e "\t\t${C_option}) ${C_var_name}=\"\$OPTARG\";;" >> $script
#	echo -e "\t\t${D_option}) ${D_var_name}=\"\$OPTARG\";;" >> $script
#	echo -e "\t\t${E_option}) ${E_var_name}=\"\$OPTARG\";;" >> $script
#	echo -e "\t\t${F_option}) ${F_var_name}=\"\$OPTARG\";;" >> $script
#	echo -e "\t\t${G_option}) ${G_var_name}=\"\$OPTARG\";;" >> $script
	fi
done
	echo -e "\t\th)" >> $script
	echo -e "\t\t\thelps" >> $script
	echo -e "\t\t\texit 0" >> $script
	echo -e "\t\t\t;;" >> $script
	echo -e "\tesac" >> $script
	echo -e "done" >> $script
	echo -e "" >> $script
	echo -e "" >> $script

j=""
for j in {A..G}
do
	if [[ $(eval echo -e "\$${j}_par") != "" ]]; then
		eval getopthelper.sh -f $script -t $(eval echo -e "\$${j}_op_type") -n $(eval echo -e "\$${j}_var_name") -o $(eval echo -e "\$${j}_option")
	fi
done

