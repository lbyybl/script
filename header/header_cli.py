# -*- coding: utf-8 -*-
import click
#import os
import os.path
import sys
import re
import time

@click.group()
def cli():
    pass

@click.command()
#@click.argument('mode')
@click.option('--yourscript', "-s", required=True ,help='R script name you want write')
@click.option('--author', "-z", default = "Boyuan_Li" ,help='your name;defult:Boyuan_Li [optional]')
@click.option('--first-par', "-a", default = "" , help='should be Var_name,Short_name,is_requrie,type [optional]')
@click.option('--second-par', "-b",  default = "" ,help='should be Var_name,Short_name,is_requrie,type [optional]')
@click.option('--third-par', "-c",  default = "" ,help='should be Var_name,Short_name,is_requrie,type [optional]')
@click.option('--fourth-par', "-d",  default = "" ,help='should be Var_name,Short_name,is_requrie,type [optional]')
@click.option('--fifth-par', "-e",  default = "" ,help='should be Var_name,Short_name,is_requrie,type [optional]')
@click.option('--sixth-par', "-f",  default = "" ,help='should be Var_name,Short_name,is_requrie,type [optional]')
@click.option('--seventh-par', "-g",  default = "" ,help='should be Var_name,Short_name,is_requrie,type [optional]')

def r(yourscript, author, first_par, second_par, third_par,fourth_par,fifth_par, sixth_par, seventh_par):
	
    """
    It's used to make header for R script;
    and every parameters is optional, and every optional should
    be the format Var_name,Short_name,is_requrie,type, such as
    Input_file,i,1,character; Input_file is the long parameters
    and var name will be used in the R script; i is the short 
    parameters; 1 means the parameters is requred or optional:
    0:no par, 1:required, 2: optional
    type is the parameters type, should be logical；integer；double；
    complex；character；numeric
	
    """
    if os.path.exists(yourscript):
        sys.exit("ERROR: Outputfile (%s) already exists; remove it or specify a new file." % yourscript)
        
    
    script = open(yourscript, 'w')
    click.echo("#!/usr/bin/env R",script)
    click.echo("# " + time.asctime(), script)
    click.echo("# " + author, script)
    click.echo("\n",script)           
    click.echo("command=matrix(c(",script)
    par_list=[]
    par_arg=[]
    def get_par(b,li,arg):
        if (b):
            par=b.split(",")
            if len(par):
                arg.append(par[0])
                for i in par:
                    li.append(i)
        return(li,arg) 
		
    for i in first_par, second_par, third_par,fourth_par,fifth_par, sixth_par, seventh_par:
        par_list,par_arg=get_par(i, par_list,par_arg)
        
    n=0   
    print('\t',end='',file=script)
    par_list.reverse()
    for c in range(20):
        i=par_list.pop()
        if par_list:            
            n=1+n
            if n % 4 != 0:
                print('"',end='',file=script)
                print(i,end='',file=script)
                print('", ',end='',file=script)
            else:
                print('"',end='',file=script)
                print(i,end='',file=script)            
                print('", ',end='\n',file=script)
                print('\t',end='',file=script)
        else:
            print('"',end='',file=script)
            print(i,end='',file=script) 
            print('" ',end='\n',file=script)
            print('\t',end='',file=script)
            break
			
    click.echo("),byrow=T,ncol=4)",script) 
    click.echo("\n",script)     
    click.echo("args=getopt::getopt(command)",script)
    click.echo("\n",script)
    print("if (",end='',file=script)
    #print("\t",end='',file=script)
    final2=''
    for i in par_arg:
        media=re.sub(r'(^)','is.null(args$',i)    
        final=re.sub(r'($)',') || ',media)
        final2=final2+final
    final3=re.sub(r'(\ \|\|\ $)','',final2)
    print(final3,end='',file=script)
    #print("\n",end='',file=script)	
    click.echo(") {",script)
    click.echo('\tcat(paste(getopt::getopt(command, usage = T), "\\n"))',script)
    #click.echo("\n",script)
    click.echo("\tq()",script)
    #click.echo("\n",script)
    click.echo("}",script)

	
    script.close()    
		
@click.command()
@click.option('--yourscript', "-s", required=True ,help='bash script name you want write')
@click.option('--author', "-z", default = "Boyuan_Li" ,help='your name;defult:Boyuan_Li [optional]')
@click.option('--first-par', "-a", default = "" , help='should be Var_name,Short_name,type [optional]')
@click.option('--second-par', "-b",  default = "" ,help='should be Var_name,Short_name,type [optional]')
@click.option('--third-par', "-c",  default = "" ,help='should be Var_name,Short_name,type [optional]')
@click.option('--fourth-par', "-d",  default = "" ,help='should be Var_name,Short_name,type [optional]')
@click.option('--fifth-par', "-e",  default = "" ,help='should be Var_name,Short_name,type [optional]')
@click.option('--sixth-par', "-f",  default = "" ,help='should be Var_name,Short_name,type [optional]')
@click.option('--seventh-par', "-g",  default = "" ,help='should be Var_name,Short_name,type [optional]')

def bash(yourscript, author, first_par, second_par, third_par,fourth_par,fifth_par, sixth_par, seventh_par):
	
    """
    It's used to make header for bash script;
    and every parameters is optional, and every optional should
    be the format Var_name,Short_name,type, such as
    Input_file,i,string; Input_file is the long parameters
    and var name will be used in the bash script; i is the short 
    parameters; type is the parameters type, should be dir string file；
    	
    """
    if os.path.exists(yourscript):
        sys.exit("ERROR: Outputfile (%s) already exists; remove it or specify a new file." % yourscript)
        
    
    script = open(yourscript, 'w')
    echo_header(script,author)
    short_list=[]
    long_list=[]
    type_list=[]
    for i in first_par, second_par, third_par,fourth_par,fifth_par, sixth_par, seventh_par:
        long_list,short_list,type_list=get_par(i, long_list,short_list,type_list)
    
    help_long=''
    help_function='' 
    getopt_short=''
    case_long=''
    for i in range(len(short_list)):
        if short_list[i]:
            help_long=help_long+"-"+short_list[i]+' ' + '<' + long_list[i]+'> '
            help_function=help_function+"\t"+'echo " -' + short_list[i] + ' ' + type_list[i] + "\t" + "[required] \""+"\n\t" + 'echo \"\"'+"\n"
            getopt_short=getopt_short+short_list[i]+':'
            case_long=case_long+"\t\t"+short_list[i]+') '+long_list[i]+'="$OPTARG";;'+"\n"
        else:            
            break
        
    getopt_short="while getopts \""+getopt_short+"\" optionName\ndo"+"\n\tcase $optionName in"
    case_long=case_long+"\t\t"+'h)'+"\n"+"\t\t\thelps\n"+"\t\t\texit 0\n"+"\t\t\t;;\n"+"\tesac\ndone\n\n"    
    print(help_long+'"',file=script)
    print("\techo \"\"",file=script)    
    print(help_function,end='',file=script)
    print("\t"+"echo \" -h\t\thelp\"",file=script)
    print('\techo ""',file=script)
    print('\techo ""',file=script)
    print('}',file=script)
    print("\n",end='',file=script)
    print("\n",end='',file=script)
    print('if [ $# -eq 0 ]; then',file=script)
    print('\thelps',file=script)
    print('\texit 0',file=script)
    print('fi',file=script)    
    print(getopt_short,file=script)
    print(case_long,file=script)
    judge=''
    for i in range(len(type_list)):
        if type_list[i]=="dir":
            dir_ju=dir_judge(long_list[i])
            judge=judge+dir_ju
        elif type_list[i]=="file":
            file_ju=file_judge(long_list[i])
            judge=judge+file_ju
        elif type_list[i]=="string":
            string_ju=string_judge(long_list[i])
            judge=judge+string_ju
            
    print(judge,end='',file=script)        
            
    script.close()
    
def echo_header(mscript,name='Boyuan_Li'):
    print('#!/bin/bash',file=mscript)
    print("\n",end='',file=mscript)
    print("# "+name+"\n",file=mscript)
    print('# '+time.asctime(),file=mscript)
    print("\n",end='',file=mscript)
    print('set -euo pipefail',file=mscript)
    print('function helps',file=mscript)
    print('{',file=mscript)
    print('\techo ""',file=mscript)
    print('\techo -e "Usage: $0 [options] ',end='',file=mscript)
#    print('{',file=mscript)

def get_par(b,list_l,list_s,list_type):
    if b:
        par=b.split(',')
        if len(par):
            list_l.append(par[0])
            list_s.append(par[1])
            list_type.append(par[2])
    return(list_l,list_s,list_type)

def file_judge(var):
    echo_file="if [[ $"+var+" = \"\" ]]; then\n\t"+"echo \"the $"+var+" file is needed \"\n\texit 1\nelif [[ ! -f $"+var+" ]]; then\n\t echo \"$"+var+":   is not found\"\n\texit 2\nfi\n\n"
    return(echo_file)
def string_judge(var):
    echo_string="if [[ $"+var+" = \"\" ]]; then\n\t"+"echo \" the $"+var+" string is needed \"\n\texit 1\nfi\n\n"
    return(echo_string)
def dir_judge(var):
    echo_dir="if [[ $"+var+" = \"\" ]]; then\n\t"+"echo \"the $"+var+" file is needed \"\n\texit 1\nelif [[ ! -d $"+var+" ]]; then\n\t echo \"$"+var+":   is not found\"\n\texit 2\nfi\n\n"
    return(echo_dir)
    
    
cli.add_command(r)
cli.add_command(bash)    
if __name__ == '__main__':
    cli()