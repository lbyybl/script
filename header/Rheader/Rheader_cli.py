# -*- coding: utf-8 -*-
import click
#import os
import os.path
import sys
import re
import time

@click.command()
#@click.argument('mode')
@click.option('--rscript', "-s", required=True ,help='R script name you want write [required]')
@click.option('--author', "-z", default = "Boyuan_Li" ,help='your name [optional]')
@click.option('--first-par', "-a", default = "" , help='should be Var_name,Short_name,is_requrie,type [optional]')
@click.option('--second-par', "-b",  default = "" ,help='should be Var_name,Short_name,is_requrie,type [optional]')
@click.option('--third-par', "-c",  default = "" ,help='should be Var_name,Short_name,is_requrie,type [optional]')
@click.option('--fourth-par', "-d",  default = "" ,help='should be Var_name,Short_name,is_requrie,type [optional]')
@click.option('--fifth-par', "-e",  default = "" ,help='should be Var_name,Short_name,is_requrie,type [optional]')
@click.option('--sixth-par', "-f",  default = "" ,help='should be Var_name,Short_name,is_requrie,type [optional]')
@click.option('--seventh-par', "-g",  default = "" ,help='should be Var_name,Short_name,is_requrie,type [optional]')

def main(rscript, author, first_par, second_par, third_par,fourth_par,fifth_par, sixth_par, seventh_par):
	
    """
    This script is used to make the header for your R script;
    and every parameters is optional, and every optional should
    be the format Var_name,Short_name,is_requrie,type, such as
    Input_file,i,1,character; Input_file is the long parameters
    and var name will be used in the R script; i is the short 
    parameters; 1 means the parameters is requred or optional:
    0:no par, 1:required, 2: optional
    type is the parameters type, should be logical；integer；double；
    complex；character；numeric
	
    """
    if os.path.exists(rscript):
        sys.exit("ERROR: Outputfile (%s) already exists; remove it or specify a new file." % Rscript)
        
    
    script = open(rscript, 'w')
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

