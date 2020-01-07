#!/home/boyuanli/tools/anaconda2/envs/python35/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 16:34:27 2019

@author: lbybyl

 this file is used to find most max 50 bp window; so the stastic to find
 sig promoter and stastic the paused sig for the 50bp windows, you should
 used R, becaused R more convenient
 input: the sig promoter has been select by R 
"""
import click
import os.path
import sys
import os
import pandas as pd
import numpy as np

@click.command()
@click.option('--inputfile', "-i", required=True , help='the input file name you want find max 50bp window')
@click.option('--outputprefix', "-p", required=True ,help='output prefix')
#@click.option('--first-par', "-a", default = "" , help='should be Var_name,Short_name,is_requrie,type [optional]')

def main(inputfile, outputprefix):
    """
    inputfile example:
    chr1    4806822 4808822 NM_008866       0       +       chr1    4806822 4806847 26      0       +       25
	it was got by the bedtools intersect for the promoter region and reads;
    it used to find the max 50kp window and it's reads num for file like up;
    """
    filname=inputfile
    outputname=outputprefix
    file = open(filname,'r')
    
    gene_dir = {}
    gene_info = {}
    line_split = []
        
    for line in file:
        line_split = line.split('\t')
        if line_split[5]=="+":
            try:
                gene_dir[line_split[3]].append(line_split[7])
            except KeyError:
                gene_info[line_split[3]] = line_split[0:6]
                gene_dir[line_split[3]] = [line_split[7]]
        elif line_split[5]=="-":
            try:
                gene_dir[line_split[3]].append(line_split[8])
    
            except KeyError:
                gene_info[line_split[3]] = line_split[0:6]
                gene_dir[line_split[3]] = [line_split[8]]
                
    file.close()            
    """        
    plus 400 num into the list coresponding to each gene
    and the find the most sig 50 bp for the active promoter
     
    now you have two dict; one contain every gene's mapping reads loc;
    another contain every gene's common info,such as chr,st,en,id,0,strand
    
    process: for each gene, append 400 elements into the list,and sort it, got
    the index of the elements you append,  and substract them with the window/strp as
    step, and got the max value for substract list;
    and got the index, the st_region=gene_st*(index+1) en_region=gene_st*(index+1)+50
    
    output: format the 50 bp region data to chr,st,en,id,0,strand,reads_num
    """

    
    sp=5 # sliding steps  
    windwos=10 # 10*sp  
    gene_window = {}
    gene_window_readsnum = {}
    gene_max_window = {}
    for i in gene_dir.keys():
        try:
            ex_list=produce_list(gene_info[i][1],gene_info[i][2],sp)
        except  IndexError:
            print('i=',i)
            print(gene_info[i])
            print(gene_info[i])
            sys.exit(3)
        try:
            gene_dir[i] = gene_dir[i] + ex_list
            gene_dir[i]=list(map(int,gene_dir[i]))
            gene_dir[i].sort()
        except AttributeError:
            print(gene_dir[i])
            sys.exit(4)
        except TypeError:
            print(type(gene_dir[i]))
            print(gene_dir[i])
            sys.exit(4.1)
        try:
            gene_window[i]=got_eleindex_for_list(gene_dir[i],ex_list)
    
        except TypeError:
            print('i=',i)
            print(gene_dir[i])
            sys.exit(5)
        try:
            gene_window_readsnum[i]=got_list_diff(gene_window[i],windwos)
        except TypeError:
            print(gene_window[i])
            print(type(gene_window[i]))
            sys.exit(7)
        gene_max_window[i]=[max(gene_window_readsnum[i])]
        max_idx=gene_window_readsnum[i].index(max(gene_window_readsnum[i]))
        wi_st=int(gene_info[i][1])+(max_idx+1)*sp
        tmp_list=[gene_info[i][0],wi_st,wi_st+50,gene_info[i][3],gene_info[i][4],gene_info[i][5]]
        gene_max_window[i]=gene_max_window[i]+tmp_list
        


    save_file(gene_max_window,outputname)

def produce_list(st,en,sp):
    st=int(st)
    en=int(en)
    li = [x for x in range(st,en+sp,sp)]
    return li

def got_eleindex_for_list(li_all,li_sub):
    """
    input a list and a sub list and return the elements index of 
    the sulist in the list
    """
    li_index=[]
    li_all=list(map(int,li_all))
    for i in li_sub:
        li_index.append(li_all.index(i))
    return li_index

def got_list_diff(li,sp):
    li_diff=[]
    len_li=len(li)
    for i in range(0,len_li):
        if i+sp < len_li:
            li_diff.append(li[i+sp]-li[i]-sp)
    return li_diff  

def save_file(dict_data,outputname):
    dict_data_df = pd.DataFrame(dict_data)
    dict_data_df = dict_data_df.T
    dict_data_df.to_csv(outputname,index=False,sep="\t")  
    
if __name__ == '__main__':
    main()
