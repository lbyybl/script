#!/usr/bin/env python
# 2018/11/8
# Boyuan_Li
# It's used to deal the file produced by trans_str.R
import sys
import argparse
parser = argparse.ArgumentParser()
#parser.add_argument("filename")
parser.add_argument("-i", "--Input_file", help="The file produced by trans_str.R")
parser.add_argument("-o", "--output", help="Output filename")

args = parser.parse_args()

if args.Input_file is None or args.output is None:
    print "parameter lacked"
    sys.exit(1)


#import os
import pandas as pd
import numpy as np
file = open(args.Input_file)
file.readline()
line_split = []
list_subtract = []

subtract_count = {}
for line in file:
    a=[]
    line_split = line.split('\t')
    if subtract_count.get(line_split[8],1)==1:
        subtract_count[line_split[8]]=[line_split[9]]
    else:
        a=subtract_count[line_split[8]]
        a.append(line_split[9])
        subtract_count[line_split[8]]=a #subtract_count[line_split[8]].append(line_split[9])
GG=[]
GN=[]
NN=[]
for i in list(subtract_count.items()):
    region, class_region = i
    if "GERE-GERE\n" in class_region:
        GG.append(region)
    elif "NONG-GERE\n" in class_region:
        GN.append(region)
    elif "GERE-NONG\n" in class_region:
        GN.append(region)
    else:
        NN.append(region)

df_interaction = pd.DataFrame([GG,GN,NN], index=['GG','GN','NN']).T

df_interaction.to_csv(args.output,index=False,sep="\t")        
#print(df_interaction["GG"])       
    
#for key in     
#print(subtract_count.items())

