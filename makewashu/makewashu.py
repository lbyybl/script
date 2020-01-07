#!/home/boyuanli/tools/anaconda2/envs/python35/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 14:46:11 2019

@author: lbybyl

This script is used to add track to washu genome browser,
You should give track type, url, track name;
"""

import json
import os
import copy
import sys
import click
import os.path
import sys
import re
import time

@click.group()
def cli():
    pass

@click.command()
#@click.argument('mode')
@click.option('--url-file', "-u", required=True ,help='The file store the url you want load to washu')
@click.option('--output-file', "-o", required=True ,help='The output file name')
@click.option('--input-file', "-i", default = "haha" , help='The input file name, should be the washu json file, can skip')

def add_washu_track(url_file, output_file, input_file):

    """
    This script is used to add track to washu genome browser,
    You should give track type, url, track name;
    url file can only include the url or name+\t+url
    """
    if os.path.exists(url_file):
        pass
    else:
        sys.exit('ERROR, no file named (%s)' % url_file)

    templete='/home/boyuanli/bashscript/bin/makewashu/washu_track_templete.json'
    url_list=url_file
    with open(templete,'r') as templete:
        templete=json.load(templete)
        
    hic_temp = copy.deepcopy(templete['tracks'][2])
    bigwig_temp = copy.deepcopy(templete['tracks'][5])
    
    
    if input_file == 'haha':
        psudo_input=copy.deepcopy(templete)
        for i in range(2,len(psudo_input['tracks'])):
            del(psudo_input['tracks'][2])
        input_f=psudo_input
    else:
        '''
        Input_file should be a json_file you have download from the ucsc
        '''
        with open(input_file) as f:
            input_f = json.load(f)
        
    '''
    Url_file should have the url of file you want to load
    if there tow col the firs col is name and the second col is url,
    and seperate by tab
    if only one col, the col must be url
    '''
    with open(url_list, 'r') as file:
        url_file = file.read()
    
    def add_track(type,name,url):
        if type=='bigwig':
            track = copy.deepcopy(bigwig_temp)
        elif type=='hic':
            track = copy.deepcopy(hic_temp)
        else:
            print('The track type you input is not support!!!')
            sys.exit(1)
        track['label']=name
        track['options']['label']=name
        track['name']=name
        track['url']=url
        return track
    
    for track_url in url_file.split('\n'):
        if '.bw' in track_url or '.bigwig' in track_url:
            file_type='bigwig'
        elif '.hic' in track_url:
            file_rype='hic'
        else:
            continue
        url_line=track_url.split('\t')
        if len(url_line) >= 2:
            file_name=track_url.split('\t')[0]
            file_url=track_url.split('\t')[1]
        elif len(url_line) == 1:
            print('You input file only have the url and \n the file name will use the default')
            prefile_name=track_url.split('/')[-1]
            file_name=prefile_name.split('.')[0]
            file_url=track_url
        track=add_track(file_type,file_name,file_url)
        input_f['tracks'].insert(2,track)
        
    with open(output_file,'w') as output:
        json.dump(input_f,output)
        
cli.add_command(add_washu_track)
if __name__ == '__main__':
    cli()
