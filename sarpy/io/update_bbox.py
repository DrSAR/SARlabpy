#!/usr/bin/env python 
# -*- coding: utf-8 -*-

"""
Created on Thu May 30 17:49:35 2013

@author: fmoosvi
"""

import sarpy.fmoosvi.analysis
import sarpy.fmoosvi.getters
import argparse
import os
import json
import re

#TODO: for some reason the required args still show up as optional

## dealing with parameters
parser = argparse.ArgumentParser()

parser.add_argument('-m', '--masterlist_name', type=str, required=True,
                   help='Name of the masterlist file. e.g. NecS3. \
                   Usage: python update_bbox.py -m NecS3 -a roi -i export ')

parser.add_argument('-s', '--roi_source', type=str, required=True,
                   help='Optional: source of scan used to draw roi') 

parser.add_argument('-a', '--adata_label', type=str,
                   help='Optional: adata label, defaults to roi') 
                   
parser.add_argument('-x', '--suffix', type=str,
                   help='Optional: only modifies labels that start with the \
                   suffix, defaults to all labels. Useful for scans done on different days, e.g., 0h- or 24h-')
                   
                 
args = parser.parse_args()
masterlist_name = args.masterlist_name
scan_label = args.roi_source

try:
    adata_label = args.adata_label
except AttributeError:
    adata_label = 'roi'
    
try:
    suffix = str(args.suffix)
except AttributeError:
    suffix = None    


## Start changing the master lists, see whether an updated masterlis exists

root = os.path.join(os.path.expanduser('~/mdata'),masterlist_name,masterlist_name)

if os.path.exists(os.path.join(root+'_updated.json')):
    
    with open(os.path.join(root+'_updated.json'),'r') as master_file:
        master_list = json.load(master_file).copy()   
else:
    with open(os.path.join(root+'.json'),'r') as master_file:
        master_list = json.load(master_file).copy()
    
    
for k,v in master_list.iteritems():
    
    try:
   
        scn = sarpy.Scan(v[scan_label][0])
        new_bbox = sarpy.fmoosvi.getters.get_roi_bbox(scn,adata_label,type='pct')
        
        for j in v:
            
            if len(v[j][1]) == 4 and re.match(suffix,str(j)): 
                v[j][1] = new_bbox
                # leave the empty ones as is...        
            
    except IOError:
        
        print('update_bbox: Coud not find {0}, {1} \n'.format(k,v))

json.dump(master_list, open(os.path.join(os.path.expanduser('~/mdata'),masterlist_name,masterlist_name+'_updated.json'),'w'), indent=4)




    