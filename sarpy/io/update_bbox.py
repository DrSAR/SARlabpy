#!/usr/bin/env python 
# -*- coding: utf-8 -*-

"""
Created on Thu May 30 17:49:35 2013

@author: fmoosvi
"""

def update_bbox(masterlist_name, scan_label, adata_label=None, prefix= None):
    
    import sarpy.fmoosvi.analysis
    import sarpy.fmoosvi.getters
    import os
    import json
    import re
    import collections

    if adata_label is None:
        adata_label = 'roi'
    
    if prefix is None:
        prefix = ''


    ## Start changing the master lists, see whether an updated masterlist exists
    
    root = os.path.expanduser(os.path.join('~/sdata',
                        masterlist_name,
                        masterlist_name))
    if os.path.exists(os.path.join(root+'_updated.json')):
        fname_to_open = root+'_updated.json'
    else: 
        fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str)    
        
    for k,v in master_list.iteritems():
        
        try:
            scn = sarpy.Scan(v[scan_label][0])
        except(IOError):
            print('update_bbox: Could not find {0}, {1} \n'.format(k,adata_label))
        
        else:
            try:
                new_bbox = sarpy.fmoosvi.getters.get_roi_bbox(scn,adata_label,
                                                              type='pct')
            except KeyError:
                print('update_bbox: Could not key {0} in scan {1} \n'.format(adata_label,
                      k))
                pass
    
            for j in v:
                
                if len(v[j][1]) == 4 and re.match(prefix,str(j)): 
                    v[j][1] = new_bbox
                    # leave the empty ones as is...        
                
    json.dump(master_list, open(os.path.join(os.path.expanduser('~/sdata'),
                                             masterlist_name,
                                             masterlist_name+'_updated.json'),'w'), 
                                             indent=4)
    
    
    print('update_bbox: Success, updated file is: {0}'.format(
                                                masterlist_name+'_updated.json'))
    
        
if __name__ == '__main__':
    
    import argparse

    #TODO: for some reason the required args still show up as optional
    
    ## dealing with parameters
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-m', '--masterlist_name', type=str, required=True,
                       help='Name of the masterlist file. e.g. NecS3. \
                       Usage: python update_bbox.py -m NecS3 -a roi -i export ')
    
    parser.add_argument('-s', '--roi_source', type=str, required=True,
                       help='Required: source of scan used to draw roi') 
    
    parser.add_argument('-a', '--adata_label', type=str,
                       help='Optional: adata label, defaults to roi') 
                       
    parser.add_argument('-x', '--prefix', type=str,
                       help='Optional: only modifies labels that start with the \
                       prefix, defaults to all labels. Useful for scans done on different days, e.g., 0h- or 24h-')
                       
                     
    args = parser.parse_args()
    masterlist_name = args.masterlist_name
    scan_label = args.roi_source
    adata_label = args.adata_label
    
    try:
        prefix = str(args.prefix)
    except AttributeError:
        prefix = None
        
    
    update_bbox(masterlist_name, scan_label, adata_label, prefix)