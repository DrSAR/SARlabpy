#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
Created on Wed Nov. 13

@author: fmoosvi
"""
def special_finder(masterlist_name, 
                   data_label, 
                   adata_label,
                   excludePatients = None):
    
    import json
    import os
    import collections
    import csv
    import scipy.stats
    import sarpy
    import time
    import getpass # used to get the current username
    import numpy
    import sarpy.fmoosvi.getters

    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)
    fname_to_open = root+'.json'
    
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str) 

    scan_list = []

    for k,v in master_list.items():

        if k in excludePatients:
            print(('Excluded {0}'.format(k)))
            continue
        try:        
            data = sarpy.Scan(v[data_label][0]).adata[adata_label].data
            scan_list.append(v[data_label][0])

        except(IOError,KeyError):
            continue 

    return scan_list            

def getCumData(day0scans, adata_label,roi_label):

    import sarpy
    import numpy
    
    dataa = []
    
    if day0scans:

        for t in day0scans:
        
            d = sarpy.Scan(t).adata[adata_label].data
            roi = sarpy.Scan(t).adata[roi_label].data
        
            droi = d*roi

            roi[numpy.isfinite(roi)] = True
            roi[numpy.isnan(roi)] = False
            droi = droi[roi.astype(int)==1]
            #droi = d[numpy.isfinite(droi)]
                      
        dataa.append(list(droi.flatten()))
    
    return [x for x in dataa]    
