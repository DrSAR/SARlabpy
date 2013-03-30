# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 05:23:58 2013

@author: firas
"""
import sarpy
import numpy
import sarpy.fmoosvi.wrappers
import json


## Script to save IR-RARE files as niftis


with open('/Volumes/Data/Dropboxes/PhD./Dropbox/Studies/NecS3/NecS3.json','r') as infile:
    master_sheet = json.load(infile)

NecS3 = sarpy.Experiment('NecS3').find_scan_by_protocol('04')

for k,v in master_sheet.iteritems():
    
    try:

        sarpy.Scan(master_sheet[k]['24h-IR_A'][0]).pdata[0].write2nii(k)
        
    except:
        print k
    
