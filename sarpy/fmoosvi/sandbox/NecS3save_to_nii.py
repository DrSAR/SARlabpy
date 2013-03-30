# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 05:23:58 2013

@author: firas
"""
import sarpy
import numpy
import sarpy.fmoosvi.wrappers


## Script to save IR-RARE files as niftis


with open('/Volumes/Data/Dropboxes/PhD./Dropbox/Studies/NecS3/NecS3.json','r') as infile:
    master_sheet = json.load(infile)


for patient in NecS3_patients:
    
    IR = [study.find_scan_by_protocol('05')[0] for study in patient.studies if len(study.find_scan_by_protocol('05'))>0]
    
    
#    if patient.patient_id in exceptions: continue
#        
#    try:     