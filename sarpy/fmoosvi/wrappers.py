# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:11:03 2013

@author: fmoosvi
"""

import sarpy
import sarpy.fmoosvi.analysis

def calculate_AUC(Bruker_object, protocol_name = '06_FLASH2D.6sl-DCE', pdata_num = 0, time = 60):
    
    scan_list = []
    auc = []
    
    if type(Bruker_object) == sarpy.io.BRUKER_classes.Experiment:
        scan_list = Bruker_object.find_scan_in_experiment(protocol_name)
        print('You put in an experiment')
        
    elif type(Bruker_object) == sarpy.io.BRUKER_classes.Study:
        scan_list = Bruker_object.find_scan(protocol_name)
        print('You put in a single study')
            
    elif type(Bruker_object) == sarpy.io.BRUKER_classes.Scan:
        scan_list = Bruker_object
        print('You put in a single scan')
            
    elif type(Bruker_object) == list:
        
        if type(Bruker_object[0]) == sarpy.io.BRUKER_classes.Study:                                  
            for study in Bruker_object:
                scan_list.append(study.find_scan(protocol_name))
            print('You put in a list of studies')

        elif type(Bruker_object[0]) == sarpy.io.BRUKER_classes.Scan:            
             scan_list = Bruker_object
             print('You put in what you were supposed to!')
            
    for scan in scan_list:
        auc.append(sarpy.fmoosvi.analysis.h_calculate_AUC(scan))
        
    return auc
    
def calculate_BSB1map(Experiment_object, protocol_name = '04_ubcLL2'):
    
    
    
    
    return
    