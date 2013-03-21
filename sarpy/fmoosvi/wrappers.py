# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:11:03 2013

@author: fmoosvi
"""

import numpy
import sarpy
import sarpy.fmoosvi.analysis

def calculate_AUC(Bruker_object, protocol_name = '06_FLASH2D+', pdata_num = 0, time = 60):
       
    try:
        scan_list = Bruker_object.find_scan_by_protocol(protocol_name)
    except:
        try: 
            scan_list = []                                
            for study in Bruker_object:
                scan_list.append(study.find_scan(protocol_name))
            print('You put in a list of studies, try to avoid that')
        except:
            raise

    ## Now to calculate the AUC

    auc = []

    for scan in scan_list:
        auc.append(sarpy.fmoosvi.analysis.h_calculate_AUC(scan))
        
    return auc
    
def calculate_BSB1map(Experiment_object, protocol_name = ''):
    
    # TODo = need to actually fix this!
    
    print ("Why are you using this!? It' not implemented yet")
    
    try:
        scan_list = Experiment_object.find_scan_by_protocol(protocol_name)
    except:
        try: 
            scan_list = []                                
            for study in Experiment_object:
                scan_list.append(study.find_scan(protocol_name))
            print('You put in a list of studies, try to avoid that')
        except:
            raise

    ## Now to calculate the B1 map

    b1map  = []
    for scan in scan_list:
        b1map.append(sarpy.fmoosvi.analysis.h_BS_B1map(scan))
        
    return b1map

#TODO: Implemen BSB1 mapwrapper around h_BS_B1map

######################### Start Example #########################
## This is for one phantom study using LL and BSB1 MSME
#APExp = sarpy.Experiment('DiLL.iI1')
#BS_scans = APExp.find_scan_by_protocol('09\-2DBSB1map\-MSME_bas\(modified\)')[-4:]
#
#LLdata = sarpy.Experiment('DiLL.iI1').find_scan_by_protocol('MOBILE*')[0]
#
#zero_BS_minus = BS_scans[0]
#zero_BS_plus = BS_scans[1]
#power_BS_minus = BS_scans[2]
#power_BS_plus = BS_scans[3]
#
#b1map = sarpy.fmoosvi.analysis.h_BS_B1map(zero_BS_minus,\
#                                          zero_BS_plus,\
#                                          power_BS_minus,\
#                                          power_BS_plus,\
#                                          LLdata)
#                                          
#pylab.imshow(b1map[:,:,0,0], aspect = 0.5)
#pylab.colorbar()
######################### End Example #########################

    return #b1map
    
def calculate_T1map(Bruker_object, protocol_name = '04_ubcLL2', flip_angle_map = 0):


#TODO: Fix this so the list/object situation is sorted!
#    try:
#        
#        if type(Bruker_object) == sarpy.io.BRUKER_classes.Scan:
#            scan_list = Bruker_object
#        else:
#            scan_list = Bruker_object.find_scan_by_protocol(protocol_name)
#    except:
#        try: 
#            scan_list = []                                
#            for study in Bruker_object:
#                scan_list.append(study.find_scan_by_protocol(protocol_name))
#            print('You put in a list of studies, try to avoid that')
#        except:
#            raise

    ## Now to calculate the AUC
    scan_list = []
    scan_list.append(Bruker_object)
    T1_map = []

    for scan in scan_list:
        T1_map.append(sarpy.fmoosvi.analysis.h_fit_T1_LL(Bruker_object,\
                                                      flip_angle_map = flip_angle_map))
    return T1_map
    
    
    
    
    
    
    
    
    
    
    
    