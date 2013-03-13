# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 14:41:54 2013

@author: fmoosvi
"""

def get_num_slices(scan_object, pdata_num = 0):
    
    """
    Returns the number of slices

    :param object scan_object: scan object from a study
    :param integer pdata_num: reconstruction number, according to python numbering.
            default reconstruction is pdata_num = 0.
    :return: integer num_slices: number of slices, detrmined from visu_pars
    """     
    
       ## Ridiculouly long (but definitely complete and correct) method of 
       #  obtaining the third dimension from visu_pars
    
    if scan_object.pdata[pdata_num].visu_pars.VisuCoreDim == 2:
        num_slices = 1 # we hope this will get updated below in the case of
              # multi-slice 2D data
    else:
        num_slices = scan_object.pdata[pdata_num].visu_pars.VisuCoreSize[2]
    try:
        VisuFGOrderDesc = scan_object.pdata[pdata_num].visu_pars.VisuFGOrderDesc
    except KeyError:
        # data missing, we should have all we need anyway
        pass
    else:
        for VisuFGOrderDesc_element in VisuFGOrderDesc:
            if VisuFGOrderDesc_element[1].strip()=='<FG_SLICE>':
                num_slices = int(VisuFGOrderDesc_element[0])
                
    return num_slices
    
def get_matching_scans(study_object, desired_protocol):

    """
    Returns array of scans that match the protocol name(s)

    :param object Experiment: Object with studies from a study, 
            e.g., say NecS3Hs01.hi1
    :param string protocol_name: Name of protocol being searched, 
            e.g., '06_FLASH2D.6sl-DCE'
    :return: integer num_slices: number of slices, detrmined from visu_pars
    """
    
    matching_scans = []    
    
    for scan in study_object.scans:
        
        protocol = scan.acqp.ACQ_protocol_name
        
        if protocol == desired_protocol:
            
    
    
    acqp.ACQ_protocol_name
    
    
