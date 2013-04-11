# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 14:41:54 2013

@author: fmoosvi
"""

import sarpy
import sarpy.ImageProcessing.resample_onto
import nibabel 

import numpy

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


def get_patients_from_experiment(Experiment_name, verbose = False):
    
    subject_list = get_unique_list_elements(sarpy.Experiment(Experiment_name).get_SUBJECT_id())

    Patients = []  
    
    for subject in subject_list:
        curr_patient = sarpy.Patient(subject)
        Patients.append(curr_patient)        
        
        if verbose:
            print('Patient {0} has {1} sessions.'.format(subject,len(curr_patient.studies)))
    
    return Patients
    
def get_unique_list_elements(list, idfun=None):

    #Get uniue items in list (maintaining order)
    # From: http://www.peterbe.com/plog/uniqifiers-benchmark
    
   # order preserving
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in list:
       marker = idfun(item)
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result         

def get_image_clims(data):
    
    min_lim = numpy.median(data)
    max_lim = data.max()
    
    return [min_lim, max_lim]
    
def get_goodness_map(data, fit_dict):
    
    goodness_map = data[:]*numpy.zeros()
    
    # TODO: gahhh can't believe I have to use nested loops here. Fix it!
    
    for x in xrange(data.shape[0]):
        
        for y in xrange(data.shape[1]):
            
            for z in xrange(data.shape[2]):
                
                goodness_map[x,y,z] = fit_dict['goodness']
    
    return goodness_map
            

#def get_resampled_data(roi_scan_object, resampled_scan_object, adata_key = 'IR_tumour_rois', pdata = 0):
#
#
#    roi_data = roi_scan_object.adata[adata_key].data
#    roi_scan_object.adata[adata_key].export2nii('res_data.nii')
#    
#    fixed_data = resampled_scan_object.pdata[pdata].data
#    resampled_scan_object.pdata[pdata].export2nii('fixed_data.nii')
#    
#    resampled_data = sarpy.ImageProcessing.resample_onto.resample_onto('res_data.nii','fixed_data.nii')
#    
#    
#    if fixed_data.shape == roi_data.shape:
#        
#        # continue
#        
#    else:
#    
##        sarpy.Scan(master_sheet[k][key_list[2]][0]).adata[key_list[6]].export2nii('LL1.nii')
##        sarpy.Scan(master_sheet[k][key_list[3]][0]).adata[key_list[6]].export2nii('LL2.nii')
##        
##        roi1 = sarpy.Scan(master_sheet[k][key_list[4]][0]).adata[key_list[7]]
##        roi1.export2nii('roi1.nii')
##        
##        roi2 = sarpy.Scan(master_sheet[k][key_list[5]][0]).adata[key_list[7]]
##        roi2.export2nii('roi2.nii')
##      
##        LLresample1 = sarpy.ImageProcessing.resample_onto.resample_onto('LL1.nii','roi1.nii')          
#    
#
#
#
#
#    
#    resampled_data = sarpy.ImageProcessing.resample_onto.resample_onto('res_data.nii','fixed_data.nii')
#    
#    
#    
#    
#    return resampled_data