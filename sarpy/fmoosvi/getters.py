# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 14:41:54 2013

@author: fmoosvi
"""

import sarpy
import sarpy.ImageProcessing.resample_onto
import nibabel
import scipy
import scipy.stats

import numpy
import copy


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
            
def get_enhancement_curve(scan_object, adata_mask, pdata_num = 0):

    try:
        norm_data = sarpy.fmoosvi.analysis.h_normalize_dce(scan_object)
        num_slices = norm_data.shape[-2]
        reps = norm_data.shape[-1]        
                
        ## THIS IS INCREDIBLY SKETCHY, AND I'M NOT SURE WHAT THE RAMIFICATIONS ARE        
        new_scan_object = copy.deepcopy(scan_object)
        new_scan_object.pdata[pdata_num].data = norm_data
        data_scan = new_scan_object.pdata[pdata_num]
        
        print data_scan
        ## END SKETCHY BIT
        
        roi_image = sarpy.ImageProcessing.resample_onto.resample_onto_pdata(adata_mask,data_scan)   

        roi_mask= sarpy.fmoosvi.analysis.h_image_to_mask(roi_image)
        
        masked_data = data_scan.data * numpy.tile(numpy.reshape(roi_mask,[
roi_mask.shape[0], roi_mask.shape[1], roi_mask.shape[2],1]),reps)


        enhancement_curve = numpy.empty(shape = [num_slices, reps])
        
        for slice in range(num_slices):

            enhancement_curve[slice,:] = scipy.stats.nanmean(scipy.stats.nanmean(masked_data[:,:,slice,:], axis=0), axis =0)

        return enhancement_curve
    except:    
        print("Perhaps you didn't pass in a valid mask or passed bad data")
        raise
        
def get_bbox(value,data_label,type=None):
    
    data = sarpy.Scan(value[data_label][0])
    shape = data.pdata[0].data.shape
    
    bbox = numpy.array([float(x) for x in value[data_label][1]])    
    bbox_px = (bbox.reshape(2,2).T*shape[0:2]).T.flatten()
    
    #TODO: this will need to be updated for python 3.x+
    bbox_px = map(int,bbox_px) # Casts all elements to ints
    
    if type is None:
        return numpy.array(bbox_px)
    
    elif type == 'pct':
        return numpy.array(bbox)
    
    
    
    
    
    
    
    
    
