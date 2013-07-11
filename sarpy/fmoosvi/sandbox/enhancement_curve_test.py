# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 22:49:33 2013

@author: fmoosvi
"""
import sarpy
import numpy
import pylab
import copy
import scipy
import sarpy.fmoosvi.analysis

def h_enhancement_curve(scan_object, adata_roi_label, pdata_num = 0):

    try:
        norm_data = sarpy.fmoosvi.analysis.h_normalize_dce(scan_object)
        num_slices = norm_data.shape[-2]
        reps = norm_data.shape[-1]
        
        #Calculating the time per rep.
        phase_encodes = scan_object.pdata[pdata_num].visu_pars.VisuAcqPhaseEncSteps
        repetition_time = scan_object.method.PVM_RepetitionTime*1E-3
        time_per_rep = repetition_time * phase_encodes
        time = numpy.linspace(0,99,num=100)*time_per_rep
                       
        ## THIS IS INCREDIBLY SKETCHY, AND I'M NOT SURE WHAT THE RAMIFICATIONS ARE        
        new_scan_object = copy.deepcopy(scan_object)
        new_scan_object.pdata[pdata_num].data = norm_data
        data_scan = new_scan_object.pdata[pdata_num]
        ## END SKETCHY BIT
        
        roi = scan_object.adata[adata_roi_label].data   
              
        masked_data = data_scan.data * numpy.tile(numpy.reshape(roi,[
roi.shape[0], roi.shape[1], roi.shape[2],1]),reps)

        enhancement_curve = numpy.empty(shape = [num_slices, 2, reps])
        
        for slice in range(num_slices):
            
            enhancement_curve[slice,0,:] = time
            enhancement_curve[slice,1,:] = scipy.stats.nanmean(scipy.stats.nanmean(masked_data[:,:,slice,:], axis=0), axis =0)            
        return enhancement_curve
    except:    
        print("Perhaps you didn't pass in a valid mask or passed bad data")
        raise


a=sarpy.Scan('NecS3Hs13.iM1/14')

ec = h_enhancement_curve(a, 'auc60_roi')