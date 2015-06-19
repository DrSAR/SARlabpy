# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 14:41:54 2013

@author: fmoosvi
"""

import sarpy
import nibabel
import scipy
import scipy.stats
import os
import numpy
import copy
import pylab

def get_roi_weights(roi):

    # Calculate the weights for each slice
    weights = numpy.zeros(shape=roi.shape[-1])

    roi[numpy.isfinite(roi)] = 1
    ROIpx = numpy.nansum(roi.flatten())

    # Single slice ROIs
    if len(roi.shape)==2:
        weights= [1]

    # Multi slice ROIs
    else:
        for slice in xrange(roi.shape[-1]):    
            weights[slice] = numpy.divide(numpy.nansum(roi[:,:,slice].flatten()),ROIpx)
            weights[numpy.isnan(weights)] = 0

    return weights

def get_tumour_volume(scan_object_name, roi_label):

    """ 
    Returns the area of a volume in units of mm^3

    : scan_object_name: scanname containing the roi mask
    : roi_label: adata label of the roi mask
    """

    scan_object = sarpy.Scan(scan_object_name)
    roi = scan_object.adata[roi_label].data

    px_count = scipy.nansum(roi)

    # Get the x and y dimensions
    px_size = scan_object.method.PVM_SpatResol

    if len(px_size) == 2:
        assert(px_size[0] == px_size[1])

        # Get the z-dimension
        px_size.append(scan_object.method.PVM_SliceThick)

    return px_size[0]*px_size[1]*px_size[2]*px_count

def get_num_slices(scan_object_name, pdata_num = 0):
    
    """
    Returns the number of slices

    :param object scan_object: scan object from a study
    :param integer pdata_num: reconstruction number, according to python numbering.
            default reconstruction is pdata_num = 0.
    :return: integer num_slices: number of slices, detrmined from visu_pars
    """     

    scan_object = sarpy.Scan(scan_object_name)

    
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


def get_patients_from_experiment(Experiment_name, verbose = False,namesOnly=False):
    
    subject_list = get_unique_list_elements(sarpy.Experiment(Experiment_name).get_SUBJECT_id())

    Patients = []  
    
    for subject in subject_list:
        curr_patient = sarpy.Patient(subject)
        Patients.append(curr_patient)        
        
        if verbose:
            print('Patient {0} has {1} sessions.'.format(subject,len(curr_patient.studies)))

    if namesOnly:
        return [a.patient_id for a in Patients]

    else:
        return list(set(Patients))
    
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

def get_image_clims(data,std_modifier=None,allowNegative=False):
    
    xd = data[numpy.isfinite(data)].flatten()
    
    mean = numpy.mean(xd)
    std = numpy.std(xd)
    
    if std_modifier is None:
        std_modifier = 2.5
               
    min_lim = mean - std_modifier*std
    max_lim = mean + std_modifier*std
    
    if min_lim < 0 and allowNegative == False:
        min_lim = 0
    
    return [min_lim, max_lim]
    
def get_goodness_map(data, fit_dict):
    
    goodness_map = data[:]*numpy.zeros()
    
    # TODO: gahhh can't believe I have to use nested loops here. Fix it!
    
    for x in xrange(data.shape[0]):
        
        for y in xrange(data.shape[1]):
            
            for z in xrange(data.shape[2]):
                
                goodness_map[x,y,z] = fit_dict['goodness']
    
    return goodness_map

def get_fid_enhancement(scan_string):
    
    try:
        scan = sarpy.Scan(scan_string)
    except:
        raise
        
    fiddir = os.path.join(scan.dirname,'fid')
    
    data = numpy.abs(sarpy.io.BRUKERIO.readfid(fiddir,resetNR=True)['data'])
    
    x = data.shape[0]
    y = data.shape[1]
    
    xmin = numpy.round(x*0.3)
    xmax = numpy.round(x*0.8)
    
    ymin = numpy.round(y*0.3)
    ymax = numpy.round(y*0.8)
    

    mean = numpy.mean(numpy.mean(data[xmin:xmax,ymin:ymax,3,:],0),0)
    return pylab.plot(mean)

def get_enhancement_curve(scan_object_name, adata_mask=None, pdata_num = 0):

    scan_object = sarpy.Scan(scan_object_name)

    import sarpy.ImageProcessing.resample_onto

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

        if adata_mask is not None:
            
            roi_mask= sarpy.fmoosvi.analysis.h_image_to_mask(roi_image, 
                                                             background=None, 
                                                             foreground=None)
        
            masked_data = data_scan.data * numpy.tile(numpy.reshape(roi_mask,[
roi_mask.shape[0], roi_mask.shape[1], roi_mask.shape[2],1]),reps)

        else:
            masked_data = data_scan.data


        enhancement_curve = numpy.empty(shape = [num_slices, reps])
        
        for slice in range(num_slices):

            enhancement_curve[slice,:] = scipy.stats.nanmean(scipy.stats.nanmean(masked_data[:,:,slice,:], axis=0), axis =0)

        return enhancement_curve
    except:    
        print("Perhaps you didn't pass in a valid mask or passed bad data")
        raise
        
def convert_bbox(scan_object_name, bbox, exporttype=None):
    
    data = sarpy.Scan(scan_object_name)

    if data.pdata[0].data.shape == data.fftfid.shape:
        shape = data.pdata[0].data.shape
    else: # Addresses #281
        shape = data.fftfid.shape

    bbox_px = (numpy.array(bbox).reshape(2,2).T*shape[0:2]).T.flatten()
    
    #TODO: this will need to be updated for python 3.x+
    bbox_px = map(int,bbox_px) # Casts all elements to ints
    
    if exporttype is None:
        return numpy.array(bbox_px)
    
    elif exporttype == 'pct':
        return numpy.array(bbox)
    
def get_roi_bbox(scan_object_name, roi_adata_label = 'roi',exporttype=None):
    
    scan_object = sarpy.Scan(scan_object_name)
    roi = scan_object.adata[roi_adata_label].data  
    
    # First check to see if the roi is filled with 0s and 1s or nans and 1s,
    # otherwise assume it's in 1s and 0s format already

    if numpy.isnan(numpy.sum(roi.flatten())):
        mask = numpy.where(numpy.isnan(roi),0,1) # where ever roi is nan, give it a value of 0, else 1
    else:
        mask = roi

    # Add a third spatial dimension if it's missing.
    if numpy.size(roi.shape) == 2:
        # add an empty dimension to make it 4D, this code appends the exta axis
        mask_sum = mask.copy()

    else: 
        # Next, sum up the array in the slice dimension and convert the aray to a float
        mask_sum = numpy.array(numpy.sum(mask,axis=-1))
        mask_sum = mask_sum.astype(float)
        
        mask_sum[mask_sum != 0 ] = 1 # undo the effects of summing
    
    # Perform calculations to get the bounds of the box in pixels for each dim
    
    a = numpy.sum(mask_sum,axis=0) 
    numpy.where(a>1E-4,a,numpy.nan)
    a_low = numpy.where(a>0)[0][0]
    a_high = numpy.where(a>0)[0][-1]
 
    b = numpy.sum(mask_sum,axis=1)
    numpy.where(b>1E-4,b,numpy.nan)

    b_low = numpy.where(b>0)[0][0]
    b_high = numpy.where(b>0)[0][-1]
    
    # Add a border around this to avoid cliping (aesthetics mostly)
    b_low=b_low -numpy.round(roi.shape[0]*0.10)
    b_high=b_high+numpy.round(roi.shape[0]*0.10)
    a_low=a_low-numpy.round(roi.shape[0]*0.10)
    a_high=a_high+numpy.round(roi.shape[0]*0.10)
    
    # Now perform some checks to make sure the value are within the orig. data
    
    if b_low < 0:
        b_low = 0
    if a_low < 0:
        a_low= 0
    if b_high > roi.shape[0]:
        b_high = roi.shape[0]
    if a_high > roi.shape[1]:
        a_high = roi.shape[1]
    
   
    bbox = [b_low, b_high, a_low, a_high]    

    # Add code segment to convert the bbox into percent as well as pixels
        
    bbox_arr =numpy.array(bbox)

    if exporttype is None:
        return bbox_arr.astype(int)
    elif exporttype == 'pct':
        shape = roi.shape
        bbox[0] = numpy.true_divide(bbox[0],shape[0])
        bbox[1] = numpy.true_divide(bbox[1],shape[0])
        bbox[2] = numpy.true_divide(bbox[2],shape[1])
        bbox[3] = numpy.true_divide(bbox[3],shape[1])
        return bbox
        
def get_rescaled_data(data,minVal=None,maxVal=None,std_modifier=None):
    
#    Max/min formula came from:
#    http://stackoverflow.com/questions/5294955/how-to-scale-down-a-range-of-numbers-with-a-known-min-and-max-value?rq=1

    if minVal is None and maxVal is None:
        [minVal,maxVal] = get_image_clims(data,std_modifier)
    
    rescaled_data =  numpy.true_divide((maxVal-minVal)*(data-minVal),
                                       (maxVal-minVal)) + minVal
    return rescaled_data