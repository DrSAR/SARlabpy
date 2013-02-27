# Copyright (C) 2012-2013 Stefan A Reinsberg and SARlab members
# full license details see LICENSE.txt
"""Collection of analysis routines for Bruker data
test

"""

from __future__ import division

import numpy
import os
import pylab
import SARlabpy as sar


def calculateAUC(study_path, series_num = 0):

    # Read in data and headers

    procdirname = os.path.join(study_path,series_num)

    print procdirname
       
    data_dict = sar.read2dseq(procdirname)
        
    repetition_time = data_dict['header']['method']['PVM_RepetitionTime']*1E-3
    num_slices = data_dict['header']['method']['PVM_SPackArrNSlices'][0]
    reps = data_dict['header']['method']['PVM_NRepetitions']
        
    # Determine point of injection by averaging one slice in the entire image
    tumour_curve = enhancementCurve(data_dict)
    pylab.plot(tumour_curve[round(num_slices/2),:], '-bx')
    
    inj_point = sar.determineInjectionPoint(data_dict)

    # Now calculate the Normalized Intesity AUC voxel by voxel

    auc_data = data_dict['data'] # initialize auc_data
    
    x_size = data_dict['data'].shape[0]
    y_size = data_dict['data'].shape[1]
    
    baseline = numean(a[:,:,0:inj_point],axis=2)
    result = a[:,:,inj_point:]/tile(baseline.reshape(x_size,y_size,1),40)    
    
    
    # Time this code block and compare it with the vectorized form to convince
        # yourelf that the above is faster
    
#    for x in range(0,data_dict['data'].shape[0]):
#    
#        for y in range(0, data_dict['data'].shape[1]):
#
#            for z in range(0,data_dict['data'].shape[2]):
#            
#                auc_data[x,y,z,inj_point:] = ( data_dict['data'][x,y,z,inj_point:] / numpy.mean(data_dict['data'][x,y,z,0:inj_point])) -1
#                
    return auc_data


def enhancementCurve(data_dict, mask=False):

    if mask:
        print "I don't know how to deal with masks quite yet"
        
        #TODO: Complete this to read in an ROI somehow
        
    else:
        
        img_mean = data_dict['data'][:,:,:,:].sum(0).sum(0)
            
        return img_mean

def determineInjectionPoint(data_dict):
    
    from collections import Counter

    try:      
        # Pool all the data together
        img_mean = data_dict['data'][:,:,:,:].sum(0).sum(0)

    except IndexError:
        print "You might only have 2D or 3D data, check data source"
        raise IndexError
        
    num_slices = data_dict['header']['method']['PVM_SPackArrNSlices'][0]
    injection_point = []

    for slice in range(0,num_slices):
        
        diff_slice = numpy.diff(img_mean[slice,:])
        std_slice =  numpy.std(diff_slice)
        
        try: 
            injection_point.append(next(i for i,v in enumerate(diff_slice) if v > 3*std_slice))
        except StopIteration:
            print "Could not find the injection point, possibly okay"
            injection_point.append(0)
            
    # look through the list of elements in injection pont and report the most common (mode)
    injection_point_counter = Counter(injection_point)
    injection_point = injection_point_counter.most_common(1)

    return injection_point[0][0]








