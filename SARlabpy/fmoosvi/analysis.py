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
import pdb

def calculateAUC(data_dict, time = 60):

    # Read in data and headers
    
    x_size = data_dict['data'].shape[0]
    y_size = data_dict['data'].shape[1]    
    num_slices = data_dict['header']['method']['PVM_SPackArrNSlices'][0]
    repetition_time = data_dict['header']['method']['PVM_RepetitionTime']*1E-3
    phase_encodes = data_dict['header']['method']['PVM_EncMatrix'][1] # Num phase encoding steps
    time_per_rep = repetition_time * phase_encodes
    reps = data_dict['header']['method']['PVM_NRepetitions']
    auc_reps = int(numpy.round(time / time_per_rep))
        
    # Determine point of injection by averaging one slice in the entire image
    tumour_curve = enhancementCurve(data_dict)
    pylab.plot(tumour_curve[round(num_slices/2),:], '-bx')
    
    inj_point = sar.determineInjectionPoint(data_dict)
    # Now calculate the Normalized Intesity voxel by voxel

    norm_data = numpy.empty([x_size,y_size,reps])
    auc_data = numpy.empty([x_size,y_size,num_slices])
    
        
    for slice in range(num_slices):
        
        #pdb.set_trace()
        baseline = numpy.mean(data_dict['data'][:,:,slice,0:inj_point],axis=2)
        norm_data = (data_dict['data'][:,:,slice,:] / numpy.tile(baseline.reshape(x_size,y_size,1),reps)) -1    

        # Now calculate the actual AUC    
        auc_data[:,:,slice] = numpy.sum(norm_data[:,:,inj_point:inj_point+auc_reps],axis=2)
    
    return auc_data
#
#import pylab
#import numpy
#x=200
#y=400
#s=3
#
#orig_data = numpy.ones([x,y,s,111])
#orig_data[5:15,15:25,:,10:] = 101
#orig_data[5:15,15:25,:,10] = 111
#
#baseline = numpy.mean(orig_data[:,:,0,0:9],axis = 2)
#
#norm_data = numpy.empty([x,y,s,111])
#auc_data = numpy.empty([x,y,s])
#
#norm_data[:,:,0,:] = (orig_data[:,:,0,:] / numpy.tile(baseline.reshape(x,y,1),111)) -1
#auc_data[:,:,0] = numpy.sum(norm_data[:,:,0,10:70],axis = 2)
#
#pylab.imshow(auc_data[:,:,0])



def enhancementCurve(data_dict, mask=False):

    if mask:
        print "I don't know how to deal with masks quite yet"
        
        #TODO: Complete this to read in an ROI somehow
        
    else:
        
        img_mean = data_dict['data'][:,:,:,:].mean(0).mean(0)
            
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

    for slice in range(num_slices):
        
        diff_slice = numpy.diff(img_mean[slice,:])
        std_slice =  numpy.std(diff_slice)
        
        try: 
            injection_point.append(next(i for i,v in enumerate(diff_slice) if v > 3*std_slice))
        except StopIteration:
            print "Could not find the injection point, possibly okay" + str(slice)
            injection_point.append(0)
            
    # look through the list of elements in injection pont and report the most common (mode)
    injection_point_counter = Counter(injection_point)
    injection_point = injection_point_counter.most_common(1)

    return injection_point[0][0]








