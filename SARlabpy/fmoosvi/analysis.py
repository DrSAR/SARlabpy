# Copyright (C) 2012-2013 Stefan A Reinsberg and SARlab members
# full license details see LICENSE.txt
"""Collection of analysis routines for Bruker data


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
    whole_tumour_curve = enhancementCurve(data_dict)
    pylab.plot(whole_tumour_curve[round(num_slices/2),:], '-bx')

    #TODO Complete this to actually get the point of injection

    # Now calculate the AUC voxel by voxel


def enhancementCurve(data_dict, mask=False):

    if mask:
        print "I don't know how to deal with masks quite yet"
    else:
        
        img_mean = data_dict['data'][:,:,:,:].sum(1).sum(1)
            
        return img_mean

def determineInjectionPoint(data_dict):

    try:      
        # Pool all the data together
        img_mean = data_dict['data'][:,:,:,:].sum(1).sum(1)

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

    return injection_point, diff_slice








