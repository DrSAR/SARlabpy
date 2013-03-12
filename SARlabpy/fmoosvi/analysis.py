# -*- coding: utf-8 -*-
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
import scipy.integrate


import SARlabpy.io.BRUKER_classes as sar
readfidExp = sar.Experiment('readfid')
for study in readfidExp.studies:
    print('-'*40+'\n'+study.subject.SUBJECT_id)
    for scan in study.scans:
        print("  "+scan.acqp.ACQ_protocol_name)

def describe_object(class_object):

    try:
        for study in class_object.studies:
            print('-'*40+'\n'+study.subject.SUBJECT_id)    
    except:
        try:
            for scan in class_object.scans:
                print("  "+scan.acqp.ACQ_protocol_name)
        except:
            print("I don't know what's going on")


def calculate_AUC(data_dict, time = 60):
    # Read in data and headers    
    x_size = data_dict['data'].shape[0]
    y_size = data_dict['data'].shape[1]    
    num_slices = data_dict['header']['method']['PVM_SPackArrNSlices'][0]
    repetition_time = data_dict['header']['method']['PVM_RepetitionTime']*1E-3
    phase_encodes = data_dict['header']['method']['PVM_EncMatrix'][1] # Num phase encoding steps
    reps = data_dict['header']['method']['PVM_NRepetitions']
    time_per_rep = repetition_time * phase_encodes
    auc_reps = int(numpy.round(time / time_per_rep))
    time_points = numpy.arange(time_per_rep,time_per_rep*auc_reps + time_per_rep,time_per_rep)
    
    # Determine point of injection by averaging one slice in the entire image
    inj_point = sar.analysis.inj_point(data_dict)
    
    # Now calculate the Normalized Intesity voxel by voxel
    norm_data = normalize_dce(data_dict)

    # Now calculate the actual AUC
    auc_data = numpy.empty([x_size,y_size,num_slices])
    
    for slice in range(num_slices):
        auc_data[:,:,slice] = scipy.integrate.simps(norm_data[:,:,slice,inj_point:inj_point+auc_reps],x=time_points)
    return auc_data

def normalize_dce(data_dict):

    x_size = data_dict['data'].shape[0]
    y_size = data_dict['data'].shape[1]
    inj_point = sar.analysis.inj_point(data_dict)
    reps = data_dict['header']['method']['PVM_NRepetitions']
    num_slices = data_dict['header']['method']['PVM_SPackArrNSlices'][0]
    
    norm_data = numpy.empty([x_size,y_size,num_slices,reps])

    for slice in range(num_slices):
        baseline = numpy.mean(data_dict['data'][:,:,slice,0:inj_point],axis=2)
        norm_data[:,:,slice,:] = (data_dict['data'][:,:,slice,:] / numpy.tile(baseline.reshape(x_size,y_size,1),reps))-1    

    return norm_data

def enhancement_curve(data_dict, mask=False):

    if mask:
        print "I don't know how to deal with masks quite yet"
        
        #TODO: Complete this to read in an ROI somehow
        
    else:
        print "Work in progress"

                
        return

def inj_point(data_dict):
    
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

    return injection_point[0][0]+1








