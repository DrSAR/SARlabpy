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
import sarpy
import pdb
import scipy.integrate
import sarpy.fmoosvi.getters as getters


def describe_object(class_object):

    try: # if patient, print studies
        for study in class_object.studies:
            print('-'*40+'\n'+study.subject.SUBJECT_id)    
    except: 
        try: # if multiple studies, print study names and scans         
            for study in class_object:
                print('-'*40+'\n'+study.subject.SUBJECT_id)    
                for scan in study.scans:
                    print("  "+scan.acqp.ACQ_protocol_name)
        except: 
            try: #if single study, print scan names             
                for scan in class_object.scans:
                    print("  "+scan.acqp.ACQ_protocol_name)      
            except: #otherwise print useless error message
                print("I don't know what's going on")

def calculate_AUC(scan_object, time = 60, pdata_num = 0):
    
    """
    Returns an area under the curve data for the scan object

    :param object scan_object: scan object from a study
    :param integer pdata_num: reconstruction number, according to python numbering.
            default reconstruction is pdata_num = 0.
    :return: array with auc data
    """ 
    
    ########### Getting and defining parameters
    
    # Visu_pars params
    x_size = scan_object.pdata[pdata_num].visu_pars.VisuCoreSize[0]
    y_size = scan_object.pdata[pdata_num].visu_pars.VisuCoreSize[1]
    num_slices = getters.get_num_slices(scan_object,pdata_num)
    phase_encodes = scan_object.pdata[pdata_num].visu_pars.VisuAcqPhaseEncSteps
  
    # Method params
    repetition_time = scan_object.method.PVM_RepetitionTime*1E-3

    # Calculated parms
    time_per_rep = repetition_time * phase_encodes
    auc_reps = int(numpy.round(time / time_per_rep))
    time_points = numpy.arange(time_per_rep,time_per_rep*auc_reps + time_per_rep,time_per_rep)
    
    ########### Start AUC code
      
    # Determine point of injection by averaging one slice in the entire image
    inj_point = sarpy.analysis.inj_point(scan_object)
    
    # Now calculate the Normalized Intesity voxel by voxel
    norm_data = normalize_dce(scan_object)

    # Now calculate the actual AUC
    auc_data = numpy.empty([x_size,y_size,num_slices])
    
    for slice in range(num_slices):
        auc_data[:,:,slice] = scipy.integrate.simps(norm_data[:,:,slice,inj_point:inj_point+auc_reps],x=time_points)
    return auc_data


def normalize_dce(scan_object, pdata_num = 0):

    ########### Getting and defining parameters
    
    # Visu_pars params
    x_size = scan_object.pdata[pdata_num].visu_pars.VisuCoreSize[0]
    y_size = scan_object.pdata[pdata_num].visu_pars.VisuCoreSize[1]
    num_slices = getters.get_num_slices(scan_object,pdata_num)

    # Method params
    reps =  scan_object.method.PVM_NRepetitions
    
    # Data
    data = scan_object.pdata[pdata_num].data
    
    # Calculated params      
    inj_point = sarpy.analysis.inj_point(scan_object)
    
    norm_data = numpy.empty([x_size,y_size,num_slices,reps])

    for slice in range(num_slices):
        baseline = numpy.mean(data[:,:,slice,0:inj_point],axis=2)
        norm_data[:,:,slice,:] = (data[:,:,slice,:] / numpy.tile(baseline.reshape(x_size,y_size,1),reps))-1    

    return norm_data

def enhancement_curve(data_dict, mask=False):

    if mask:
        print "I don't know how to deal with masks quite yet"
        
        #TODO: Complete this to read in an ROI somehow
        
    else:
        print "Work in progress"

                
        return

def inj_point(scan_object, pdata_num = 0):

    from collections import Counter   


    # Method params    
    num_slices = getters.get_num_slices(scan_object,pdata_num)
          
     # Data
    data = scan_object.pdata[pdata_num].data  

    try:      
        # Pool all the data together
        img_mean = data[:,:,:,:].sum(0).sum(0)

    except IndexError:
        print "You might only have 2D or 3D data, check data source! "
        raise IndexError
        

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








