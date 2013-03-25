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
import scipy.optimize
import scipy.fftpack
import sarpy.fmoosvi.getters as getters
import math


def h_describe_object(class_object):

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

def h_calculate_AUC(scan_object, time = 60, pdata_num = 0):
    
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
    inj_point = sarpy.fmoosvi.analysis.h_inj_point(scan_object)
    
    # Now calculate the Normalized Intesity voxel by voxel
    norm_data = h_normalize_dce(scan_object)

    # Now calculate the actual AUC
    auc_data = numpy.empty([x_size,y_size,num_slices])
    
    for slice in range(num_slices):
        auc_data[:,:,slice] = scipy.integrate.simps(norm_data[:,:,slice,inj_point:inj_point+auc_reps],x=time_points)
    return auc_data


def h_normalize_dce(scan_object, pdata_num = 0):

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
    inj_point = sarpy.analysis.h_inj_point(scan_object)
    
    norm_data = numpy.empty([x_size,y_size,num_slices,reps])

    for slice in range(num_slices):
        baseline = numpy.mean(data[:,:,slice,0:inj_point],axis=2)
        norm_data[:,:,slice,:] = (data[:,:,slice,:] / numpy.tile(baseline.reshape(x_size,y_size,1),reps))-1    

    return norm_data

def h_enhancement_curve(scan_object, pdata_num = 0, mask=False):

    if mask:
        print "I don't know how to deal with masks quite yet"
        
        #TODO: Complete this to read in an ROI somehow
        
    else:
        print "Work in progress"

                
        return

def h_inj_point(scan_object, pdata_num = 0):

    from collections import Counter   


    # Method params    
    num_slices = getters.get_num_slices(scan_object,pdata_num)
          
     # Data
    data = scan_object.pdata[pdata_num].data  

    try:      
        # Pool all the data together
        img_mean = data[:,:,:,:].sum(0).sum(0)

    except IndexError:
        print "You might only have 2D or 3D data, need 4D data check data source! "
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
    
def h_calculate_KBS(scan_object):
    
    KBS = 71.16*2
    # print('Still working on it')
    
    return KBS


    
def h_BS_B1map(zero_BSminus, zero_BSplus, high_BSminus, high_BSplus, scan_with_POI):
    
    try:
        TPQQ_POI = scan_with_POI.method.ExcPulse[3] #11.3493504066491 # 5.00591 #11.3493504066491
        pulse_width_POI = scan_with_POI.method.ExcPulse[0]*1E-3
    except:
        print('Please use a scan that has a valid power level for the pulse \
                of interest. scan_with_POI.method.ExcPulse[3]')
    
    TPQQ_BS = high_BSminus.method.BSPulse[3]
    
    integral_ratio = high_BSminus.method.ExcPulse[10] #0.071941 default from AY

    #TODO: Write function to calculateKBS
    KBS = h_calculate_KBS(high_BSminus)
    gamma = 267.513e6
    
    # Get phase data from fid
    offset = h_phase_from_fid(zero_BSplus) - h_phase_from_fid(zero_BSminus)
    phase_diff = h_phase_from_fid(high_BSplus) - h_phase_from_fid(high_BSminus) + offset
    
    # Calculate B1 peak
    B1peak = numpy.sqrt(numpy.absolute(phase_diff)/(2*KBS))
    
    # Calculate Flip Angle for the pulse of interest
    alpha_BS = (gamma*B1peak/10000) * (math.pow(10,(TPQQ_BS-TPQQ_POI)/20)) *\
                integral_ratio*pulse_width_POI


    return alpha_BS

def h_fit(x_data, y_data, fit_function, initial_params):    

    # Distance to the target function    
    errfunc = lambda param, x_data, y_data: fit_function(param, x_data) - y_data

    # Fit function to data    
    final_params, success = scipy.optimize.leastsq(errfunc, \
                            initial_params[:], args = (x_data, y_data))
    
    return final_params
        
def h_fit_T1_LL(scan_object, flip_angle_map = 0, pdata_num = 0):
       
    num_slices = getters.get_num_slices(scan_object, pdata_num)
    repetition_time = scan_object.method.PVM_RepetitionTime
    inversion_time = scan_object.method.PVM_InversionTime
    
    if type(flip_angle_map) != numpy.ndarray:
        flip_angle_map = math.radians(scan_object.acqp.ACQ_flip_angle)
   
    # Visu_pars params 
        
    data = scan_object.pdata[pdata_num].data 
    # data_for_fitting = numpy.zeros( [data.shape[0],data.shape[1],data.shape[2],data.shape[3]] )
    #for slice in range(num_slices):
        #data_for_fitting[:,:,slice,:] = numpy.mean(numpy.mean(data[:,:,slice,:],0),0)
    
    data_after_fitting = numpy.zeros( [data.shape[0],\
                                       data.shape[1],\
                                       data.shape[2]] )
            
    x_data = numpy.linspace(inversion_time,\
        scan_object.pdata[pdata_num].data.shape[3]*repetition_time,\
        scan_object.pdata[pdata_num].data.shape[3])
                                                       
                                                       
    fitfunc = lambda param, x_data: numpy.abs(param[0]*\
                                (1 - param[1]*numpy.exp(-x_data/param[2])))
    initial_params = [3E5, 2, 350]
    
    for x in xrange(data.shape[0]):
        for y in range(data.shape[1]):
            for slice in range(num_slices):
                
                y_data = data[x,y,slice,:]
                
                data_after_fitting[x,y,slice] = h_fit(x_data, y_data, \
                                                fitfunc, initial_params)[2]
    
    # Need to conert T1_eff to T1
    T1 = 1 / (( (1 / data_after_fitting) + numpy.log(numpy.cos(flip_angle_map))/repetition_time))
                                        
    return T1
                    
        #TODO: Implement code to deal with other methos of calculating T1
        # e.g., IR, VFA
        # NecS3Exp= sarpy.Experiment('NecS3')
        # scan_object = NecS3Exp.studies[0].find_scan_by_protocol('04_ubcLL2')
        
def h_phase_from_fid(scan_object):
    
    phase_data = numpy.angle(scipy.fftpack.fftshift(scipy.fftpack.fftn(scipy.fftpack.fftshift(scan_object.fid))))
    
    return phase_data
    
def h_mag_from_fid(scan_object):

    mag_data = numpy.abs(scipy.fftpack.fftshift(scipy.fftpack.fftn(scipy.fftpack.fftshift(scan_object.fid))))
    
    return mag_data
    
    
    







