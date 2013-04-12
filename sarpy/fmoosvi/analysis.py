# -*- coding: utf-8 -*-
# Copyright (C) 2012-2013 Stefan A Reinsberg and SARlab members
# full license details see LICENSE.txt
"""Collection of analysis routines for Bruker data
test

"""

from __future__ import division

import numpy
import sarpy
import scipy.integrate
import scipy.optimize
import scipy.fftpack
import sarpy.fmoosvi.getters as getters
import math



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
    auc_data = numpy.empty([norm_data.shape[0],norm_data.shape[1],num_slices])
    
    for slice in range(num_slices):
        auc_data[:,:,slice] = scipy.integrate.simps(norm_data[:,:,slice,inj_point:inj_point+auc_reps],x=time_points)
    return auc_data


def h_normalize_dce(scan_object, pdata_num = 0):

    ########### Getting and defining parameters
    
    # Data
    data = scan_object.pdata[pdata_num].data
        
    # Visu_pars params
    x_size = data.shape[0]
    y_size = data.shape[1]
    num_slices = getters.get_num_slices(scan_object,pdata_num)

    # Method params
    reps =  scan_object.method.PVM_NRepetitions

    # Calculated params      
    inj_point = sarpy.fmoosvi.analysis.h_inj_point(scan_object)
    
    norm_data = numpy.empty([x_size,y_size,num_slices,reps])

    for slice in range(num_slices):
        baseline = numpy.mean(data[:,:,slice,0:inj_point],axis=2)
        norm_data[:,:,slice,:] = (data[:,:,slice,:] / numpy.tile(baseline.reshape(x_size,y_size,1),reps))-1    

    return norm_data
 
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
            injection_point.append(next(i for i,v in enumerate(diff_slice) if v > 2*std_slice))
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
    
    
    
##############
    #
    # Fitting section
    #
    #
    #
    #
    #
    #
    #
##############    

def h_func_T1(params,t):
    M,B,T1_eff = params
    return numpy.abs(M*(1-B*numpy.exp(-t/T1_eff))) #edit 2: funcion didn't have an abs!! FAIL!
    
def h_within_bounds(params,bounds):
    try:        
        for p in xrange(len(params)):
            if params[p] >= bounds[p,0] or params[p] <= bounds[p,1] :
                return True
        return False
    except:
        print('You have some funky inputs for the bounds, you fail.')
        return False # edit 1 to fix fitting
            

def h_residual_T1(params, y_data, t):
    
    bounds = numpy.zeros(shape=[3,2])
    bounds[0,0] = 1e2
    bounds[0,1] = 1e8
    bounds[1,0] = 0
    bounds[1,1] = 5
    bounds[2,0] = 50
    bounds[2,1] = 10000

    if h_within_bounds(params,bounds):
        return y_data - h_func_T1(params, t)
    else:
        return 1e9

def h_fit_T1_LL(scan_object, flip_angle_map = 0, pdata_num = 0, 
                params = []):
    
    if len(params) == 0:      
        params = [3E5, 2, 350]

    num_slices = getters.get_num_slices(scan_object, pdata_num)
    repetition_time = scan_object.method.PVM_RepetitionTime
    inversion_time = scan_object.method.PVM_InversionTime
    
    if type(flip_angle_map) != numpy.ndarray:
        flip_angle_map = math.radians(scan_object.acqp.ACQ_flip_angle)
   
    # Visu_pars params 
        
    data = scan_object.pdata[pdata_num].data[:]

    data_after_fitting = numpy.zeros( [data.shape[0],\
                                       data.shape[1],\
                                       data.shape[2]] )
                                       
    fit_results = numpy.array(data_after_fitting[:], dtype=dict)                       
            
    t_data = numpy.linspace(inversion_time,\
        scan_object.pdata[pdata_num].data.shape[3]*repetition_time,\
        scan_object.pdata[pdata_num].data.shape[3])
  
    for x in xrange(data.shape[0]):
        for y in range(data.shape[1]):
            for slice in range(num_slices):
                
                y_data = data[x,y,slice,:]
                fit_dict = {}
                
                fit_params,cov,infodict,mesg,ier = scipy.optimize.leastsq(h_residual_T1,params,args=(y_data,t_data), full_output = True,maxfev = 200)

                goodness_of_fit = h_goodness_of_fit(y_data,infodict)
                
                [M,B,T1_eff] = fit_params
                fit_dict = {
                            'fit_params': fit_params,
                            'cov' : cov,
                            'infodict' : infodict,
                            'mesg' : mesg,
                            'ier' : ier,
                            'goodness': goodness_of_fit
                            }
                            
                data_after_fitting[x,y,slice] = T1_eff
                fit_results[x,y,slice] = fit_dict
    
    # Need to convert T1_eff to T1
    T1 = 1 / (( (1 / data_after_fitting) + numpy.log(numpy.cos(flip_angle_map)) / repetition_time))
    
    # Make absurd values nans to make my life easer:
    T1[T1<0] = numpy.nan
    T1[T1>1e4] = numpy.nan

    return T1, fit_results
                    
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
    
def h_image_to_mask(scan_object, adata_key):
    
    img_data = scan_object.adata[adata_key].data.get_data()

    masked_data = img_data[:] 
    
    mask_val = scipy.percentile(masked_data.flatten(),95)
    masked_data[masked_data == mask_val] = numpy.nan
    masked_data[numpy.isfinite(masked_data)] = 1
    
    return masked_data    

def h_goodness_of_fit(data,infodict, indicator = 'rsquared'):
    
    if indicator == 'rsquared':
        ss_err=(infodict['fvec']**2).sum()
        ss_tot=((data-data.mean())**2).sum()
        rsquared=1-(ss_err/ss_tot)
        
        print type(rsquared)
        
        return rsquared
        
    else:
        print ('There is no code to produce that indicator. Do it first.')
        raise Exception
        
    
