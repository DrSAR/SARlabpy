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
import scipy.stats
import sarpy.fmoosvi.getters as getters
import sarpy.ImageProcessing.resample_onto
import sarpy.io
import math
import copy
import os
import json
import nibabel
import datetime
import collections
import random
import copy

def h_calculate_AUC(scn_to_analyse=None, 
                    bbox = None, 
                    time = 60, 
                    pdata_num = 0,
                    **kwargs):
    
    """
    Returns an area under the curve data for the scan object

    :param object scan_object: scan object from a study
    :param integer pdata_num: reconstruction number, according to python numbering.
            default reconstruction is pdata_num = 0.
    :return: array with auc data
    """ 
    import sarpy.fmoosvi.analysis            

    
    scan_object = sarpy.Scan(scn_to_analyse)

    ########### Getting and defining parameters
    
    # Visu_pars params
    num_slices = getters.get_num_slices(scn_to_analyse,pdata_num)
    
    reps =  scan_object.method.PVM_NRepetitions
    
    if reps != scan_object.pdata[pdata_num].data.shape[-1]:
        reps = scan_object.pdata[pdata_num].data.shape[-1]
        print('\n \n ***** Warning **** \n \n !!! Incomplete dce data for {0}'.format(scan_object.shortdirname) )
    
    # there are problems with using phase encodes for certain cases (maybe 3D)
    # so now I have to use the tuid time
    total_time = scan_object.method.PVM_ScanTimeStr
    format = "%Hh%Mm%Ss%fms"
    t=datetime.datetime.strptime(total_time,format)
    total_time = (3600*t.hour) + (60*t.minute) + (t.second) + t.microsecond*1E-6

    time_per_rep = numpy.round(numpy.divide(total_time,reps))
    
    try:  
        auc_reps = int(numpy.round(time / time_per_rep))
        
        if auc_reps == 0:
            raise ZeroDivisionError      
    except ZeroDivisionError:
        print('h_calculate_auc: Insufficient data for AUC (0 reps) in scan {0}'.format(scan_object.shortdirname))
        raise ZeroDivisionError
        
    time_points = numpy.arange(time_per_rep,time_per_rep*auc_reps + time_per_rep,time_per_rep)

    ########### Start AUC code
      
    # Determine point of injection by averaging one slice in the entire image
    inj_point = sarpy.fmoosvi.analysis.h_inj_point(scn_to_analyse)
    
    # Now calculate the Normalized Intesity voxel by voxel
    norm_data = h_normalize_dce(scn_to_analyse)

    # Size info
    x_size = norm_data.shape[0]
    y_size = norm_data.shape[1]
    num_slices = norm_data.shape[2]    

    # Now calculate the actual AUC
    auc_data = numpy.empty([x_size,y_size,num_slices])
    
    for slice in range(num_slices):
        if inj_point+auc_reps <= reps:
            auc_data[:,:,slice] = scipy.integrate.simps(norm_data[:,:,slice,inj_point:inj_point+auc_reps],x=time_points)
        else:
            auc_data[:,:,slice] = numpy.nan

       
    # Deal with bounding boxes

    if bbox is None:        
        bbox = numpy.array([0,x_size-1,0,y_size-1])    

    else:
        bbox = sarpy.fmoosvi.getters.convert_bbox(scn_to_analyse,bbox) 

    if bbox.shape == (4,):            
    
        bbox_mask = numpy.empty([x_size,y_size])
        bbox_mask[:] = numpy.nan        
        bbox_mask[bbox[0]:bbox[1],bbox[2]:bbox[3]] = 1
    
        # First tile for slice
        bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size,1),num_slices) 
     
    return {'':auc_data*bbox_mask}


def h_normalize_dce(scn_to_analyse=None, bbox = None, pdata_num = 0):

    scan_object = sarpy.Scan(scn_to_analyse)

    ########### Getting and defining parameters
    
    # Data
    data = scan_object.pdata[pdata_num].data

    x_size = data.shape[0]
    y_size = data.shape[1]
    num_slices = getters.get_num_slices(scn_to_analyse,pdata_num)
    
    # Method params
    #TODO: change this so it doesn't require method file WIHOUT BREAKING IT!
    reps =  scan_object.method.PVM_NRepetitions
    
    if reps != scan_object.pdata[pdata_num].data.shape[-1]:
        reps = scan_object.pdata[pdata_num].data.shape[-1]
        print('\n \n ***** Warning **** \n \n !!! Incomplete dce data for {0} \n \n'.format(scan_object.shortdirname) )

    ## Check for bbox traits and create bbox_mask to output only partial data

    if bbox is None:        
        bbox = numpy.array([0,x_size-1,0,y_size-1])
    else:      
        bbox = sarpy.fmoosvi.getters.convert_bbox(scn_to_analyse,bbox) 

    if bbox.shape == (4,):            
    
        bbox_mask = numpy.empty([x_size,y_size])

        bbox_mask[:] = numpy.nan        
        bbox_mask[bbox[0]:bbox[1],bbox[2]:bbox[3]] = 1
    
        # First tile for slice
        bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size,1),num_slices)
        # Next tile for reps
        bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size,num_slices,1),reps)


    # Calculated params      
    inj_point = sarpy.fmoosvi.analysis.h_inj_point(scn_to_analyse)
    
    norm_data = numpy.empty([x_size,y_size,num_slices,reps])


    for slice in range(num_slices):
        baseline = numpy.mean(data[:,:,slice,0:inj_point],axis=2)
        norm_data[:,:,slice,:] = (data[:,:,slice,:] / numpy.tile(baseline.reshape(x_size,y_size,1),reps))-1

    return norm_data*bbox_mask
 
def h_enhancement_curve(scn_to_analyse=None, 
                        adata_roi_label=None, 
                        pdata_num = 0,
                        **kwargs):

    scan_object = sarpy.Scan(scn_to_analyse)

    norm_data = sarpy.fmoosvi.analysis.h_normalize_dce(scn_to_analyse)
    num_slices = norm_data.shape[-2]
    reps = norm_data.shape[-1]
    
    if reps != scan_object.pdata[pdata_num].data.shape[-1]:
        reps = scan_object.pdata[pdata_num].data.shape[-1]
        print('\n \n ***** Warning **** \n \n !!! Incomplete dce data for {0}'.format(scan_object.shortdirname) )

    # there are problems with using phase encodes for certain cases (maybe 3D)
    # so now I have to use the tuid time
    total_time = scan_object.method.PVM_ScanTimeStr
    format = "%Hh%Mm%Ss%fms"
    t=datetime.datetime.strptime(total_time,format)
    total_time = (3600*t.hour) + (60*t.minute) + (t.second) + t.microsecond*1E-6

    time_per_rep = numpy.divide(total_time,reps)
        
    #Calculating the time per rep.
    time = numpy.linspace(0,reps-1,num=reps)*time_per_rep
                   
    ## THIS IS INCREDIBLY SKETCHY, AND I'M NOT SURE WHAT THE RAMIFICATIONS ARE        
    new_scan_object = copy.deepcopy(scan_object)
    new_scan_object.pdata[pdata_num].data = norm_data
    data_scan = new_scan_object.pdata[pdata_num]
    ## END SKETCHY BIT
    
    roi = scan_object.adata[adata_roi_label].data
     
    masked_data = data_scan.data * numpy.tile(numpy.reshape(roi,[roi.shape[0], roi.shape[1], roi.shape[2],1]),reps)

    enhancement_curve = numpy.empty(shape = [num_slices, 2, reps])
    
    for slice in range(num_slices):
        
        enhancement_curve[slice,0,:] = time
        enhancement_curve[slice,1,:] = scipy.stats.nanmean(scipy.stats.nanmean(masked_data[:,:,slice,:], axis=0), axis =0)   

    return {'':enhancement_curve}


def h_inj_point(scn_to_analyse=None, pdata_num = 0):

    scan_object = sarpy.Scan(scn_to_analyse)

    from collections import Counter   

    # Method params    
    num_slices = getters.get_num_slices(scn_to_analyse,pdata_num)
          
     # Data
    data = scan_object.pdata[pdata_num].data  

    try:      
        # Pool all the data together
        img_mean = data[:,:,:,:].sum(0).sum(0)

    except IndexError:
        print('h_inj_point: Scan {0}: You might only have 2D or 3D data, need 4D data check data source!'.format(scan_object.shortdirname))
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
            
    # look through the list of elements in injection point and report the most common (mode)
    injection_point_counter = Counter(injection_point)
    injection_point = injection_point_counter.most_common(1)

    return injection_point[0][0]+1

def h_calculate_AUGC(scn_to_analyse=None, 
                     adata_label=None, 
                     bbox = None, 
                     time = 60, 
                     pdata_num = 0,
                     **kwargs):
    
    """
    Returns an area under the gadolinium concentration curve adata for the scan object

    :param object scan_object: scan object from a study
    :param str adata_label: indicates the label with which gd_conc is stored
    :param integer pdata_num: reconstruction number, according to python numbering.
            default reconstruction is pdata_num = 0.
    :return: array with augc data
    """
    scan_object = sarpy.Scan(scn_to_analyse)   

    # Get the concentration data stored as an adata

    try:
        data = scan_object.adata[adata_label].data
    except KeyError:
        print('h_caculate_AUGC: Source data {0} does not exist yet.'.format(adata_label))
        raise KeyError
    
    ########### Getting and defining parameters
    
    # Visu_pars params
    num_slices = getters.get_num_slices(scn_to_analyse,pdata_num)
    phase_encodes = scan_object.pdata[pdata_num].visu_pars.VisuAcqPhaseEncSteps
  
    # Method params
    repetition_time = scan_object.method.PVM_RepetitionTime*1E-3
    
    reps =  scan_object.method.PVM_NRepetitions
    
    if reps != scan_object.pdata[pdata_num].data.shape[-1]:
        reps = scan_object.pdata[pdata_num].data.shape[-1]
        print('\n \n ***** Warning **** \n \n !!! Incomplete dce data for {0} \n \n '.format(scan_object.shortdirname) )
    
    # there are problems with using phase encodes for certain cases (maybe 3D)
    # so now I have to use the tuid time
    total_time = scan_object.method.PVM_ScanTimeStr
    format = "%Hh%Mm%Ss%fms"
    t=datetime.datetime.strptime(total_time,format)
    total_time = (3600*t.hour) + (60*t.minute) + (t.second) + t.microsecond*1E-6

    time_per_rep = numpy.divide(total_time,reps)

    # Calculated parms
    augc_reps = int(numpy.round(time / time_per_rep))
    time_points = numpy.arange(time_per_rep,time_per_rep*augc_reps + time_per_rep,time_per_rep)

    ########### Start AUGC code
      
    # Determine point of injection by averaging one slice in the entire image
    inj_point = sarpy.fmoosvi.analysis.h_inj_point(scn_to_analyse)
    
    # Size info
    x_size = data.shape[0]
    y_size = data.shape[1]
    num_slices = data.shape[2]
    
    # Now calculate the actual AUGC
    augc_data = numpy.empty([x_size,y_size,num_slices])
    
    for slice in range(num_slices):
        augc_data[:,:,slice] = scipy.integrate.simps(data[:,:,slice,inj_point:inj_point+augc_reps],x=time_points)
    
    # Deal with bounding boxes

    if bbox is None:        
        bbox = numpy.array([0,x_size-1,0,y_size-1])

    else:      
        bbox = sarpy.fmoosvi.getters.convert_bbox(scn_to_analyse,bbox)           
       
    if bbox.shape == (4,):            
    
        bbox_mask = numpy.empty([x_size,y_size])
        bbox_mask[:] = numpy.nan        
        bbox_mask[bbox[0]:bbox[1],bbox[2]:bbox[3]] = 1
    
        # First tile for slice
        bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size,1),num_slices)

        # If this gives a value error about operands not being broadcast together, go backand change your adata to make sure it is squeezed
    return {'':augc_data*bbox_mask}
    
def h_conc_from_signal(scn_to_analyse=None, 
                       other_scan_name=None, 
                       adata_label = None,
                       bbox = None,
                       relaxivity=4.3e-3, 
                       pdata_num = 0,
                       masterlist_name=None,
                       other_dce_label = None,
                       **kwargs):

    scan_object = sarpy.Scan(scn_to_analyse)
    scan_object_T1map = sarpy.Scan(other_scan_name)
    ########### Getting and defining parameters
    
    # Data
    data = scan_object.pdata[pdata_num].data
    
    # resample the t1map onto the dce
   
    data_t1map = sarpy.ImageProcessing.resample_onto.resample_onto_pdata(scan_object_T1map.adata[adata_label],
                                                                         scan_object.pdata[pdata_num],
                                                                         use_source_dims=True,
                                                                         replace_nan=0)
    x_size = data.shape[0]
    y_size = data.shape[1]
    num_slices = getters.get_num_slices(scn_to_analyse,pdata_num)
    
    # Method params
    #TODO: change this so it doesn't require method file WIHOUT BREAKING IT!
    reps =  scan_object.method.PVM_NRepetitions

    if reps != scan_object.pdata[pdata_num].data.shape[-1]:
        reps = scan_object.pdata[pdata_num].data.shape[-1]
        print('\n \n ***** Warning **** \n \n !!! Incomplete dce data for {0} \n \n'.format(scan_object.shortdirname) )
    
    total_time = scan_object.pdata[pdata_num].visu_pars.VisuAcqScanTime / 1000.
    time_per_rep = numpy.round(numpy.divide(total_time,reps))

    TR = scan_object.method.PVM_RepetitionTime
    FA = scan_object.acqp.ACQ_flip_angle
    
    inj_point = sarpy.fmoosvi.analysis.h_inj_point(scn_to_analyse, pdata_num = 0)    
    
    # Deal with bounding boxes
    if bbox is None:        
        bbox = numpy.array([0,x_size-1,0,y_size-1])    
    else:      
        bbox = sarpy.fmoosvi.getters.convert_bbox(scn_to_analyse,bbox) 

    if bbox.shape == (4,):            
    
        bbox_mask = numpy.empty([x_size,y_size])
        bbox_mask[:] = numpy.nan        
        bbox_mask[bbox[0]:bbox[1],bbox[2]:bbox[3]] = 1
    
        # First tile for slice
        bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size,1),num_slices)

    if other_dce_label:

        dic = h_stitch_dce_scans(masterlist_name=masterlist_name,
                               bbox=bbox,
                               scn_to_analyse=scn_to_analyse, 
                               other_dce_label=other_dce_label,
                               pdata_num = 0)

        source_data = dic['']
        time = dic['_time']
        fulltime = dic['_fulltime']
        idx = dic['_idx']

    else:
        source_data = data
        time = numpy.linspace(0,reps-1,num=reps)*time_per_rep
        fulltime = time
        idx = numpy.arange(0,time.shape[-1])

    T1 = numpy.empty_like(source_data) + numpy.nan
    
    for x in xrange(bbox[0],bbox[1]):
        for y in xrange(bbox[2],bbox[3]):
            for slice in xrange(num_slices):

                baseline_s = numpy.mean(source_data[x,y,slice,0:inj_point])
                E1 = numpy.exp(-TR/data_t1map[x,y,slice])
                c = numpy.cos(numpy.radians(FA))
                T1[x,y,slice,0:inj_point] = data_t1map[x,y,slice]

                # Use the SPGR equation twice, once before agent admin. and once after. If you
                # divide the two, the M0s cancel out, and so do the sin thetas. what remains is:
                # (1-E1)*(1-E0*cos) / ((1-E0) * (1-E1*cos)). use wolfram alpha to solve this:
                # http://www.wolframalpha.com/input/?i=solve+%28%281-x%29*%281-b*c%29%29+%2F+%28%281-b%29*%281-x*c%29%29+%3D+r+for+x

                for rep in xrange(inj_point,source_data.shape[-1]):                    
                    s = source_data[x,y,slice,rep] / baseline_s
                    E2 = (-E1*c + E1*s - s + 1) / (E1*s*c - E1*c - s*c +1)
                    T1[x,y,slice,rep] = -TR / numpy.log(E2)


    # If this gives a value error about operands not being broadcast together, go backand change your adata to make sure it is squeezed

    T1baseline = numpy.squeeze(data_t1map)*bbox_mask
    T1baseline = numpy.tile(T1baseline.reshape(x_size,y_size,num_slices,1),source_data.shape[-1])
    conc = (1/relaxivity) * ( (1/T1) - (1/T1baseline) )
   
    #conc[conc<0] = 0

    times = {'time':time,
            'fulltime':fulltime,
            'idx': idx.astype(int)}
       
    return {'':conc,
            '_times':times}
    
def h_stitch_dce_scans(masterlist_name,
                        bbox,
                       scn_to_analyse=None, 
                       other_dce_label=None,
                       pdata_num = 0,
                       **kwargs):

    # Sigh, here we have to break a cardinal rule not to have anything to do with masterlists in my 
    # analysis functions, but alas.

    # Masterlist reading

    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)

    fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str)

    scanA = sarpy.Scan(scn_to_analyse)
    timeA = scanA.acqp.ACQ_time # time dcea finished
    date_objectA = datetime.datetime.strptime(timeA, '%H:%M:%S %d %b %Y')
    patient = scanA.pdata[pdata_num].visu_pars.VisuSubjectId
    total_timeA = scanA.pdata[pdata_num].visu_pars.VisuAcqScanTime / 1000. # total scan time in seconds
    repsA = scanA.pdata[pdata_num].data.shape[-1]
    time_per_repA = numpy.divide(total_timeA,repsA)

    scanB = sarpy.Scan(master_list[patient][other_dce_label][0])
    timeB = scanB.acqp.ACQ_time # time dceb finished
    date_objectB = datetime.datetime.strptime(timeB, '%H:%M:%S %d %b %Y')
    total_timeB = scanB.pdata[pdata_num].visu_pars.VisuAcqScanTime / 1000. # total scan time in seconds
    repsB = scanB.pdata[pdata_num].data.shape[-1]
    duration = scanB.pdata[pdata_num].visu_pars.VisuAcqScanTime / 1000. #duration of dceb scan
    time_per_repB = numpy.divide(total_timeB,repsB)    
    deltaT = date_objectB-date_objectA
    added_reps = numpy.round(numpy.divide(deltaT.seconds - duration,time_per_repA))

    dataShape = scanA.pdata[pdata_num].data.shape

    timeA = numpy.linspace(0,repsA-1,num=repsA)*time_per_repA
    timeB = numpy.linspace(0,repsB-1,num=repsB)*time_per_repB + (repsA-1+added_reps)*time_per_repA
    timeAB = numpy.hstack((timeA,timeB))

    timeAtotB = numpy.linspace(0,(repsA+repsB+added_reps)-1,num=repsA+repsB+added_reps)*time_per_repA

    if time_per_repA != time_per_repB:
        raise NotImplementedError('The time per reps of both scans need to be the same')

    gd_conc_stitched = numpy.empty(shape=[dataShape[0],dataShape[1],dataShape[2],repsA+repsB]) + numpy.nan

    idx = numpy.hstack((numpy.arange(0,repsA),numpy.arange(repsA+added_reps,repsA+added_reps+repsB)))

    for x in xrange(bbox[0],bbox[1]):
        for y in xrange(bbox[2],bbox[3]):
            for slice in xrange(dataShape[2]):

                dA = scanA.pdata[pdata_num].data[x,y,slice,:]
                dB = scanB.pdata[pdata_num].data[x,y,slice,:]

                gd_conc_stitched[x,y,slice,:] = numpy.hstack((dA,dB))

    return {'':gd_conc_stitched,
            '_time': timeAB,
            '_fulltime': timeAtotB,
            '_idx':idx}
    
### T1 fitting section
    
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

def h_within_bounds(params,bounds):

    for p in xrange(len(params)):
        if params[p] >= bounds[p,0] or params[p] <= bounds[p,1] :
            return True
    return False


def h_func_T1_FAind(params,t_data):
    a,b,T1_eff,phi = params
    #print type(a), type(b), type(n), type(phi), type(T1_eff)

    res = a*(1-(1-b)*numpy.exp(-t_data/T1_eff))*numpy.exp(1j*phi)

    if numpy.isfinite(res).all():
        return res
    else:
        return [1e30]*t_data.shape[-1]

def T1eff_to_T1(T1,Td, Tp, tau, T1_eff, b, c):    
    return (1-2*numpy.exp(-Td/T1) + numpy.exp(-(Tp+Td)/T1))*((1-numpy.exp(-tau/T1_eff))/(1-numpy.exp(-tau/T1))) -c*(numpy.exp(-tau/T1_eff)/(numpy.exp(-tau/T1)))*numpy.exp(-(Tp+Td)/T1) - b 
             
def h_residual_T1_FAind(params, y_data, t_data):
    
    bounds = numpy.zeros(shape=[4,2])
    bounds[0,0] = 1e2
    bounds[0,1] = 1e8
    bounds[1,0] = -1e3
    bounds[1,1] = 1e3
    bounds[2,0] = 20
    bounds[2,1] = 3000
    bounds[3,0] = -5
    bounds[3,1] = +5

    if h_within_bounds(params,bounds):
        return numpy.abs(y_data - h_func_T1_FAind(params,t_data))
    else:
        return 1e30


def h_fitpx_T1_LL_FAind(scn_to_analyse=None,
                        y_data=None,
                        slice=None,
                        initial_params = None,
                        fit_algorithm=None,
                        **kwargs):

    scan_object = sarpy.Scan(scn_to_analyse)

    if fit_algorithm is None:
        fit_algorithm = 'leastsq'

    # Get the parameters necessary for fitting the individual pixel

    tau = float(scan_object.method.PVM_RepetitionTime)
    total_TR = float(scan_object.method.Inv_Rep_time)
    Nframes = float(scan_object.method.Nframes)
    delay = float(scan_object.method.InterSliceDelay)
    FA = float(numpy.radians(scan_object.acqp.ACQ_flip_angle))    
    inv_time = float(scan_object.method.PVM_InversionTime)
    # The following are slice dependent parameters

    Td = delay*(slice+1.0)
    Tp = total_TR - (Nframes -1.0)*tau - Td    
    #TODO: make this better; slice dependent
    t_data = numpy.arange(0,Nframes) * tau + inv_time     

    fit_dict = {}

    if initial_params is None:

        params = h_fit_T1_LL_FAind_initial_params_guess(y_data,t_data)

    else:
        params = numpy.ones(4)
        try:
            params = initial_params
        except ValueError:
            print('Please supply a list of 4 starting parameters')
            raise

    # Step 1: Fit Eq.1 from Koretsky paper for a,b,T1_eff, and 
    # phi (phase factor to fit real data)

    if fit_algorithm == 'leastsq':

        try:

            fit_params,cov,infodict,mesg,ier = scipy.optimize.leastsq(
                                                    h_residual_T1_FAind,
                                                    params,
                                                    args=(y_data, t_data), 
                                                    full_output = True,
                                                    maxfev = 200)
        except ValueError:
            # This is going to fail because when i moved the innards into a different function, x and y are no longer defined
            #TODO fix this because this might actually be useful.
            #print "Firas, look here fix this because it might actually be useful to know where I'm failing."
            #print x,y,slice
            raise
    else:
        print('fit type {0} not implemented yet'.format(fit_algorithm))

    [a,b,T1_eff,phi] = fit_params

    # Step 2: Calculate M(N-1)/M(inf) from Eq.1 from Koretsky paper
    # Divide both sides by M(inf) and then use M(0)/M(inf) -> step1

    c = 1- (1-b)*numpy.exp(-(Nframes-1)*tau/T1_eff)

    # Step 3: Solve Equation 6 to get T1 from it. Try using Newton-
    # Rhapsod method
   
    calc_params = (Td, Tp, tau, T1_eff,b,c)

    try:                    
        
        # Using Bisect method (other possibilities include newton, brentq:
        T1 = scipy.optimize.bisect(T1eff_to_T1,5,4000,
                                   args = (calc_params))

    except ValueError:
        T1=numpy.nan
        pass
    
    # Make absurd values nans to make my life easer:
    if (T1<0 or T1>=1e4):
        T1 = numpy.nan

    if fit_algorithm == 'leastsq':
        return infodict,mesg,ier,fit_params, T1,t_data

    else:
        return a,b,T1_eff,phi,T1, t_data

def h_fit_T1_LL_FAind_initial_params_guess(y_data,t_data):

    params = numpy.ones(4)  

    ## Guesses at parameters if none provided
    params[3] = numpy.angle(numpy.mean(y_data[-5:])) #phi
    params[0] = numpy.abs(numpy.mean(y_data[-5:])) #a
    params[1] = numpy.real(numpy.divide(y_data[0]*numpy.exp(-1j*params[3]),params[0])) #b

    # Getting a good guess for T1eff by finding the zero crossing
    yp = y_data[0:5]
    xp = t_data[0:5]

    A = numpy.array([xp, numpy.ones(xp.shape[0])])
    w=numpy.linalg.lstsq(A.T,yp)[0]
    params[2] = numpy.real(numpy.divide(-w[1],w[0]))

    # To prevent negative T1_eff
    if params[2] < 0:
        params[2] = 50

    return params


def h_fit_T1_LL_FAind(scn_to_analyse=None, 
                      bbox = None, 
                      pdata_num = 0, 
                      params = [],
                      fit_algorithm = None,
                      **kwargs):

    scan_object = sarpy.Scan(scn_to_analyse)
   
    ## Setting parameters
    data = scan_object.fftfid
    num_slices = getters.get_num_slices(scn_to_analyse,pdata_num)                                        

    x_size = data.shape[0]
    y_size = data.shape[1]
                                       
    # Initializations

    if fit_algorithm is None:
        fit_algorithm = 'leastsq'

    data_after_fitting = numpy.nan + numpy.empty( [x_size,y_size,num_slices] )

    fit_params1 = numpy.empty_like(data_after_fitting,dtype=dict)

    a_arr = numpy.empty_like(data_after_fitting)+numpy.nan
    b_arr = numpy.empty_like(data_after_fitting)+numpy.nan
    T1_eff_arr = numpy.empty_like(data_after_fitting)+numpy.nan
    phi_arr = numpy.empty_like(data_after_fitting)+numpy.nan
    T1_arr = numpy.empty_like(data_after_fitting)+numpy.nan

    infodict1 = numpy.empty_like(data_after_fitting,dtype=dict)
    mesg1 = numpy.empty_like(data_after_fitting,dtype=str)
    ier1 = numpy.empty_like(data_after_fitting)
    goodness_of_fit1 = numpy.empty_like(data_after_fitting)

    ## Check for bbox traits and create bbox_mask to output only partial data
    if bbox is None:        
        bbox = numpy.array([0,x_size-1,0,y_size-1])
    else:
        bbox = sarpy.fmoosvi.getters.convert_bbox(scn_to_analyse,bbox)  
        
    # Start the fitting process        

    for slice in xrange(num_slices):          
        for x in xrange(bbox[0],bbox[1]):
            for y in xrange(bbox[2],bbox[3]):


                # Deal with scans that have only one slice
                if num_slices > 1:
                    y_data = data[x,y,slice,:]
                else:
                    y_data = data[x,y,:]

                [infodict,mesg,ier,res_fit_params,T1,t_data] = h_fitpx_T1_LL_FAind(scn_to_analyse,
                                                                                   y_data=y_data,
                                                                                   slice=slice,
                                                                                   fit_algorithm=fit_algorithm)

                [a,b,T1_eff,phi] = res_fit_params

                a_arr[x,y,slice] = a
                b_arr[x,y,slice] = b
                T1_eff_arr[x,y,slice] = T1_eff
                phi_arr[x,y,slice] = phi
                T1_arr[x,y,slice] = T1

                data_after_fitting[x,y,slice] = T1

                if fit_algorithm == 'leastsq':   
                    infodict1[x,y,slice] = infodict
                    mesg1[x,y,slice] = mesg
                    goodness_of_fit1[x,y,slice] = h_goodness_of_fit(y_data,infodict)  

                ier1[x,y,slice] = ier


### this will need some serious debugging

## Moved the entire T1 fitting process into a separate function

    fit_params1 = {'a': a_arr,
                   'b': b_arr,
                   'T1_eff': T1_eff_arr,
                   'phi': phi_arr}

    if fit_algorithm == 'leastsq':

        result = {'rawdata':data,
                  't_data':t_data,
                  'params':fit_params1,
                  'infodict':infodict1,
                  'ier':ier1,
                  'goodness':goodness_of_fit1,
                  'mesg':mesg1}

        return {'':numpy.squeeze(data_after_fitting),
                '_fit':result}

    elif fit_algorithm == 'fmin':

        result = {'params': fit_params1,
                'iter':ier1}         

        return {'':numpy.squeeze(data_after_fitting),
                '_fit':result}    

### Other Helpers

def h_phase_from_fid(scn_to_analyse=None):

    scan_object = sarpy.Scan(scn_to_analyse)

    phase_data = numpy.angle(scipy.fftpack.fftshift(scipy.fftpack.fftn(scipy.fftpack.fftshift(scan_object.fid))))
    
    return phase_data
    
def h_mag_from_fid(scn_to_analyse=None):

    scan_object = sarpy.Scan(scn_to_analyse)

    mag_data = numpy.abs(scipy.fftpack.fftshift(scipy.fftpack.fftn(scipy.fftpack.fftshift(scan_object.fid))))
    
    return mag_data
    
def h_image_to_mask(roi_data, background=None, foreground=None,  peaks = None):
    
    if background is None:
        background = numpy.nan
    if foreground is None:
        foreground = 1        
    if peaks is None:
        peaks = 1
    else: 
        peaks =numpy.int(peaks)

    roi_mask = copy.deepcopy(roi_data)
    
    if peaks == 1:
        for slice in xrange(roi_mask.shape[2]):
        
            curr_slice = roi_mask[:,:,slice]
            
            # the most common value will be the background; We assume the ROI 
            # occupies only a small region in the image (less than 50%). 
            # By choosin the median we could have a few pixel values higher or 
            # lower than the very common background pixel intensity.
            mask_val = scipy.median(curr_slice.flatten())
            places = numpy.where(curr_slice == mask_val)
            notplaces = numpy.where(curr_slice != mask_val)
            curr_slice[places] = background
            curr_slice[notplaces] = foreground

        return roi_mask    
                
    else:
                          
        for slice in xrange(roi_mask.shape[2]):
            
            for x in xrange(peaks):
        
                curr_slice = roi_mask[:,:,slice]
              
                if x==0: # On the first slice, do the conventionl method
                
                    # the most common value will be the background; We assume the ROI 
                    # occupies only a small region in the image (less than 50%). 
                    # By choosin the median we could have a few pixel values higher or 
                    # lower than the very common background pixel intensity.
                    mask_val = scipy.median(curr_slice[numpy.isfinite(curr_slice)].flatten())
                    places = numpy.where(curr_slice == mask_val)
                    curr_slice[places] = numpy.nan
                else:
                    # the conventional method doesn't work for multiple peaks 
                    # because the next highest peak is to close to the other
                    # values. This is best, but is susceptible to throwing away
                    # data if th peak are set too high

                    # If/else block is needed to skip all the situations where
                    # the array is completely empty

                    if curr_slice[numpy.isfinite(curr_slice)].flatten() !=[]:
                        mask_val = scipy.stats.mode(curr_slice[numpy.isfinite(curr_slice)].flatten())[0]
                        places = numpy.where(curr_slice == mask_val)
                        curr_slice[places] = numpy.nan

                    else:
                        pass
                    
        nanwhere = numpy.where(numpy.isnan(roi_mask))
        notnanwhere = numpy.where(numpy.isfinite(roi_mask))
        roi_mask[nanwhere] = background
        roi_mask[notnanwhere] = foreground

        return roi_mask    

def h_goodness_of_fit(data,infodict, indicator = 'rsquared'):
    
    if indicator == 'rsquared':
        ss_err=(infodict['fvec']**2).sum()
        ss_tot=((data-data.mean())**2).sum()
        rsquared=1-(ss_err/ss_tot)
               
        return numpy.abs(rsquared)
        
    else:
        print ('There is no code to produce that indicator. Do it first.')
        raise Exception

def h_generate_VTC(scn_to_analyse=None, 
                   bbox = None, 
                   pdata_num = 0, 
                   **kwargs):

    scan_object = sarpy.Scan(scn_to_analyse)

    # Normalize data
    ndata = sarpy.fmoosvi.analysis.h_normalize_dce(scn_to_analyse)
    
    #Get useful params        
    x_size = ndata.shape[0]
    y_size = ndata.shape[1]
    num_slices = ndata.shape[2]
    reps = ndata.shape[3]
    
    mask = numpy.empty([x_size,y_size,num_slices,reps])
    mask[:] = numpy.nan

    # Deal with bounding boxes

    if bbox is None:        
        bbox = numpy.array([0,x_size-1,0,y_size-1])

    else:      
        bbox = sarpy.fmoosvi.getters.convert_bbox(scn_to_analyse,bbox)

    # Set bounding boxes and get ready to join together
    mask[bbox[0]:bbox[1],bbox[2]:bbox[3],:,:] = 1
        
    #ndata[:,:,:,-1] = numpy.nan
    
    ndata = mask * ndata
    # Reshape it  to stitch together all the data
    nrdata = numpy.empty([x_size,y_size*reps,num_slices])
    
    for s in xrange(num_slices):
        nrdata[:,:,s] = ndata[:,:,s,:].reshape([x_size,y_size*reps])
        
    return {'':nrdata}

## Not working, or of unknown reliability

    
# def h_calculate_KBS(scan_object):
    
#     KBS = 71.16*2
#     # print('Still working on it')
    
#     return KBS

# def h_BS_B1map(zero_BSminus, zero_BSplus, high_BSminus, high_BSplus, scan_with_POI):
    
#     try:
#         TPQQ_POI = scan_with_POI.method.ExcPulse[3] #11.3493504066491 # 5.00591 #11.3493504066491
#         pulse_width_POI = scan_with_POI.method.ExcPulse[0]*1E-3
#     except:
#         print('Please use a scan that has a valid power level for the pulse \
#                 of interest. scan_with_POI.method.ExcPulse[3]')
    
#     TPQQ_BS = high_BSminus.method.BSPulse[3]
    
#     integral_ratio = high_BSminus.method.ExcPulse[10] #0.071941 default from AY

#     #TODO: Write function to calculateKBS
#     KBS = h_calculate_KBS(high_BSminus)
#     gamma = 267.513e6
    
#     # Get phase data from fid
#     offset = h_phase_from_fid(zero_BSplus) - h_phase_from_fid(zero_BSminus)
#     phase_diff = h_phase_from_fid(high_BSplus) - h_phase_from_fid(high_BSminus) + offset
    
#     # Calculate B1 peak
#     B1peak = numpy.sqrt(numpy.absolute(phase_diff)/(2*KBS))
    
#     # Calculate Flip Angle for the pulse of interest
#     alpha_BS = (gamma*B1peak/10000) * (math.pow(10,(TPQQ_BS-TPQQ_POI)/20)) *\
#                 integral_ratio*pulse_width_POI


#     return alpha_BS
    
    



#######