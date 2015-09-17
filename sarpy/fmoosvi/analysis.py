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
                    adata_label=None,
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
    
    # WHY IS THIS HERE!? This makes no sense!
    if adata_label is None:
        # Now calculate the Normalized Intesity voxel by voxel
        norm_data = h_normalize_dce(scn_to_analyse)
    else:
        norm_data = numpy.squeeze(scan_object.adata[adata_label].data)

    # Size info
    x_size = norm_data.shape[0]
    y_size = norm_data.shape[1]
    num_slices = getters.get_num_slices(scn_to_analyse,pdata_num)

    # Now calculate the actual AUC
    auc_data = numpy.squeeze(numpy.empty([x_size,y_size,num_slices]))
    
    if num_slices > 1:
        for slc in range(num_slices):
            if inj_point+auc_reps <= reps:
                auc_data[:,:,slc] = scipy.integrate.simps(norm_data[:,:,slc,inj_point:inj_point+auc_reps],x=time_points)
            else:
                auc_data[:,:,slc] = numpy.nan

    else:
        if inj_point+auc_reps <= reps:
            auc_data[:,:] = scipy.integrate.simps(norm_data[:,:,inj_point:inj_point+auc_reps],x=time_points)
        else:
            auc_data[:,:] = numpy.nan        
       
    # Deal with bounding boxes

    try:
        bbox = scan_object.adata['bbox'].data
    except KeyError:       
        bbox = numpy.array([0,x_size-1,0,y_size-1])    

    if bbox.shape == (4,):            
    
        bbox_mask = numpy.empty([x_size,y_size])
        bbox_mask[:] = numpy.nan        
        bbox_mask[bbox[0]:bbox[1],bbox[2]:bbox[3]] = 1
    
        # First tile for slice
        if num_slices > 1:
            bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size,1),num_slices) 
        else:
            # Do nothing because the bbox max is of one shape
            bbox_mask = numpy.squeeze(bbox_mask.reshape(x_size,y_size,1))

    return {'':auc_data*bbox_mask}


def h_normalize_dce(scn_to_analyse=None, pdata_num = 0):

    scan_object = sarpy.Scan(scn_to_analyse)

    ########### Getting and defining parameters
    
    # Data
    rdata = scan_object.pdata[pdata_num].data

    x_size = rdata.shape[0]
    y_size = rdata.shape[1]
    num_slices = getters.get_num_slices(scn_to_analyse,pdata_num)
    
    # Method params
    #TODO: change this so it doesn't require method file WIHOUT BREAKING IT!
    reps =  scan_object.method.PVM_NRepetitions
    
    if reps != scan_object.pdata[pdata_num].data.shape[-1]:
        reps = scan_object.pdata[pdata_num].data.shape[-1]
        print('\n \n ***** Warning **** \n \n !!! Incomplete dce data for {0} \n \n'.format(scan_object.shortdirname) )

    # Add a third spatial dimension if it's missing.
    if numpy.size(rdata.shape) == 3:

        # add an empty dimension to make it 4D, this code appends the exta axis
        data = sarpy.ImageProcessing.resample_onto.atleast_4d(rdata.copy()) 

        # Move the appended dimension to position 2 to keep data formats the same
        data = data.reshape([data.shape[0], data.shape[1],
                             data.shape[3], data.shape[2]])
    else:
        data = rdata

    # Deal with bounding boxes
    ## Check for bbox traits and create bbox_mask to output only partial data

    try:
        bbox = scan_object.adata['bbox'].data
    except KeyError:       
        bbox = numpy.array([0,x_size-1,0,y_size-1])   

    if bbox.shape == (4,):            
    
        bbox_mask = numpy.empty([x_size,y_size])
        bbox_mask[:] = numpy.nan        
        bbox_mask[bbox[0]:bbox[1],bbox[2]:bbox[3]] = 1
    
        if num_slices > 1:
            # First tile for slice
            bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size,1),num_slices)
            # Next tile for reps
            bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size,num_slices,1),reps)
        else:
            # First tile for slice
            bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size),num_slices)
            # Next tile for reps
            bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size,1),reps)            


    # Calculated params      
    inj_point = sarpy.fmoosvi.analysis.h_inj_point(scn_to_analyse)

    if num_slices > 1:
        norm_data = numpy.empty([x_size,y_size,num_slices,reps])

        for slc in range(num_slices):
            baseline = numpy.mean(data[:,:,slc,0:inj_point],axis=2)
            norm_data[:,:,slc,:] = (data[:,:,slc,:] / numpy.tile(baseline.reshape(x_size,y_size,1),reps))-1

        return norm_data*bbox_mask
    else:
        norm_data = numpy.empty([x_size,y_size,reps])

        baseline = numpy.mean(data[:,:,0,0:inj_point],axis=2)
        norm_data[:,:,:] = (data[:,:,0,:] / numpy.squeeze(numpy.tile(baseline.reshape(x_size,y_size,1),reps)) )-1

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
    
    for slc in range(num_slices):
        
        enhancement_curve[slc,0,:] = time
        enhancement_curve[slc,1,:] = scipy.stats.nanmean(scipy.stats.nanmean(masked_data[:,:,slc,:], axis=0), axis =0)   

    return {'':enhancement_curve}


def h_inj_point(scn_to_analyse=None, pdata_num = 0):

    import sarpy
    import sarpy.ImageProcessing.resample_onto
    from collections import Counter    
    import copy
    import sarpy.helpers

    scan_object = sarpy.Scan(scn_to_analyse)

    # Method params    
    num_slices = sarpy.fmoosvi.getters.get_num_slices(scn_to_analyse,pdata_num)

     # Data
    rawdata = scan_object.pdata[pdata_num].data

    if numpy.size(rawdata.shape) == 3:

        # add an empty dimension to make it 4D, this code appends the exta axis
        data = sarpy.ImageProcessing.resample_onto.atleast_4d(rawdata.copy())

        # Move the appended dimension to position 2 to keep data formats the same
        data = data.reshape([data.shape[0], data.shape[1],
                             data.shape[3], data.shape[2]])
    else:
        data = rawdata

    try:
        # Pool all the data together

        dcelimit = 150

        if dcelimit > rawdata.shape[-1]:
            dcelimit = rawdata.shape[-1]

        global_sum = numpy.sum(numpy.sum(data[:,:,:,0:dcelimit],axis=0),axis=0)

    except IndexError:
        # If it's only one slice, maybe this will still work
        global_sum = data[:,:,:].sum(0).sum(0)

    injection_point = []

    for slc in xrange(num_slices):

        try:
            filtered_sum = smooth_SG(global_sum[slc,:],11,3)
        except IndexError:
            filtered_sum = smooth_SG(global_sum,11,3)

        diff_slice = numpy.diff(filtered_sum)
        std_slice =  numpy.std(diff_slice)

        try:
            injection_point.append(next(i for i,v in enumerate(diff_slice) if v > 2*std_slice))
        except StopIteration:
            print "Could not find the injection point, possibly okay" + str(slc)
            injection_point.append(0)

    # look through the list of elements in injection point and report the most common (mode)
    injection_point_counter = Counter(injection_point)
    final_injection_point = injection_point_counter.most_common(1)

    if scn_to_analyse == 'SARgalbpiB.PG1/6':
        final_injection_point[0] = (10,0)

    return final_injection_point[0][0]+1

def h_calculate_AUGC(scn_to_analyse=None, 
                     adata_label=None,
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
    num_slices = getters.get_num_slices(scn_to_analyse,pdata_num)
    
    # Now calculate the actual AUGC
    augc_data = numpy.empty([x_size,y_size,num_slices])
    
    for slc in range(num_slices):
        augc_data[:,:,slc] = scipy.integrate.simps(data[:,:,slc,inj_point:inj_point+augc_reps],x=time_points)
    
    # Deal with bounding boxes

    try:
        bbox = scan_object.adata['bbox'].data
    except KeyError:       
        bbox = numpy.array([0,x_size-1,0,y_size-1])   
          
    if bbox.shape == (4,):            
    
        bbox_mask = numpy.empty([x_size,y_size])
        bbox_mask[:] = numpy.nan        
        bbox_mask[bbox[0]:bbox[1],bbox[2]:bbox[3]] = 1
    
        # First tile for slice
        bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size,1),num_slices)

        # If this gives a value error about operands not being broadcast together, go back and change your adata to make sure it is squeezed
    return {'':augc_data*bbox_mask}
    
def h_conc_from_signal(scn_to_analyse=None, 
                       other_scan_name=None, 
                       adata_label = None,
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
    try:
        bbox = scan_object.adata['bbox'].data
    except KeyError:       
        bbox = numpy.array([0,x_size-1,0,y_size-1])   

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
            for slc in xrange(num_slices):

                baseline_s = numpy.mean(source_data[x,y,slc,0:inj_point])
                E1 = numpy.exp(-TR/data_t1map[x,y,slc])
                c = numpy.cos(numpy.radians(FA))
                T1[x,y,slc,0:inj_point] = data_t1map[x,y,slc]

                # Use the SPGR equation twice, once before agent admin. and once after. If you
                # divide the two, the M0s cancel out, and so do the sin thetas. what remains is:
                # (1-E1)*(1-E0*cos) / ((1-E0) * (1-E1*cos)). use wolfram alpha to solve this:
                # http://www.wolframalpha.com/input/?i=solve+%28%281-x%29*%281-b*c%29%29+%2F+%28%281-b%29*%281-x*c%29%29+%3D+r+for+x

                for rep in xrange(inj_point,source_data.shape[-1]):                    
                    s = source_data[x,y,slc,rep] / baseline_s
                    E2 = (-E1*c + E1*s - s + 1) / (E1*s*c - E1*c - s*c +1)
                    T1[x,y,slc,rep] = -TR / numpy.log(E2)


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

    ######################
    # WARNING
    # This will likely break in v4 of sarpy where new masterlists were introduced, and bbox behaviour was changed
    # Guess I shouldn't have broken that cardinal rule afterall
    ######################
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
            for slc in xrange(dataShape[2]):

                dA = scanA.pdata[pdata_num].data[x,y,slc,:]
                dB = scanB.pdata[pdata_num].data[x,y,slc,:]

                gd_conc_stitched[x,y,slc,:] = numpy.hstack((dA,dB))

    return {'':gd_conc_stitched,
            '_time': timeAB,
            '_fulltime': timeAtotB,
            '_idx':idx}
    
### T1 fitting section
    
##############
    #
    # Fitting section
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

def T1eff_to_T1(T1,     # true T1
                Td,     # inv pulse to first excitation pulse
                Tp,     # last excitation to next inversion pulse
                tau,    # time of one frame (long is better for long T1)
                T1_eff, # T1*
                b,      # M(0)/M(infty)
                c       # M(N-1)/M(infty)
               ):    
    ''' re-arranged equation [6] from Chuang and Koretsky 'Improved MEMRI Tract
    Tracing by Fast T1 Mapping' essentially a LL analysis method without the
    need for knowing FAmaps '''
    return ((1-2*numpy.exp(-Td/T1) + 
            numpy.exp(-(Tp+Td)/T1)) * ((1-numpy.exp(-tau/T1_eff)) /
                                       (1-numpy.exp(-tau/T1))) - 
            c*(numpy.exp(-tau/T1_eff)/(numpy.exp(-tau/T1))) * numpy.exp(-(Tp+Td)/T1)
            - b) 

def T1simple(T1,     # true T1
             tau,    # time of one frame (long is better for long T1)
             T1_eff, # T1*
             FA      # Flip angle
             ):    
    return 1.0/T1 - 1.0/T1_eff - numpy.divide(numpy.log(numpy.cos(FA)),tau) 
             
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

def h_fitpx_T1_LL_FAassumed(scn_to_analyse=None,
                        y_data=None,
                        slc=None,
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
    delay = float(scan_object.acqp.ACQ_inversion_time[0]) # same as .method.PVM_InversionTime
    FA = float(numpy.radians(scan_object.acqp.ACQ_flip_angle))
    interSliceDelay = float(scan_object.method.InterSliceDelay)
    sliceAcqTime = float(scan_object.method.SliceSeqTime)
    sliceOrder = scan_object.method.PVM_ObjOrderList

    # The following are slice dependent parameters
    slicedelay = (interSliceDelay + sliceAcqTime) * sliceOrder.index(slc)
    #t_data = delay + numpy.arange(0,Nframes)*tau + slicedelay 
    t_data = numpy.arange(0,Nframes)*tau
    
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
            #print x,y,slc
            raise
    else:
        print('fit type {0} not implemented yet'.format(fit_algorithm))

    [a,b,T1_eff,phi] = fit_params

    # Step 2: Calculate T1 from T1_eff using equation relating T1, T1_eff and cos(a)

    calc_params = (tau,T1_eff,FA)

    try:                    
        # Using Bisect method (other possibilities include newton, brentq:
        T1 = scipy.optimize.bisect(T1simple,5,4000,
                                   args = (calc_params))

    except ValueError:
        T1=numpy.nan
        pass
    
    # Make absurd values nans to make my life easer:
    if (T1<0 or T1>=1e4):
        T1 = numpy.nan

    return infodict,mesg,ier,fit_params, T1,t_data

def h_fit_T1_LL_FAassumed(scn_to_analyse=None,
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
    t_data_arr = numpy.empty_like(data_after_fitting,dtype='object')

    infodict1 = numpy.empty_like(data_after_fitting,dtype=dict)
    mesg1 = numpy.empty_like(data_after_fitting,dtype=str)
    ier1 = numpy.empty_like(data_after_fitting)
    goodness_of_fit1 = numpy.empty_like(data_after_fitting)

    ## Check for bbox traits and create bbox_mask to output only partial data
    try:
        bbox = scan_object.adata['bbox'].data
    except KeyError:       
        bbox = numpy.array([0,x_size-1,0,y_size-1])   
        
    # Start the fitting process        

    for slc in xrange(num_slices):          
        for x in xrange(bbox[0],bbox[1]):
            for y in xrange(bbox[2],bbox[3]):

                # Deal with scans that have only one slice
                if num_slices > 1:
                    y_data = data[x,y,slc,:]
                else:
                    y_data = data[x,y,:]

                [infodict,mesg,ier,res_fit_params,T1,t_data] = h_fitpx_T1_LL_FAassumed(scn_to_analyse,
                                                                                   y_data=y_data,
                                                                                   slc=slc,
                                                                                   fit_algorithm=fit_algorithm)

                [a,b,T1_eff,phi] = res_fit_params

                a_arr[x,y,slc] = a
                b_arr[x,y,slc] = b
                T1_eff_arr[x,y,slc] = T1_eff
                phi_arr[x,y,slc] = phi
                T1_arr[x,y,slc] = T1
                t_data_arr[x,y,slc] = t_data

                data_after_fitting[x,y,slc] = T1

                if fit_algorithm == 'leastsq':   
                    infodict1[x,y,slc] = infodict
                    mesg1[x,y,slc] = mesg
                    goodness_of_fit1[x,y,slc] = h_goodness_of_fit(y_data,infodict)  

                ier1[x,y,slc] = ier

## Moved the entire T1 fitting process into a separate function

    fit_params1 = {'a': a_arr,
                   'b': b_arr,
                   'T1_eff': T1_eff_arr,
                   'phi': phi_arr}

    if fit_algorithm == 'leastsq':

        result = {'rawdata':data,
                  't_data':t_data_arr,
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

def h_fitpx_T1_LL_FAind(scn_to_analyse=None,
                        y_data=None,
                        slc=None,
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
    delay = float(scan_object.acqp.ACQ_inversion_time[0]) # same as .method.PVM_InversionTime
    FA = float(numpy.radians(scan_object.acqp.ACQ_flip_angle))
    interSliceDelay = float(scan_object.method.InterSliceDelay)
    sliceAcqTime = float(scan_object.method.SliceSeqTime)
    sliceOrder = scan_object.method.PVM_ObjOrderList

    # The following are slice dependent parameters
    slicedelay = (interSliceDelay + sliceAcqTime) * sliceOrder.index(slc)
    #t_data = delay + numpy.arange(0,Nframes)*tau + slicedelay
    t_data = numpy.arange(0,Nframes)*tau

    Td = delay + slicedelay
    Tp = total_TR - (Nframes -1.0)*tau - Td    

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
            #print x,y,slc
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
    t_data_arr = numpy.empty_like(data_after_fitting,dtype='object')

    infodict1 = numpy.empty_like(data_after_fitting,dtype=dict)
    mesg1 = numpy.empty_like(data_after_fitting,dtype=str)
    ier1 = numpy.empty_like(data_after_fitting)
    goodness_of_fit1 = numpy.empty_like(data_after_fitting)

    ## Check for bbox traits and create bbox_mask to output only partial data
    try:
        bbox = scan_object.adata['bbox'].data
    except KeyError:       
        bbox = numpy.array([0,x_size-1,0,y_size-1])   
        
    # Start the fitting process        

    for slc in xrange(num_slices):          
        for x in xrange(bbox[0],bbox[1]):
            for y in xrange(bbox[2],bbox[3]):

                # Deal with scans that have only one slice
                if num_slices > 1:
                    y_data = data[x,y,slc,:]
                else:
                    y_data = data[x,y,:]

                [infodict,mesg,ier,res_fit_params,T1,t_data] = h_fitpx_T1_LL_FAind(scn_to_analyse,
                                                                                   y_data=y_data,
                                                                                   slc=slc,
                                                                                   fit_algorithm=fit_algorithm)

                [a,b,T1_eff,phi] = res_fit_params

                a_arr[x,y,slc] = a
                b_arr[x,y,slc] = b
                T1_eff_arr[x,y,slc] = T1_eff
                phi_arr[x,y,slc] = phi
                T1_arr[x,y,slc] = T1
                t_data_arr[x,y,slc] = t_data
                test = numpy.ndarray([10,10,10],dtype='object')

                data_after_fitting[x,y,slc] = T1

                if fit_algorithm == 'leastsq':   
                    infodict1[x,y,slc] = infodict
                    mesg1[x,y,slc] = mesg
                    goodness_of_fit1[x,y,slc] = h_goodness_of_fit(y_data,infodict)  

                ier1[x,y,slc] = ier

## Moved the entire T1 fitting process into a separate function

    fit_params1 = {'a': a_arr,
                   'b': b_arr,
                   'T1_eff': T1_eff_arr,
                   'phi': phi_arr}

    if fit_algorithm == 'leastsq':

        result = {'rawdata':data,
                  't_data':t_data_arr,
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
        for slc in xrange(roi_mask.shape[2]):
        
            curr_slice = roi_mask[:,:,slc]
            
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
                          
        for slc in xrange(roi_mask.shape[2]):
            
            for x in xrange(peaks):
        
                curr_slice = roi_mask[:,:,slc]
              
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
                   vtc_type = None,
                   roi_label = None,
                   pdata_num = 0, 
                   **kwargs):

    scan_object = sarpy.Scan(scn_to_analyse)

    # Hack for the moment

    ndata = scan_object.pdata[0].data

    #Get useful params        
    x_size = ndata.shape[0]
    y_size = ndata.shape[1]
    num_slices = getters.get_num_slices(scn_to_analyse,pdata_num)
    reps = ndata.shape[-1]

    # Deal with bounding boxes
    try:
        bbox = scan_object.adata['bbox'].data
    except KeyError:       
        bbox = numpy.array([0,x_size-1,0,y_size-1])    

    # Specify VTC data

    if vtc_type is None:
        ndata = sarpy.fmoosvi.analysis.h_normalize_dce(scn_to_analyse)
    elif vtc_type == 'gd_KM':
        ndata = numpy.squeeze(sarpy.Scan(scn_to_analyse).adata['gd_KM'].data)
    elif vtc_type =='CEST':

        import cest.analysis as cest
        reload(cest)

        # Get the ndata to figure out the shape of the array (since some data gets thrown out)
        # make it nan
        throwaway, ndata = cest.cest_spectrum(scn_to_analyse,
                                    bbox[0],bbox[2],
                                    shift_water_peak = True,
                                    normalize=True,
                                    normalize_to_ppm = 200,
                                    ppm_limit_min = -5,
                                    ppm_limit_max = 5,
                                    exclude_ppm = 66.6,
                                    pdata_num = 0)

        reps = len(ndata)

        ndata = numpy.empty(shape=[x_size,y_size,reps]) + numpy.nan

        for xx in xrange(bbox[0],bbox[1]):
            for yy in xrange(bbox[2],bbox[3]):

                xvals, res = cest.cest_spectrum(scn_to_analyse,
                                                            xx,yy,
                                                            shift_water_peak = True,
                                                            normalize=True,
                                                            normalize_to_ppm = 200,
                                                            ppm_limit_min = -5,
                                                            ppm_limit_max = 5,
                                                            exclude_ppm = 66.6,
                                                            pdata_num = 0)

                # Flip the x-axis to deal with stupid spectroscopy convention of + then -
                ndata[xx,yy,:] = res[::-1]

    if roi_label is None:
        mask = numpy.squeeze(numpy.empty([x_size,y_size,num_slices,reps]) + numpy.nan)
    else:
        mask = scan_object.adata[roi_label].data
        mask = numpy.tile(mask.reshape(x_size,y_size,num_slices,1),reps)

    # Set bounding boxes and get ready to join together
    if num_slices > 1:
        mask[bbox[0]:bbox[1],bbox[2]:bbox[3],:,:] = 1
        ndata = mask * ndata

        if roi_label is None: # In case you want a bboxed VTC map
            mask[bbox[0]:bbox[1],bbox[2]:bbox[3],:,:] = 1
        else: # In the case where you want JUST the roi VTCs (not bboxed)
            mask = numpy.where(mask<1,numpy.nan,1)
            mask = numpy.squeeze(mask)        

        # Reshape it  to stitch together all the data
        nrdata = numpy.empty([x_size,y_size*reps,num_slices]).astype(float)  
        for s in xrange(num_slices):
            tmp = ndata[:,:,s,:].flatten()
            tmp[::reps] = numpy.nan
            # Sets the last point of each enhancement curve to numpy.nan
            nrdata[:,:,s] = tmp.reshape([x_size,y_size*reps])
    
    else: # Account for single slice data

        if roi_label is None: # In case you want a bboxed VTC map
            mask[bbox[0]:bbox[1],bbox[2]:bbox[3],:] = 1
        else: # In the case where you want JUST the roi VTCs (not bboxed)
            mask = numpy.where(mask<1,numpy.nan,1)
            mask = numpy.squeeze(mask)

        ndata = mask * ndata
        nrdata = numpy.empty([x_size,y_size*reps]) 
        # Incomprehensible list comprehenshion
        # Sets the last point of each enhancement curve to numpy.nan 
        tmp = ndata.flatten()
        tmp[::reps] = numpy.nan   

        nrdata = tmp.reshape([x_size,y_size*reps])
    return {'':nrdata}

def createSaveVTC(scn_to_analyse=None,
                  adata_label=None,
                  roi_label = 'roi',
                  pdata_num = 0, 
                  **kwargs):

    import pylab

    scan_object = sarpy.Scan(scn_to_analyse)
    fig = pylab.figure()
    G = pylab.matplotlib.gridspec.GridSpec(1,1, wspace=0.0, hspace=0.0)   
    reps = scan_object.pdata[0].data.shape[-1]
    dat = scan_object.pdata[0].data

    x_size = dat.shape[0]
    y_size = dat.shape[1]

    # Deal with bounding boxes
    try:
        bbox = scan_object.adata['bbox'].data
    except KeyError:       
        bbox = numpy.array([0,x_size-1,0,y_size-1])   


    imgdata = numpy.mean(dat,axis=2)
    vtcdata = scan_object.adata[adata_label].data

    axs = fig.add_subplot(G[0, 0])
    #aspect = scn.method.PVM_SpatResol[0]*(bbox[1]-bbox[0]) / scn.method.PVM_SpatResol[1]*(bbox[3]-bbox[2])
    aspect= (scan_object.method.PVM_FovCm[0]/scan_object.method.PVM_Matrix[0])/ \
            (scan_object.method.PVM_FovCm[1]/scan_object.method.PVM_Matrix[1])
    #aspect
    axs.imshow(imgdata[bbox[0]:bbox[1],\
                       bbox[2]:bbox[3]],\
                       cmap='gray', 
                       interpolation='None',
                       alpha=1.0,
                       zorder=0,
                       aspect=aspect,
                       vmin=0,
                       vmax=1)
    axs.set_axis_off()
    fig.canvas.draw()
    pylab.axis('off')

    #axs=fig.add_subplot(G[0, 0])
    box = axs.get_position().bounds
    height = box[3] / (bbox[1]-bbox[0])

    for ht,i in enumerate(xrange(bbox[0], bbox[1])):

        #, simply use the add_axes() method which takes a list of 
        # [left, bottom, width, height] values in 0-1 relative figure coordinates
        
        #The add_axes method takes a list of four values, which are 
        #xmin, ymin, dx, and dy for the subplot, where xmin and ymin 
        #are the coordinates of the lower left corner of the subplot, 
        #and dx and dy are the width and height of the subplot, 
        #with all values specified in relative units 
        #(where 0 is left/bottom and 1 is top/right)

        tmpax = fig.add_axes([box[0], box[1]+ht*height,
                             box[2], height])


        tmpax.plot(vtcdata[i,((bbox[2])*reps):((bbox[3])*reps)],
                           color='g', 
                           linewidth=.1,
                           zorder=1)
        tmpax.set_axis_off()
        #pylab.ylim([35000,65000])
        pylab.xlim([0,((bbox[3])*reps)-(bbox[2])*reps])
    
    pylab.savefig('{0}.png'.format(scan_object.shortdirname.split('/')[0]),dpi=600)
    
    pylab.close(fig)
    fig = pylab.figure()
    G = pylab.matplotlib.gridspec.GridSpec(1,1, wspace=0.0, hspace=0.0)
    axs=fig.add_subplot(G[0, 0])
    axs.imshow(numpy.flipud(imgdata[bbox[0]:bbox[1],\
                           bbox[2]:bbox[3]]),\
                           cmap='gray', 
                           interpolation='None',
                           alpha=1.0,
                           zorder=0,
                           aspect=aspect)  
    axs.set_axis_off()    
    pylab.savefig('{0}.png'.format(scan_object.shortdirname.split('/')[0]+'_anatomy'),dpi=300)    

######################
###
###
### Fitting T1 IR data
###
###
#####################

def h_func_T1_IR(params,TI,parallelExperiment=False):

    a,b,T1 = params
    # Function to fit:
    # Intensity = a*(1-b*numpy.exp(-TI / T1))

    res = a*(1-b*numpy.exp(-TI/T1))

    if numpy.isfinite(res).all():
    	# in case you're fitting pdata, you'll need to take abs

    	if parallelExperiment:
            return numpy.abs(res)
    	else:
        	return res
    else:
        return [1e30]*TI.shape[-1]

def h_residual_T1_IR(params, y_data, TI, parallelExperiment=False):
    
    bounds = numpy.zeros(shape=[3,2])
    bounds[0,0] = 1.E2
    bounds[0,1] = 5.E10
    bounds[1,0] = 0.1
    bounds[1,1] = 3.0
    bounds[2,0] = 0.0
    bounds[2,1] = 4000.0

    if h_within_bounds(params,bounds):
        return numpy.abs(y_data - h_func_T1_IR(params,TI,parallelExperiment))
    else:
        return 1e30

def h_fitpx_T1_IR(y_data,TI,initial_guess = None,parallelExperiment=False):

    if initial_guess is None:

        # Get the TI when Y-data is at minimum
        guess = TI[list(numpy.real(y_data)).index(min(numpy.real(y_data),key=abs))]
        params = [max(y_data),2.0,guess]
    else:
        params = initial_guess

    # Not sure if I need this, but I needed to resolve this error:
    # TypeError: Cannot cast array data from dtype('complex128') to dtype('float64') 
    # according to the rule 'safe'

    fit_params,cov,infodict,mesg,ier = scipy.optimize.leastsq(
                                                        h_residual_T1_IR,
                                                        params,
                                                        args=(y_data, TI,parallelExperiment), 
                                                        full_output = True,
                                                        maxfev = 600)

    return fit_params,cov,infodict,mesg,ier


def h_fit_T1_IR(scan_name_list,parallelExperiment = False):

    # Convert scan_name_list to scan_list:

    scan_list = []
    for hh in scan_name_list:
        scan_list.append(sarpy.Scan(hh))

    data_shape = scan_list[0].fftfid.shape
    data_shape = numpy.atleast_3d(data_shape)[0]


    # Collect Inversion Times
    TI = numpy.array([float(x.method.PVM_SelIrInvTime) for x in scan_list])

    # Check to make sure all the scans have the same data
    assert [x.pdata[0].data.shape == data_shape for x in scan_list], "Data shapes don't match"

    # Definitions
    data_after_fitting = numpy.nan + numpy.empty( data_shape )
    fit_params1 = numpy.empty_like(data_after_fitting,dtype=dict)
    a_arr = numpy.empty_like(data_after_fitting)+numpy.nan
    b_arr = numpy.empty_like(data_after_fitting)+numpy.nan
    T1_arr = numpy.empty_like(data_after_fitting)+numpy.nan    
    raw_data = numpy.empty_like(data_after_fitting,dtype='object') 

    infodict1 = numpy.empty_like(data_after_fitting,dtype=dict)
    mesg1 = numpy.empty_like(data_after_fitting,dtype=str)
    ier1 = numpy.empty_like(data_after_fitting)
    goodness_of_fit1 = numpy.empty_like(data_after_fitting)

    # Loop through the data

    # Hack to get it to work, fix this properly later
    if len(data_shape) == 2:
        hack_shape = 1
    else:
        hack_shape = data_shape[2]

    for x in xrange(data_shape[0]):
        for y in xrange(data_shape[1]):
            for slc in xrange(hack_shape):

                # Collect the data
                if hack_shape> 2:
                    if parallelExperiment:
                        # in case of a Parallel Experiment as in the case of 
                        # patient DevRF1106.jr1 - for RF1106 experiment          
                        y_data = numpy.array([sc.pdata[0].data[x,y,slc] for sc in scan_list])                                  	
                    else:
                        y_data = numpy.array([sc.fftfid[x,y,slc] for sc in scan_list])

                else:
                    if parallelExperiment:
                        # in case of a Parallel Experiment as in the case of 
                        # patient DevRF1106.jr1 - for RF1106 experiment  
                        y_data = numpy.array([sc.pdata[0].data[x,y] for sc in scan_list])                                          	
                    else:
                        y_data = numpy.array([sc.fftfid[x,y] for sc in scan_list])

                y_data = numpy.real(y_data)

            	fit_params,cov,infodict,mesg,ier = h_fitpx_T1_IR(numpy.squeeze(y_data),
                                                                 TI,
                                                                 parallelExperiment=parallelExperiment)

                [a,b,T1] = fit_params

                # Make absurd values nans to make my life easer:
                if (T1<0 or T1>=1e4):
                    T1 = numpy.nan

                # Hacked to get it working with one slice data
                if hack_shape > 2:
                    
                    data_after_fitting[x,y,slc]  = T1

                    a_arr[x,y,slc] = a
                    b_arr[x,y,slc] = b
                    T1_arr[x,y,slc] = T1
                    raw_data[x,y,slc] = y_data         

                    infodict1[x,y,slc] = infodict
                    mesg1[x,y,slc] = mesg
                    goodness_of_fit1[x,y,slc] = h_goodness_of_fit(y_data,infodict)  
                    ier1[x,y,slc] = ier

                else:

                    data_after_fitting[x,y]  = T1

                    a_arr[x,y] = a
                    b_arr[x,y] = b
                    T1_arr[x,y] = T1       
                    raw_data[x,y] = y_data         

                    infodict1[x,y] = infodict
                    mesg1[x,y] = mesg
                    goodness_of_fit1[x,y] = h_goodness_of_fit(y_data,infodict)  
                    ier1[x,y] = ier


                fit_params1 = {'a': a_arr,
                               'b': b_arr,
                               'T1': T1_arr}                

    # Returning the fit results
    result = {'rawdata':raw_data,
              'TI':TI,
              'params':fit_params1,
              'infodict':infodict1,
              'ier':ier1,
              'goodness':goodness_of_fit1,
              'mesg':mesg1}

    return {'':numpy.squeeze(data_after_fitting),
            '_fit':result}         

def checkSNR(scn_to_analyse=None,
             pdata_num = 0,
             limits = (0,100)):

    import pylab

    scan_object = sarpy.Scan(scn_to_analyse)

    dcedata = scan_object.pdata[pdata_num].data

    assert(dcedata.shape[-1] >3),"Need more than 3 reps to check the SNR"

    # First check to see if there is an injection point, if so consider the points before it only:
    inj_point = h_inj_point(scn_to_analyse)

    if  inj_point > 1:
        snrcheck_data = dcedata[:,...,1:inj_point]
        print('Injection, points before {0} used'.format(inj_point))

    else:
        snrcheck_data = dcedata[:,:,:,1:]

    # Mean of each pixel divided by the std dev of each pixel
    # SAR's special SNR definition
    snrMap = numpy.divide(numpy.mean(snrcheck_data,axis=-1),
                          numpy.std(snrcheck_data,axis=-1))

    # Conventional SNR definition
    # Take the mean of the top snrMap voxels in the map above, and divide them by the sqrt(std in the bottom)
    hist = numpy.histogram(snrMap)

    lowEstimate = numpy.mean(snrcheck_data[snrMap>hist[1][1]])/numpy.sqrt(numpy.std(snrcheck_data[snrMap<hist[1][1]]))
    highEstimate= numpy.mean(snrcheck_data[snrMap>hist[1][4]])/numpy.sqrt(numpy.std(snrcheck_data[snrMap<hist[1][1]]))


    gridSize = numpy.ceil(numpy.sqrt(dcedata.shape[-2]))

    pylab.suptitle('Estimated Conventional SNR: {0} - {1}'.format(int(lowEstimate),int(highEstimate)))


    if len(dcedata.shape)>3:
        for s in xrange(dcedata.shape[-2]):
            pylab.subplot(gridSize,gridSize,s+1)
            pylab.imshow(snrMap[:,:,s],vmin=limits[0],vmax=limits[1])
            pylab.title('Slice {0}'.format(s))
            pylab.colorbar()
            pylab.axis('off')
    else:
            pylab.imshow(snrMap,vmin=limits[0],vmax=limits[1])
            pylab.colorbar()
            pylab.axis('off')

    return snrMap

def MMCA_model_fit(timecourse, start_fit=0, dt=1):
    '''Fits a linear model to a time course starting some way through the injection.
    Ostensibly, the jump from baseline to the data just after injection is related to 
    blood volume (100% BV should jump to a concentration of v_c). The slope of the curve
    after that initial line is related to extravasation (as a result of leakage from vessels?)'''

    # Removing the nans before calculating the BV and leakage
    timeaxis = dt*numpy.arange(len(timecourse))
    start_time = timeaxis[start_fit]

    nan_whereabout = numpy.where(numpy.isfinite(timecourse))
    timecourse = timecourse[nan_whereabout]
    timeaxis = timeaxis[nan_whereabout]
    timecourse = timecourse[numpy.where(timeaxis>=start_time)]
    timeaxis = timeaxis[numpy.where(timeaxis>=start_time)]
    

    try: 
        (slope,intercept,r,tt,stderr) = scipy.stats.linregress(timeaxis, timecourse)
        rel_bloodvol = intercept
        leakage_param = slope
    except ValueError:
        rel_bloodvol = numpy.nan
        leakage_param = numpy.nan   

    return (rel_bloodvol,leakage_param)

def h_fit_ec(scn_to_analyse=None,
            roi_label=None,
            adata_label=None,
            pdata_num = 0, 
            **kwargs):
    '''Fits a linear model to a time course starting some way through the injection.
    Ostensibly, the jump from baseline to the data just after injection is related to 
    blood volume (100% BV should jump to a concentration of v_c). The slope of the curve
    after that initial line is related to extravasation (as a result of leakage from vessels?)'''
    
    import scipy.stats
    scan_object = sarpy.Scan(scn_to_analyse)
    

    if adata_label is None:
        norm_data = h_normalize_dce(scn_to_analyse)
        #norm_data = sarpy.Scan(scn_to_analyse).pdata[0].data
    else:
        norm_data = numpy.squeeze(scan_object.adata[adata_label].data)

    # Size info
    x_size = norm_data.shape[0]
    y_size = norm_data.shape[1]
    num_slices = getters.get_num_slices(scn_to_analyse,pdata_num)
    reps =  scan_object.method.PVM_NRepetitions

    if roi_label is None:
        norm_data =norm_data
    else:

        msk = numpy.squeeze(scan_object.adata[roi_label].data)

        if num_slices > 1:
            # First tile for slice
            # Next tile for reps
            bbox_mask = numpy.tile(msk.reshape(x_size,y_size,num_slices,1),reps)
        else:
            # First tile for slice
            bbox_mask = numpy.tile(msk.reshape(x_size,y_size),num_slices)
            # Next tile for reps
            bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size,1),reps)   

        norm_data = norm_data*bbox_mask

    # there are problems with using phase encodes for certain cases (maybe 3D)
    # so now I have to use the tuid time
    total_time = scan_object.method.PVM_ScanTimeStr
    format = "%Hh%Mm%Ss%fms"
    t=datetime.datetime.strptime(total_time,format)
    total_time = (3600*t.hour) + (60*t.minute) + (t.second) + t.microsecond*1E-6
    
    time_per_rep = numpy.round(numpy.divide(total_time,reps))
    inj_point = sarpy.fmoosvi.analysis.h_inj_point(scn_to_analyse)

    print inj_point
    # Deal with bounding boxes
    try:
        bbox = scan_object.adata['bbox'].data
    except KeyError:    
        bbox = numpy.array([0,x_size,0,y_size])
        

    # Now calculate the actual Blood Volume & Leakage Param
    bloodvolume_data = numpy.squeeze(numpy.empty([x_size,y_size,num_slices])+numpy.nan)
    leakage_data = numpy.nan*numpy.squeeze(numpy.empty([x_size,y_size,num_slices])+numpy.nan)
    
    # Start the fitting process        
    for slc in xrange(num_slices):          
        for x in xrange(bbox[0],bbox[1]):
            for y in xrange(bbox[2],bbox[3]):

                if num_slices > 2:
                    bloodvolume, leakage = MMCA_model_fit(norm_data[x,y,slc,:],
                                                          start_fit=inj_point+2,
                                                          dt = time_per_rep)
                    bloodvolume_data[x,y,slc] = bloodvolume
                    leakage_data[x,y,slc] = leakage
                else:
                    bloodvolume, leakage = MMCA_model_fit(norm_data[x,y],
                                                          start_fit=inj_point+2,
                                                          dt = time_per_rep)                    
                    bloodvolume_data[x,y] = bloodvolume
                    leakage_data[x,y] = leakage

    return {'bloodvolume':bloodvolume_data, 
            'leakage': leakage_data}

######################
###
###
### Calculating Bolus Arrival Time (BAT)
###
###
#####################  

def bolus_arrival_time(scn_to_analyse=None,
                       pdata_num = 0,
                       bbox = None,
                       verbose=False,
                       **kwargs):
    
    def runningSum(x, N):
        output = numpy.empty_like(x)
        
        if len(x.shape)>3: #Multi Slice
            for i in xrange(x.shape[0]):
                for j in xrange(x.shape[1]):
                    for k in xrange(x.shape[2]):
                        output[i,j,k,:] = numpy.convolve(x[i,j,k,:],numpy.ones(N))[N-1:]
        elif len(x.shape) ==3: #Single Slice
            for i in xrange(x.shape[0]):
                for j in xrange(x.shape[1]):
                        output[i,j,:] = numpy.convolve(x[i,j,:],numpy.ones(N))[N-1:]            
            
        return output

    def first4D(a):  # Written by SAR to return the first non-zero entry in a 4D array in the 4th axis
        
        if len(a.shape)>3:
            di = numpy.zeros(a.shape[0:3])+numpy.Inf
            for i, j, k, l in zip(*numpy.where(a>0)):
                if l<di[i,j,k]:
                    di[i,j,k] = l
                    
        elif len(a.shape)==3:
            
            di = numpy.zeros(a.shape[0:2])+numpy.Inf
            for i, j, l in zip(*numpy.where(a>0)):
                if l<di[i,j]:
                    di[i,j] = l            
        return di

    import sarpy
    import sarpy.fmoosvi.getters as getters
    import datetime
    
    ## Get the data and various sizes
    
    scan_object = sarpy.Scan(scn_to_analyse)
    num_slices = getters.get_num_slices(scn_to_analyse,pdata_num)
    data = scan_object.pdata[0].data
    x_size = data.shape[0]
    y_size = data.shape[1]
    
    ## Deal with bounding boxes
    if bbox is None:        
        bbox = numpy.array([0,x_size-1,0,y_size-1])    
    else:      
        bbox = sarpy.fmoosvi.getters.convert_bbox(scn_to_analyse,bbox) 

    if bbox.shape == (4,):            

        bbox_mask = numpy.empty([x_size,y_size])
        bbox_mask[:] = numpy.nan        
        bbox_mask[bbox[0]:bbox[1],bbox[2]:bbox[3]] = 1

        # First tile for slice
        bbox_mask = numpy.squeeze(numpy.tile(bbox_mask.reshape(x_size,y_size,1),num_slices))
        
    ############################################################   
    #
    # Begin BAT Code here
    #
    ############################################################
    
    ## Determine the injection point using Firas' function
    # This averages and flattens the entire DCE-MRI dataset to a single 
    # Time course and determines the point of injection in a 
    # very simple way (using the largest change)
    
    #inj_point =  sarpy.fmoosvi.analysis.h_inj_point(scn_to_analyse) + 1
    inj_point =  h_inj_point(scn_to_analyse) + 1

    # Check for number of slices; this could be made more efficient, 
    # there is a lot of repeated code here
    
    # only consider means and std for data at baseline but exclude
    # 1st point (that one is always wonky (TM))
    if num_slices >1:
        sd = numpy.std(data[:,:,:,1:inj_point],axis=-1)*bbox_mask
        mean = numpy.mean(data[:,:,:,1:inj_point],axis=-1)*bbox_mask
        
    elif num_slices ==1:
        sd = numpy.std(data[:,:,1:inj_point],axis=-1)*bbox_mask
        mean = numpy.mean(data[:,:,1:inj_point],axis=-1)*bbox_mask        
        
    ## Compute the SD and Mean to be used for Analysis
    sd = numpy.repeat(sd,data.shape[-1]).reshape(data.shape)
    mean = numpy.repeat(mean,data.shape[-1]).reshape(data.shape)
    
    ## Setting Conditions
    cond1 = numpy.where(data >= (mean+3*sd),1,0)

    # Condition 2
    condaux = numpy.where(data >= mean+2*sd,1,0)
    cond2 = numpy.where(runningSum(condaux,3)>=2,1,0) 

    # Condition 3
    condaux = numpy.where(data >= mean+sd,1,0)
    cond3 = numpy.where(runningSum(condaux,5)>=4,1,0)

    # Condition 4
    condaux = numpy.where(data >= mean,1,0)
    cond4 = numpy.where(runningSum(condaux,8)>=8,1,0)

    # Condition 5 - super condition
    allcond = cond1 + cond2 + cond3 + cond4
    #cond = numpy.where(runningSum(allcond,3) >=3)
    cond = numpy.where(runningSum(allcond,3) >=3,runningSum(allcond,3),0)

    # Find the first finite entry in the array
    BAT = first4D(cond)

    # Condition 6 - compare cond to inj_point 
    # (throws away if it's more than 10 points away from inj_point)
    #BAT = numpy.where(BAT > inj_point+30,numpy.Inf,BAT)
    BAT = numpy.where(BAT==0,numpy.Inf,BAT)
    
    # Switch from Index to a time in seconds
    
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

    time_per_rep = numpy.divide(total_time,reps)
    
    if verbose: 
        return {'':[time_per_rep*(BAT+1-inj_point), (cond1,cond2,cond3,cond4,1.*cond/4)]}
    else:
        return {'':time_per_rep*(BAT+1-inj_point)}


def h_threshold(data_matrix,
                min_val = 0,
                max_val = 100):


    # Set all points below min_val to be nan
    data_matrix = numpy.where(data_matrix<min_val,numpy.nan,data_matrix)

    # Set all points above max_val to be nan
    data_matrix = numpy.where(data_matrix>max_val,numpy.nan,data_matrix)
    
    return data_matrix

def h_make_binary_mask(data_matrix,
                       min_val = 0,
                       max_val = 100):

    new_data_matrix = h_threshold(data_matrix,min_val,max_val)

    return numpy.where(numpy.isfinite(new_data_matrix),1,new_data_matrix)

def smooth_SG(y, window_size, order, deriv=0, rate=1):

    # Implementation of a Savitzky-Golay filter
    # Acquired from: 
    # http://stackoverflow.com/questions/22988882/how-to-smooth-a-curve-in-python

    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')        