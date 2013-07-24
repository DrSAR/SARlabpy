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
import math
import copy
import os
import json
import nibabel
import datetime
import collections

def h_calculate_AUC(scan_object, bbox = None, time = 60, pdata_num = 0):
    
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
    inj_point = sarpy.fmoosvi.analysis.h_inj_point(scan_object)
    
    # Now calculate the Normalized Intesity voxel by voxel
    norm_data = h_normalize_dce(scan_object)

    # Size info
    x_size = norm_data.shape[0]
    y_size = norm_data.shape[1]
    num_slices = norm_data.shape[2]    

    # Now calculate the actual AUC
    auc_data = numpy.empty([x_size,y_size,num_slices])
    
    for slice in range(num_slices):
        auc_data[:,:,slice] = scipy.integrate.simps(norm_data[:,:,slice,inj_point:inj_point+auc_reps],x=time_points)
       
    # Deal with bounding boxes

    if bbox is None:        
        bbox = numpy.array([0,x_size-1,0,y_size-1])    
       
    if bbox.shape == (4,):            
    
        bbox_mask = numpy.empty([x_size,y_size])
        bbox_mask[:] = numpy.nan        
        bbox_mask[bbox[0]:bbox[1],bbox[2]:bbox[3]] = 1
    
        # First tile for slice
        bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size,1),num_slices)

    else:      
        raise ValueError('Please supply a bbox for h_calculate_AUC')   
     
    return auc_data*bbox_mask


def h_normalize_dce(scan_object, bbox = None, pdata_num = 0):

    ########### Getting and defining parameters
    
    # Data
    data = scan_object.pdata[pdata_num].data

    x_size = data.shape[0]
    y_size = data.shape[1]
    num_slices = getters.get_num_slices(scan_object,pdata_num)
    
    # Method params
    #TODO: change this so it doesn't require method file WIHOUT BREAKING IT!
    reps =  scan_object.method.PVM_NRepetitions
    
    if reps != scan_object.pdata[pdata_num].data.shape[-1]:
        reps = scan_object.pdata[pdata_num].data.shape[-1]
        print('\n \n ***** Warning **** \n \n !!! Incomplete dce data for {0} \n \n'.format(scan_object.shortdirname) )

    ## Check for bbox traits and create bbox_mask to output only partial data

    if bbox is None:        
        bbox = numpy.array([0,x_size-1,0,y_size-1])

    if bbox.shape == (4,):            
    
        bbox_mask = numpy.empty([x_size,y_size])

        bbox_mask[:] = numpy.nan        
        bbox_mask[bbox[0]:bbox[1],bbox[2]:bbox[3]] = 1
    
        # First tile for slice
        bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size,1),num_slices)
        # Next tile for reps
        bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size,num_slices,1),reps)

    else:      
        raise ValueError('Please supply a bbox for h_normalize_dce')

    # Calculated params      
    inj_point = sarpy.fmoosvi.analysis.h_inj_point(scan_object)
    
    norm_data = numpy.empty([x_size,y_size,num_slices,reps])


    for slice in range(num_slices):
        baseline = numpy.mean(data[:,:,slice,0:inj_point],axis=2)
        norm_data[:,:,slice,:] = (data[:,:,slice,:] / numpy.tile(baseline.reshape(x_size,y_size,1),reps))-1

    return norm_data*bbox_mask
 
def h_enhancement_curve(scan_object, adata_roi_label, pdata_num = 0):

    try:
        norm_data = sarpy.fmoosvi.analysis.h_normalize_dce(scan_object)
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
         
        masked_data = data_scan.data * numpy.tile(numpy.reshape(roi,[
roi.shape[0], roi.shape[1], roi.shape[2],1]),reps)

        enhancement_curve = numpy.empty(shape = [num_slices, 2, reps])
        
        for slice in range(num_slices):
            
            enhancement_curve[slice,0,:] = time
            enhancement_curve[slice,1,:] = scipy.stats.nanmean(scipy.stats.nanmean(masked_data[:,:,slice,:], axis=0), axis =0)            
        return enhancement_curve
    except:    
        print("Perhaps you didn't pass in a valid mask or passed bad data")
        raise


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

def h_calculate_AUGC(scan_object, adata_label, bbox = None, time = 60, pdata_num = 0):
    
    """
    Returns an area under the gadolinium concentration curve adata for the scan object

    :param object scan_object: scan object from a study
    :param str adata_label: indicates the label with which gd_conc is stored
    :param integer pdata_num: reconstruction number, according to python numbering.
            default reconstruction is pdata_num = 0.
    :return: array with augc data
    """

    # Get the concentration data stored as an adata

    try:
        data = scan_object.adata[adata_label].data
    except KeyError:
        print('h_caculate_AUGC: Source data {0} does not exist yet.'.format(adata_label))
        raise KeyError
    
    ########### Getting and defining parameters
    
    # Visu_pars params
    num_slices = getters.get_num_slices(scan_object,pdata_num)
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
    inj_point = sarpy.fmoosvi.analysis.h_inj_point(scan_object)
    
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
       
    if bbox.shape == (4,):            
    
        bbox_mask = numpy.empty([x_size,y_size])
        bbox_mask[:] = numpy.nan        
        bbox_mask[bbox[0]:bbox[1],bbox[2]:bbox[3]] = 1
    
        # First tile for slice
        bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size,1),num_slices)

    else:      
        raise ValueError('Please supply a bbox for h_calculate_AUGC')   

        # If this gives a value error about operands not being broadcast together, go backand change your adata to make sure it is squeezed
    return augc_data*bbox_mask     
    
def h_conc_from_signal(scan_object, scan_object_T1map, 
                       adata_label = 'T1map_LL', bbox = None,
                       relaxivity=4.3e-3, pdata_num = 0):

    ########### Getting and defining parameters
    
    # Data
    data = scan_object.pdata[pdata_num].data
    
    # resample the t1map onto the dce
    data_t1map_pre = scan_object_T1map.adata[adata_label]
    
    data_t1map = sarpy.ImageProcessing.resample_onto.resample_onto_pdata(data_t1map_pre,scan_object.pdata[pdata_num],use_source_dims=True)

    x_size = data.shape[0]
    y_size = data.shape[1]
    num_slices = getters.get_num_slices(scan_object,pdata_num)
    
    # Method params
    #TODO: change this so it doesn't require method file WIHOUT BREAKING IT!
    reps =  scan_object.method.PVM_NRepetitions

    if reps != scan_object.pdata[pdata_num].data.shape[-1]:
        reps = scan_object.pdata[pdata_num].data.shape[-1]
        print('\n \n ***** Warning **** \n \n !!! Incomplete dce data for {0} \n \n'.format(scan_object.shortdirname) )
    
    
    TR = scan_object.method.PVM_RepetitionTime
    FA = scan_object.acqp.ACQ_flip_angle
    
    inj_point = sarpy.fmoosvi.analysis.h_inj_point(scan_object, pdata_num = 0)    
    
    # Deal with bounding boxes
    if bbox is None:        
        bbox = numpy.array([0,x_size-1,0,y_size-1])    
       
    if bbox.shape == (4,):            
    
        bbox_mask = numpy.empty([x_size,y_size])
        bbox_mask[:] = numpy.nan        
        bbox_mask[bbox[0]:bbox[1],bbox[2]:bbox[3]] = 1
    
        # First tile for slice
        bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size,1),num_slices)

    else:      
        raise ValueError('Please supply a bbox for h_conc_from_signal')   

    T1 = numpy.empty([x_size,y_size,num_slices,reps])
    T1[:] = numpy.nan
    
    for x in xrange(bbox[0],bbox[1]):
        for y in xrange(bbox[2],bbox[3]):
            for slice in xrange(num_slices):
                               
                baseline_s = numpy.mean(data[x,y,slice,0:inj_point])
                E1 = numpy.exp(-TR/data_t1map[x,y,slice])
                c = numpy.cos(numpy.radians(FA))
                
                T1[x,y,slice,0:inj_point] = data_t1map[x,y,slice]
                
# Use the SPGR equation twice, once before agent admin. and once after. If you
# divide the two, the M0s cancel out, and so do the sin thetas. what remains is:
# (1-E1)*(1-E0*cos) / ((1-E0) * (1-E1*cos)). use wolfram alpha to solve this:
# http://www.wolframalpha.com/input/?i=solve+%28%281-x%29*%281-b*c%29%29+%2F+%28%281-b%29*%281-x*c%29%29+%3D+r+for+x

                for rep in xrange(inj_point,reps):                    
                    s = data[x,y,slice,rep] / baseline_s
                    E2 = (-E1*c + E1*s - s + 1) / (E1*s*c - E1*c - s*c +1)
                    T1[x,y,slice,inj_point:] = -TR / numpy.log(E2)

    # If this gives a value error about operands not being broadcast together, go backand change your adata to make sure it is squeezed

    T1baseline = numpy.squeeze(data_t1map)*bbox_mask
    T1baseline = numpy.tile(T1baseline.reshape(x_size,y_size,num_slices,1),reps)
    conc = (1/relaxivity) * ( (1/T1) - (1/T1baseline) )
    
    conc[conc<0] = 0
       
    return conc
    
    
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

def h_func_T1(params,t):
    M,B,T1_eff = params
    return numpy.abs(M*(1-B*numpy.exp(-t/T1_eff)))
    
def h_within_bounds(params,bounds):
    try:        
        for p in xrange(len(params)):
            if params[p] >= bounds[p,0] or params[p] <= bounds[p,1] :
                return True
        return False
    except:
        print('You have some funky inputs for the bounds, you fail.')
        return False
            
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

def h_fit_T1_LL(scan_object, bbox = None, flip_angle_map = 0, pdata_num = 0, 
                params = []):
    
    if len(params) == 0:      
        params = [3E5, 2, 350]
  
    if type(flip_angle_map) != numpy.ndarray:
        flip_angle_map = math.radians(scan_object.acqp.ACQ_flip_angle)
        
    data = scan_object.pdata[pdata_num].data[:]

    data_after_fitting = numpy.zeros( [data.shape[0],\
                                       data.shape[1],\
                                       data.shape[2]] )
                                       
    repetition_time = scan_object.method.PVM_RepetitionTime
    inversion_time = scan_object.method.PVM_InversionTime   
    x_size = data.shape[0]
    y_size = data.shape[1]
    num_slices = getters.get_num_slices(scan_object,pdata_num)                                        
                                       
                                       
    fit_results = numpy.array(data_after_fitting[:], dtype=dict)                       
            
    t_data = numpy.linspace(inversion_time,\
        scan_object.pdata[pdata_num].data.shape[3]*repetition_time,\
        scan_object.pdata[pdata_num].data.shape[3])
   
    ## Check for bbox traits and create bbox_mask to output only partial data

    if bbox is None:        
        bbox = numpy.array([0,x_size-1,0,y_size-1])

    if bbox.shape != (4,):    
        raise ValueError('Please supply a bbox for h_fit_T1_LL')      
        
    # Start the fitting process        
        
    data_after_fitting = numpy.empty([x_size,y_size,num_slices])
    data_after_fitting[:] = numpy.nan

#    fit_results = numpy.empty([x_size,y_size,num_slices])
    
    for x in xrange(bbox[0],bbox[1]):
        for y in range(bbox[2],bbox[3]):
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

    return numpy.squeeze(T1), fit_results
                    
        #TODO: Implement code to deal with other methos of calculating T1
        # e.g., IR, VFA
        # NecS3Exp= sarpy.Experiment('NecS3')
        # scan_object = NecS3Exp.studies[0].find_scan_by_protocol('04_ubcLL2')



### Other Helpers


def h_phase_from_fid(scan_object):
    
    phase_data = numpy.angle(scipy.fftpack.fftshift(scipy.fftpack.fftn(scipy.fftpack.fftshift(scan_object.fid))))
    
    return phase_data
    
def h_mag_from_fid(scan_object):

    mag_data = numpy.abs(scipy.fftpack.fftshift(scipy.fftpack.fftn(scipy.fftpack.fftshift(scan_object.fid))))
    
    return mag_data
    
def h_image_to_mask(roi_data, background=None, foreground=None):
    
    if background is None:
        background = numpy.nan
    if foreground is None:
        foreground = 1

    roi_mask = copy.deepcopy(roi_data)
 
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
            


def h_goodness_of_fit(data,infodict, indicator = 'rsquared'):
    
    if indicator == 'rsquared':
        ss_err=(infodict['fvec']**2).sum()
        ss_tot=((data-data.mean())**2).sum()
        rsquared=1-(ss_err/ss_tot)
               
        return rsquared
        
    else:
        print ('There is no code to produce that indicator. Do it first.')
        raise Exception

def h_generate_VTC(scan, bbox = None, pdata_num = 0):

    # Normalize data
    #ndata = sarpy.fmoosvi.analysis.h_normalize_dce(scan)
    
    # Try without normalization
    ndata = scan.pdata[pdata_num].data
    
   # Get useful params        
    x_size = ndata.shape[0]
    y_size = ndata.shape[1]
    num_slices = ndata.shape[2]
    reps = ndata.shape[3]
    
    mask = numpy.empty([x_size,y_size,num_slices,reps])
    mask[:] = numpy.nan
    # Set bounding boxes and get ready to join together
    mask[bbox[0]:bbox[1],bbox[2]:bbox[3],:,:] = 1
        
    ndata[:,:,:,-1] = numpy.nan
    
    ndata = mask * ndata
    # Reshape it  to stitch together all the data
    nrdata = numpy.empty([x_size,y_size*reps,num_slices])
    
    for s in xrange(num_slices):
        nrdata[:,:,s] = ndata[:,:,s,:].reshape([x_size,y_size*reps])
        
    return nrdata

def h_generate_ROI(masterlist_name, data_label, adata_label = None, 
                 ioType = None, path = None, forceVal = False):

    root = os.path.join(os.path.expanduser('~/mdata'),
                        masterlist_name,
                        masterlist_name)
  
    fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str) 
    if path is None:
        path = os.path.expanduser(os.path.join('~','mdata',masterlist_name))
    
    for k,v in master_list.iteritems():
              
        try:
            scan = sarpy.Scan(v[data_label][0])
            sdir = scan.shortdirname
            sdir = sdir.replace('/','_')
            
            if ioType == 'export':
                fname = os.path.join(path, sdir + '.nii')
                scan.pdata[0].export2nii(fname)
            
            elif ioType == 'import':
                    
                if adata_label is None:
                    adata_label = 'roi'
                    
                if (not adata_label in scan.adata.keys()) or forceVal is True:
                    
                    fname = os.path.join(path, sdir + 'a.nii')
                    roi = nibabel.load(fname).get_data()[:,:,:,0]
                    
                    # the default foreground and background in h_image_to_mask
                    # will result in a roi_m that has NaN and 1 only (aka
                    # 'proper mask')
                    roi_m = sarpy.fmoosvi.analysis.h_image_to_mask(roi)               
                    
                    scan.store_adata(key=adata_label, data = roi_m, force = forceVal)
                    print('h_generate_roi: saved {0} roi label'.format(scan.shortdirname))
                else:
    
                    print('{0}: adata already exists {1} '.format(
                    adata_label,scan.shortdirname))
                    pass 
    
            else:
                
                print("Please specify either 'import' or 'export' \
                for the ioType!")
        
        except IOError:
            
            print('\n \n ** WARNING ** \n \n Not found: {0} and {1} \n'.format(k,data_label) )
            pass

        except:
            raise
            
    print('Nifti images were processed in {0}'.format(path))


## Not working, or of unknown reliability

    
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
    
    
