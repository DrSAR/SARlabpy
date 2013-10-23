# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:11:03 2013

@author: fmoosvi
"""

import sarpy
import sarpy.fmoosvi.analysis
import sarpy.fmoosvi.getters
import pylab
import numpy
import os
import json
import re
import sarpy.ImageProcessing.resample_onto
import collections
import scipy
import scipy.stats

 


def store_deltaT1(masterlist_name,
                  T1map_1,
                  T1map_2,
                  adata_label1 = 'T1map_LL',
                  adata_label2 = None,
                  analysis_label='deltaT1',
                  forceVal = False):

    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)

    fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str)

    if adata_label2 is None:
        adata_label2 = adata_label1

    for k,v in master_list.iteritems():

        try:
            scan1 = sarpy.Scan(v[T1map_1][0])
            scan2 = sarpy.Scan(v[T1map_2][0])
            bbox = sarpy.fmoosvi.getters.get_bbox(v, T1map_2)

    	except IOError:
    	    print('{0}: Not found adata {1} or {2} in patient {3}'.format(
    	          analysis_label,adata_label1,adata_label2,k) )
    	    continue

        if (not analysis_label in scan2.adata.keys()) or forceVal is True:
                
            deltaT1 = scan1.adata[adata_label1].data - scan2.adata[adata_label2].data
            scan2.store_adata(key=analysis_label, data = deltaT1, force = forceVal)

            print('{0}: Success. Saved {1}'.format(analysis_label,
                                  scan2.shortdirname))
                                  
        else:
            print('{0}: adata already exists {1} '.format(
            analysis_label,scan2.shortdirname))
            pass

### ROI based code

def roi_average(masterlist_name,
                data_label,
                adata_label,
                roi_label,
                analysis_label=None,
                forceVal = False):
                   
    if analysis_label is None:
        analysis_label = adata_label + '_avg'                   

    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)

    fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str)

    for k,v in master_list.iteritems():
        
        try:
            scan = sarpy.Scan(v[data_label][0])           
            data = scan.adata[adata_label].data
            roi = sarpy.Scan(v[data_label][0]).adata[roi_label].data
        except(KeyError, IOError):
            print('{0}: Not found adata {1} or {2} in patient {3}'.format(
                                    analysis_label,adata_label,roi_label,k) )
            continue
       
        roi_data = data * roi
        
        # I think this might be necessary to ensure that values don't get carried
        # over from previous iterations?
        try:
            del avg_T1,weights
        except NameError:
            pass
        avg_T1 = []
        weights = []
       
        for slice in xrange(roi_data.shape[-1]):
            avg_T1.append(scipy.stats.nanmean(roi_data[:,:,slice].flatten()))
            weights.append(numpy.nansum(roi[:,:,slice]))
            
        if (not analysis_label in scan.adata.keys()) or forceVal is True:
                
            scan.store_adata(key=analysis_label, 
                             data = numpy.array(avg_T1), 
                             force = forceVal)

            scan.store_adata(key=analysis_label+'_weights', 
                             data = numpy.array(weights), 
                             force = forceVal)                             

            print('{0}: Success. Saved {1}'.format(analysis_label,
                                  scan.shortdirname))
                                  
        else:
            print('{0}: adata already exists {1} '.format(
            analysis_label,scan.shortdirname))
            pass          

def bulk_transfer_roi(masterlist_name, dest_adata_label, roi_src_adata_label = None, tag = None, forceVal = False):
    '''
    Move an ROI from one scan to another. E.g., Moving an roi from an anatomy 
    scan to a LL scan. 
    
    dest_adata_label: Label specifying the destination of the roi transfer
    
    tag: string to prepend before ROI. useful if auc60_roi exists and you want 
            auc60_old_roi
    
    Note: this routine looks for the special 'roi' label. This is the canonical
    label used to specify ROIs and other scans will 'inherit' - by way of resample
    - this roi. 

    Example:
        
    sarpy.fmoosvi.wrappers.bulk_transfer_roi('NecS3',dest_adata_label,
                                             forceVal = False)
    '''
    
    # Get the name of the source roi
    if roi_src_adata_label is None:
        roi_src_adata_label = 'roi'

    # Set the name of the destination scan

    if tag is None:
        dest_label = dest_adata_label + '_roi'
    else:
        dest_label = dest_adata_label + str(tag) + '_roi'

#    # Get all the studies in the experiment

#    exp = sarpy.Experiment(exp_name)
    
    
    
    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)

    fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str)

    for k in master_list.iteritems():
        
        for study in sarpy.Patient(k[0]).studies:
            adata_dict = study.find_adata_scans()
        
            try:
                src = sarpy.Scan(adata_dict[roi_src_adata_label][0])
    
            except(IOError):
                print('bulk_transfer_roi: IO Problem with roi src scan for study {0}'.format(study.shortdirname))
                raise
            except KeyError:
                print('bulk_transfer_roi: Key roblem with roi src scan for study {0}, maybe this study is extra'.format(study.shortdirname))


    
            if dest_adata_label not in adata_dict:
                # Skip this Patient
                print('bulk_transfer_roi: {0} not found in study {1}'.format(dest_adata_label,study.shortdirname))                
                continue
                
                
            for d in adata_dict[dest_adata_label]:
                dest = sarpy.Scan(d)
           
                if (not dest_label in dest.adata.keys()) or forceVal is True:
                    
                    dest_pd = dest.pdata[0]
                    roi = src.adata[roi_src_adata_label]
                    
                    resampled_roi = sarpy.ImageProcessing.resample_onto.\
                        resample_onto_pdata(roi, dest_pd, replace_nan=0)
                        
                    places = numpy.where(resampled_roi < .5)
                    other_places = numpy.where(resampled_roi >= .5)
                    resampled_roi[places] = numpy.nan
                    resampled_roi[other_places] = 1
    
                    dest.store_adata(key = dest_label, data = resampled_roi, 
                                     force = forceVal)   
                                     
                    print('bulk_transfer_roi: Success. Saved {0} in {1}'.format(dest_label, 
                                  dest.shortdirname))                                 
                else:
                    print('bulk_transfer_roi: {0} adata already exists {1}'.format(dest_label, 
                          dest.shortdirname))
                    pass      

    #TODO: Fix this function so it uses the masterlist    
    print('Firas: you really need to fix this function so it uses the masterlist')

def plotVTC(masterlist_name, key, data_label, adata_label = None):

    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)
    if os.path.exists(os.path.join(root+'_updated.json')):
        fname_to_open = root+'_updated.json'
    else: 
        fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str)
    value = master_list[key]
        
    data = sarpy.Scan(value[data_label][0])
    bbox_px = sarpy.fmoosvi.getters.get_bbox(value,data_label)
    
    if adata_label is None:
        nrdata = sarpy.Scan(value[data_label][0]).adata['vtc'].data        
    else:      
        nrdata = sarpy.Scan(value[data_label][0]).adata[adata_label].data
        
    num_slices = nrdata.shape[2]
 
    x_bbox_px = bbox_px[1] - bbox_px[0] +1

    # Handle the special case of a corner bbox at origin
    if bbox_px[0]==0: x_bbox_px - 1
    
#    pylab.figure(figsize = [10,data.method.PVM_FovCm[1]/data.method.PVM_FovCm[0]*10])
#    pylab.imshow(data.pdata[0].data[:,:,3,80])
#   
    for slice in xrange(num_slices):
        G = pylab.matplotlib.gridspec.GridSpec(x_bbox_px, 1)
        fig = pylab.figure(figsize = [10,data.method.PVM_FovCm[1]/data.method.PVM_FovCm[0]*10])       
    
        i = 0
        
        for x in xrange(bbox_px[0],bbox_px[1]):
         
            fig.add_subplot(G[i],frameon=False,xticks=[],yticks=[])
        
            d = nrdata[x,:,slice]
        
            pylab.plot(d,'-',linewidth=0.1)
            pylab.ylim([-0.2,2])
            
            i+=1 # increment counter for the subplot (x)             