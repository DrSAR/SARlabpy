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

def store_deltaT1(masterlist_name=None,
                  T1map_1=None,
                  T1map_2=None,
                  adata_label1 = 'T1_LL',
                  adata_label2 = None,
                  adata_label3 = 'T1_LL_roi',
                  adata_save_label='deltaT1',
                  force_overwrite = False):

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

    	except IOError:
    	    print('{0}: Not found adata {1} or {2} in patient {3}'.format(
    	          adata_save_label,adata_label1,adata_label2,k) )
    	    continue

        if (not adata_save_label in scan2.adata.keys()) or force_overwrite is True:

            if 'R1' in adata_save_label and 'multiday_' not in adata_save_label:

                deltaR1 = (1/scan2.adata[adata_label2].data - 1/scan1.adata[adata_label1].data)
                scan2.store_adata(key=adata_save_label, data = deltaR1, force = force_overwrite)
                print('{0}: Saved {1}'.format(adata_save_label,
                      scan2.shortdirname))

            # This is the part of the code to subtract T1 voxel by voxel with a baseline T1 from the day before            

            elif 'multiday_' in adata_save_label:

                # Only take the mean of the tumour ROI, not the whole tumour
                baseline = scan1.adata[adata_label1].data * scan1.adata[adata_label3].data
                baseline = scipy.stats.nanmean(baseline,axis=None)

                # Apply the roi to the treated tumour as well
                treated = scan2.adata[adata_label2].data * scan2.adata[adata_label3].data

                if 'T1' in adata_save_label:

                    dat = treated - baseline
                else:
                    dat = 1/treated - 1/baseline

                scan2.store_adata(key=adata_save_label, data = dat, force = force_overwrite)

                print('{0}: Saved {1}'.format(adata_save_label,
                                      scan2.shortdirname))

            else:
                
                deltaT1 = scan2.adata[adata_label2].data - scan1.adata[adata_label1].data
                scan2.store_adata(key=adata_save_label, data = deltaT1, force = force_overwrite)

                print('{0}: Saved {1}'.format(adata_save_label,
                                      scan2.shortdirname))
                                  
        else:
            print('{0}: adata already exists {1} '.format(
            adata_save_label,scan2.shortdirname))
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
            
            return roi_data,avg_T1
            # scan.store_adata(key=analysis_label, 
            #                  data = numpy.array(avg_T1), 
            #                  force = forceVal)

            # scan.store_adata(key=analysis_label+'_weights', 
            #                  data = numpy.array(weights), 
            #                  force = forceVal)                             

            print('{0}: Success. Saved {1}'.format(analysis_label,
                                  scan.shortdirname))
                                  
        else:
            print('{0}: adata already exists {1} '.format(
            analysis_label,scan.shortdirname))
            pass

def bulk_transfer_roi(masterlist_name,
                      dest_adata_label, 
                      roi_src_adata_label = 'roi', 
                      forceVal = False):

    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)

    # Step 1: Load up the masterlist    

    fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str)

    # Step 2: Find all the roi-XXX labels in the masterlist
    # Assumes each patient has all the labels

    roi_scan_labels = []

    for scl in master_list[master_list.keys()[0]].keys():
        if 'roi' in scl:
            roi_scan_labels.append(scl) #.split('-',1)[-1]

    # Step 3: Find all the scans that need to have an roi transferred to them

    slist = sarpy.Experiment(masterlist_name).find_adata_scans(flatten=True)[dest_adata_label]

    if not slist:
        print('No scans found for adata {0}'.format(dest_adata_label))
        return None

    # Step 4: iterate over the search strings
    for search_string in roi_scan_labels:

        # Step 5: Iterate over the slist (scan list where adata was found)
        for scn in slist:

            # Step 6: Check which data label that scan has

            try:
                patname = sarpy.Scan(scn).pdata[0].visu_pars.VisuSubjectId
                src_roi = sarpy.Scan(master_list[patname][search_string][0]).adata[roi_src_adata_label]

            except(AttributeError, IOError):
                print 'adata or scan data not found for {0}'.format(scn,roi_src_adata_label)
                continue
            except KeyError:
                print 'adata or scan data not found for {0}'.format(scn,roi_src_adata_label)
                continue
                
            for v in master_list[patname].keys():
                if master_list[patname][v][0] == scn and search_string.split('-',1)[-1] in v:

                    dest_pd = sarpy.Scan(scn).pdata[0]

                    resampled_roi = sarpy.ImageProcessing.resample_onto.\
                    resample_onto_pdata(src_roi, dest_pd, replace_nan=0)

                    ## WHY DO I HAVE THIS? WHAT IS THIS FOR ?
                    #TODO: FIRAS, PLEASE FIGURE THIS OUT ASAP
                    
                    places = numpy.where(resampled_roi < .5)
                    other_places = numpy.where(resampled_roi >= .5)
                    resampled_roi[places] = numpy.nan
                    resampled_roi[other_places] = 1

                    adata_save_label = dest_adata_label+'_roi'

                    sarpy.Scan(scn).store_adata(key = adata_save_label, data = resampled_roi, 
                                     force = forceVal)   
                                     
                    print('\t Saved {0} in {1}'.format(adata_save_label,scn))   



# def bulk_transfer_roi(masterlist_name,
#                       dest_adata_label, 
#                       roi_src_adata_label = None, 
#                       tag = None, 
#                       forceVal = False):
#     '''
#     Move an ROI from one scan to another. E.g., Moving an roi from an anatomy 
#     scan to a LL scan. 
    
#     dest_adata_label: Label specifying the destination of the roi transfer
    
#     tag: string to prepend before _roi. useful if auc60_roi exists and you want 
#             auc60_old_roi
    
#     Note: this routine looks for the special 'roi' label. This is the canonical
#     label used to specify ROIs and other scans will 'inherit' - by way of resample
#     - this roi. 

#     Example:
        
#     sarpy.fmoosvi.wrappers.bulk_transfer_roi('NecS3',dest_adata_label,
#                                              forceVal = False)
#     '''
    
#     # Set the name of the destination scan

#     if tag is None:
#         dest_label = dest_adata_label + '_roi'
#     else:
#         dest_label = dest_adata_label + str(tag) + '_roi'
   
#     root = os.path.join(os.path.expanduser('~/sdata'),
#                         masterlist_name,
#                         masterlist_name)

#     fname_to_open = root+'.json'
#     with open(os.path.join(fname_to_open),'r') as master_file:
#         json_str = master_file.read()
#         master_list = json.JSONDecoder(
#                            object_pairs_hook=collections.OrderedDict
#                            ).decode(json_str)

#     # Open up the tumourboard and do it tumourboard style
#     for patname,v in master_list.iteritems():

#         # First get all the labels and put it in a list
#         lbl_list = master_list[patname].keys()

#         # Search for all the labels that have roi and create an roi list 
#         roi_labels = [r for r in lbl_list if 'roi' in r]

#         # Iterate over the roi list, get the updated bbox, check for same day-ness
#         for r_lbl in roi_labels:

#             if master_list[patname][r_lbl][0]: # check if this isn't blank

#                 search_string = r_lbl.split('-',1)[-1]

#                 # Iterate over the label list, transfer the roi into the other adata
#                 # Unfortunately this might be a bit slow as you iterate through all 
#                 # the scans
#                 #TODO: Figure out if this is a bottleneck and fix it!
#                 for lbl in lbl_list:

#                     # Check if the 0h/24h is in the label, and if a scan exists
#                     # In other words ignore labels with no scans present

#                     if search_string in lbl and master_list[patname][lbl][0] != '':
#                         scn = sarpy.Scan(master_list[patname][lbl][0])
#                         ad = sarpy.Scan(master_list[patname][lbl][0]).adata.keys()

#                         if dest_adata_label in ad:
#                             dest_pd = scn.pdata[0]

#                             resampled_roi = sarpy.ImageProcessing.resample_onto.\
#                             resample_onto_pdata(src_roi, dest_pd, replace_nan=0)

#                             ## WHY DO I HAVE THIS? WHAT IS THIS FOR ?
#                             #TODO: FIRAS, PLEASE FIGURE THIS OUT ASAP
                            
#                             places = numpy.where(resampled_roi < .5)
#                             other_places = numpy.where(resampled_roi >= .5)
#                             resampled_roi[places] = numpy.nan
#                             resampled_roi[other_places] = 1
        
#                             scn.store_adata(key = dest_label, data = resampled_roi, 
#                                              force = forceVal)   
                                             
#                             print('Saved {0} in {1}'.format(dest_label, 
#                                           scn.shortdirname))    

def plotVTC(masterlist_name, key, data_label, adata_label = None):

    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)
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