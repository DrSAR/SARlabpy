#!/usr/bin/env python 
# -*- coding: utf-8 -*-

"""
Created on Thu May 30 17:49:35 2013

@author: fmoosvi
"""

import sarpy.analysis.getters

def export_roi(masterlist_name,
               roi_scan_label = None,
               roi_suffix = None,
               path = None,
               rescale = None,
               std_modifier = None,
               forceVal = False):
    '''
    e.g., roi_suffix = 'b' or 'c' or 'd' or 'e' etc...

    '''
    import json
    import nibabel
    import collections
    import sarpy.analysis.analysis
    import os
    import scipy
    import numpy

    # Quick check to see if the name of the roi matches convention
    if roi_suffix is not None:
        assert(type(roi_suffix) is str and len(roi_suffix)==1)
    else:
        roi_suffix = '' # AA Don't want to set this to 'a' for backwards compatibility

    # Set the path and assume a default one if it doesn't exist
    if path is None:
        if roi_suffix == '':
            roi_key_name = 'roi'
        else:
            roi_key_name = 'roi' + '_' + roi_suffix        
        path = os.path.expanduser(os.path.join('~','sdata',masterlist_name,roi_key_name))
    
    # Create the experiment
    exp = sarpy.Experiment.from_masterlist(masterlist_name+'.config')

    for k in list(exp.patients.keys()):

        if roi_scan_label is None:
            # Search for all the labels that have roi and create a scan label list of rois
            roi_labels = [r for r in list(exp.labels.keys()) if 'roi' in r]
        else:
            # This overrides in case there wasn't an explict roi_scan acquired.
            roi_labels = [roi_scan_label]

        # Iterate over the roi list, import or export the rois
        for r_lbl in roi_labels:
                  
            try:
                scn_name = exp.patients[k][r_lbl]
                scan = sarpy.Scan(scn_name)

                sdir = scn_name.replace('/','_')  

                fname = os.path.join(path, sdir + '.nii')

                # Export the image as a niftii file            

                # Note: std_modifier is a value that controls scaling of the image 
                # Look at code in getters.get_image_clims(data,std_modifier=None):
                # Basically controls the upper limit of the default image scaling

                scan.pdata[0].export2nii(fname,rescale,std_modifier)                
                
            except(IOError,KeyError):    
                print(('\n \n ** WARNING ** \n \n Not found (check permissions): {0} and {1} \n'.format(k,r_lbl) ))
                pass

def import_roi(masterlist_name,
               roi_suffix = None,
               roi_scan_label = None,
               path = None,
               rescale = None,
               std_modifier = None,
               peaks = None,
               forceVal = False):
    '''
    e.g., roi_suffix = 'b' or 'c' or 'd' or 'e' etc...

    '''
    import json
    import nibabel
    import collections
    import sarpy.analysis.analysis
    import os
    import scipy
    import numpy

    # Quick check to see if the name of the roi matches convention
    if roi_suffix is not None:
        assert(type(roi_suffix) is str and len(roi_suffix)==1)
    else:
        roi_suffix = '' # AA Don't want to set this to 'a' for backwards compatibility

    # Set the path and assume a default one if it doesn't exist
    if path is None:
        if roi_suffix == '':
            roi_key_name = 'roi'
        else:
            roi_key_name = 'roi' + '_' + roi_suffix        
        path = os.path.expanduser(os.path.join('~','sdata',masterlist_name,roi_key_name))
    
    # Create the experiment
    exp = sarpy.Experiment.from_masterlist(masterlist_name+'.config')

    for k in list(exp.patients.keys()):
        if roi_scan_label is None:
            # Search for all the labels that have roi and create a scan label list of rois
            roi_labels = [r for r in list(exp.labels.keys()) if 'roi' in r]
        else:
            # This overrides in case there wasn't an explict roi_scan acquired.
            roi_labels = [roi_scan_label]

        # Iterate over the roi list, import or export the rois
        for r_lbl in roi_labels:
                  
            try:
                scn_name = exp.patients[k][r_lbl]
                scan = sarpy.Scan(scn_name)                
            except KeyError:    
                logger.error('\n \n ** WARNING ** \n \n Not found: {0} and {1} \n'.format(k,r_lbl),
                             exc_info=True)
                continue
            except IOError:    
                logger.error('\n \n ** WARNING ** \n \n Not found: {0} and {1} \n'.format(k,r_lbl),
                             exc_info=True)
                continue

            sdir = scn_name.replace('/','_')  

            # AA cont'd This is to ensure that the first roi in the folder 'roi' 
            # gets saved with adata label 'roi' --> for backwards compatibility

            if roi_suffix == '':
                fname = os.path.join(path, sdir + 'a.nii')
            else:
                fname = os.path.join(path, sdir + roi_suffix + '.nii')

            try:
                roi = nibabel.load(fname).get_data()[:,:,:,0]
            except IndexError:
                roi = nibabel.load(fname).get_data()[:,:,:]

            # Hack for CEST EPI scans because you cannot really draw ROI on any other image
            # So basically I had to export the CEST EPI scan with 1 slice and ~80 freqs as the ROI source
            # This line takes just the first slice where the tumour is drawn. If multislice CEST ever becomes a 
            # TODO: reality, that will need to be taken care of.

            if sarpy.analysis.getters.get_num_slices(scn_name) ==1 and len(roi.shape) > 2:
                roi = roi[:,:,0]

            # the default foreground and background in h_image_to_mask
            # will result in a roi_m that has NaN and 1 only (aka
            # 'proper mask')
            roi_m = sarpy.analysis.analysis.h_image_to_mask(roi,peaks=peaks)

            # Calculate the weights for each slice
            weights = numpy.zeros(shape=roi_m.shape[-1])
            ROIpx = numpy.nansum(roi_m.flatten())

            if len(roi_m.shape) < 3:
                weights = [1.0]
            else:
                for sl in range(roi_m.shape[-1]):
                    weights[sl] = numpy.divide(numpy.nansum(roi_m[:,:,sl].flatten()),ROIpx)

            try:
                scan.store_adata(key=roi_key_name, data = roi_m, force = forceVal)
                scan.store_adata(key=roi_key_name+'_weights', data = weights, force = forceVal)


                # get the bbox for this roi and scan
                bbox = sarpy.analysis.getters.get_roi_bbox(scn_name,
                                                      roi_adata_label=roi_key_name)
                scan.store_adata(key='bbox',data=bbox,force=forceVal)

                print(('h_generate_roi: saved {0} roi label as {1}'.format(scan.shortdirname,roi_key_name)))
            except AttributeError:
                print(('h_generate_roi: force save of scan {0} is required'.format(scan.shortdirname)))

    # Save adata in all scans 
    for k in list(exp.patients.keys()):

        pat = exp.patients[k]

        # The code assumes that there is a scanlabel called roi, fixed so that it can be overridden
        # Still need to fix this !

        if roi_scan_label is None:
            roi_scan_label = 'roi'
        try:
            #bbox = sarpy.Scan(pat[roi_scan_label]).adata['bbox'].data
            bbox = sarpy.analysis.getters.get_roi_bbox(scan.shortdirname,roi_scan_label)

            for lbl,scn in list(pat.items()):

                scn = sarpy.Scan(scn)
                scn.store_adata(key='bbox', data = bbox,force=True)
        except:
            raise
            #print(('No ROI exists for: {0}'.format(k)))

    print(('Nifti images were processed in {0}'.format(path)))

def transfer_roi(masterlist_name,
                 roi_lbl, 
                 scn_lbl, 
                 roi_adata = None, 
                 forceVal = False,
                 adata_save_lbl = None):

    if adata_save_lbl is None:
        adata_save_lbl = 'transferred_roi'

    if roi_adata is None:
        roi_adata = 'roi'

    exp = sarpy.Experiment.from_masterlist(masterlist_name+'.config')

    for scn,roi_scn in zip(exp.labels[scn_lbl],exp.labels[roi_lbl]) :

        try:
        
            roi= sarpy.Scan(roi_scn).adata[roi_adata].data
            scn = sarpy.Scan(scn)
            
            scn.store_adata(key=adata_save_lbl,data = roi,forceVal=forceVal) 

            print(('Successfully transferred roi from {0} to {1}'.format(roi.shortdirname, scn.shortdirname)))
        except:
            print(('FAILED to transfer roi from {0} to {1}'.format(scn, roi_scn)))
            pass