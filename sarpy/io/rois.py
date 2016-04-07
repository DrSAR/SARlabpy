#!/usr/bin/env python 
# -*- coding: utf-8 -*-

"""
Created on Thu May 30 17:49:35 2013

@author: fmoosvi
"""

import sarpy.fmoosvi.getters

def export_roi(masterlist_name,
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
    import sarpy.fmoosvi.analysis
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

    for k in exp.patients.keys():
        # Search for all the labels that have roi and create a scan label list of rois

        roi_labels = [r for r in exp.labels.keys() if 'roi' in r]

        # Iterate over the roi list, import or export the rois
        for r_lbl in roi_labels:
                  
            try:
                scn_name = exp.patients[k][r_lbl]
                scan = sarpy.Scan(scn_name)
                
            except IOError:    
                print('\n \n ** WARNING ** \n \n Not found: {0} and {1} \n'.format(k,r_lbl) )
                continue

            sdir = scn_name.replace('/','_')  

            fname = os.path.join(path, sdir + '.nii')

            # Export the image as a niftii file            

            # Note: std_modifier is a value that controls scaling of the image 
            # Look at code in getters.get_image_clims(data,std_modifier=None):
            # Basically controls the upper limit of the default image scaling

            scan.pdata[0].export2nii(fname,rescale,std_modifier)

def import_roi(masterlist_name,
               roi_suffix = None,
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
    import sarpy.fmoosvi.analysis
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

    for k in exp.patients.keys():
        # Search for all the labels that have roi and create a scan label list of rois

        roi_labels = [r for r in exp.labels.keys() if 'roi' in r]

        # Iterate over the roi list, import or export the rois
        for r_lbl in roi_labels:
                  
            try:
                scn_name = exp.patients[k][r_lbl]
                scan = sarpy.Scan(scn_name)
                
            except IOError:    
                print('\n \n ** WARNING ** \n \n Not found: {0} and {1} \n'.format(k,r_lbl) )
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

            # the default foreground and background in h_image_to_mask
            # will result in a roi_m that has NaN and 1 only (aka
            # 'proper mask')
            roi_m = sarpy.fmoosvi.analysis.h_image_to_mask(roi,peaks=peaks)

            # Calculate the weights for each slice
            weights = numpy.zeros(shape=roi_m.shape[-1])
            ROIpx = numpy.nansum(roi_m.flatten())

            for sl in xrange(roi_m.shape[-1]):
                weights[sl] = numpy.divide(numpy.nansum(roi_m[:,:,sl].flatten()),ROIpx)

            try:
                scan.store_adata(key=roi_key_name, data = roi_m, force = forceVal)
                scan.store_adata(key=roi_key_name+'_weights', data = weights, force = forceVal)


                # get the bbox for this roi and scan
                bbox = sarpy.fmoosvi.getters.get_roi_bbox(scn_name,
                                                      roi_adata_label=roi_key_name)
                scan.store_adata(key='bbox',data=bbox,force=forceVal)

                print('h_generate_roi: saved {0} roi label as {1}'.format(scan.shortdirname,roi_key_name))
            except AttributeError:
                print('h_generate_roi: force save of scan {0} is required'.format(scan.shortdirname))

    # Save adata in all scans 
    for k in exp.patients.keys():

        pat = exp.patients[k]

        try:
            bbox = sarpy.Scan(pat['roi']).adata['bbox'].data

            for lbl,scn in pat.iteritems():

                scn = sarpy.Scan(scn)
                scn.store_adata(key='bbox', data = bbox,force=True)
        except:
            raise

    print('Nifti images were processed in {0}'.format(path))