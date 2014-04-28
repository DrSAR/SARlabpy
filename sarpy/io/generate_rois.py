#!/usr/bin/env python 
# -*- coding: utf-8 -*-

"""
Created on Thu May 30 17:49:35 2013

@author: fmoosvi
"""

#TODO: for some reason the required args still show up as optional


def generate_rois(masterlist_name, 
                  ioType, 
                  roi_suffix = None, 
                  path = None,
                  rescale = None,
                  std_modifier = None, 
                  peaks = None, 
                  forceVal = False):
    # e.g., roi_suffix = 'b' or 'c' or 'd' or 'e' etc...

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
        roi_suffix = ''

    root = os.path.join(os.path.expanduser('~/sdata'),
                    masterlist_name,
                    masterlist_name)
  
    fname_to_open = root+'.json'
    
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str) 
    if path is None:
        if roi_suffix == '':
            roi_key_name = 'roi'
        else:
            roi_key_name = 'roi' + '_' + roi_suffix        
        path = os.path.expanduser(os.path.join('~','sdata',masterlist_name,roi_key_name))
    
    for k,v in master_list.iteritems():

        # Search for all the labels that have roi and create an roi list 
        roi_labels = [r for r in master_list[k].keys() if 'roi' in r]

        # Iterate over the roi list, import or export the rois
        for r_lbl in roi_labels:
                  
            try:
                scan = sarpy.Scan(v[r_lbl][0])
                
            except IOError:    
                print('\n \n ** WARNING ** \n \n Not found: {0} and {1} \n'.format(k,r_lbl) )
                continue
        
            sdir = scan.shortdirname
            sdir = sdir.replace('/','_')        
            
            if ioType == 'export':
                fname = os.path.join(path, sdir + '.nii')

                # Note: std_modifier is a value that controls scaling of the image 
                # Look at code in getters.get_image_clims(data,std_modifier=None):
                # Basically controls the upper limit of the default image scaling

                scan.pdata[0].export2nii(fname,rescale,std_modifier)
            
            elif ioType == 'import':
                                                        
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
                    print('h_generate_roi: saved {0} roi label as {1}'.format(scan.shortdirname,roi_key_name))
                except AttributeError:
                    print('h_generate_roi: force save of scan {0} is required'.format(scan.shortdirname))

            else:          
                print("Please specify either 'import' or 'export' \
                for the ioType!")
                
    print('Nifti images were processed in {0}'.format(path))


if __name__ == '__main__':

    import argparse
    import os
    import numpy

    # we have been called as a script
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-m', '--masterlist_name', type=str, required=True,
                       help='Name of the masterlist file. e.g. NecS3 \
                       Usage: python generate_rois.py -m HerP2 -d axref -a roi -i export ')

    parser.add_argument('-s', '--suffix', type=str, required=True,
                       help='Optional: adds _b, _c, etc... to rois. only insert the letter after _')                                       
    
    parser.add_argument('-i', '--iotype', type=str, required=True,
                       choices = ['import','export'], help='IOType, import or export')    
                                         
    parser.add_argument('-p', '--path', type=str,
                       help='Optional: Path to/from export/import data, default to sdata')

    parser.add_argument('-k', '--peaks', type=str,
                       help='Optional: Number of peaks in the roi')
                       
    parser.add_argument('-f', '--force', type=str, choices = ['True','False'],
                        help='Optional: Replace data? True or False, defaults to False')   
                       
    args = parser.parse_args()
    masterlist_name = args.masterlist_name
    ioType = args.iotype
    
    #TODO: Integrate this intothe parser for cleaner code!

    try:
        peaks = numpy.int(args.peaks)
    except AttributeError:
        peaks = None
    
    try:
        path = os.path.expanduser(os.path.join(args.path))
    except AttributeError:
        path = os.path.expanduser(os.path.join('~','sdata',masterlist_name))
    
    try:
        force = args.force
        force = True
    except AttributeError:
        force = False

    generate_rois(masterlist_name, ioType, path, peaks = peaks, forceVal = force)