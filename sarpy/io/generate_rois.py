#!/usr/bin/env python 
# -*- coding: utf-8 -*-

"""
Created on Thu May 30 17:49:35 2013

@author: fmoosvi
"""

#TODO: for some reason the required args still show up as optional


def generate_rois(masterlist_name, ioType, path = None, 
                  rescale = None, std_modifier = None, 
                  peaks = None, forceVal = False):

    import json
    import nibabel
    import collections
    import sarpy.fmoosvi.analysis
    import os

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
        path = os.path.expanduser(os.path.join('~','sdata',masterlist_name,'rois'))
    
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
                scan.pdata[0].export2nii(fname,rescale,std_modifier)
            
            elif ioType == 'import':
                                        
                fname = os.path.join(path, sdir + 'a.nii')
                roi = nibabel.load(fname).get_data()[:,:,:,0]
                
                # the default foreground and background in h_image_to_mask
                # will result in a roi_m that has NaN and 1 only (aka
                # 'proper mask')
                roi_m = sarpy.fmoosvi.analysis.h_image_to_mask(roi,peaks=peaks)  

                try:             
                    scan.store_adata(key='roi', data = roi_m, force = forceVal)
                    print('h_generate_roi: saved {0} roi label'.format(scan.shortdirname))
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