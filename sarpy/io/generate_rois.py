#!/usr/bin/env python 
# -*- coding: utf-8 -*-

"""
Created on Thu May 30 17:49:35 2013

@author: fmoosvi
"""

#TODO: for some reason the required args still show up as optional


def generate_rois(masterlist_name, data_label, 
                  adata_label, ioType, path, forceVal = False):

    import json
    import nibabel
    import collections
    import sarpy.fmoosvi.analysis
    import argparse
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
              
        try:
            scan = sarpy.Scan(v[data_label][0])
            
        except IOError:    
            print('\n \n ** WARNING ** \n \n Not found: {0} and {1} \n'.format(k,data_label) )
            continue
    
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
            
    print('Nifti images were processed in {0}'.format(path))


if __name__ == '__main__':


    import argparse
    import os

    # we have been called as a script
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-m', '--masterlist_name', type=str, required=True,
                       help='Name of the masterlist file. e.g. NecS3 \
                       Usage: python generate_rois.py -m HerP2 -d axref -a roi -i export ')
    
    parser.add_argument('-d', '--data_label', type=str, required=True,
                       help='Data label, usually anatomy or IR_A')                   
    
    parser.add_argument('-i', '--iotype', type=str, required=True,
                       choices = ['import','export'], help='IOType, import or export')    
    
    parser.add_argument('-a', '--adata_label', type=str,
                       help='Optional: adata label, defaults to roi') 
                                         
    parser.add_argument('-p', '--path', type=str,
                       help='Optional: Path to/from export/import data, default to sdata')
                       
    parser.add_argument('-f', '--force', type=str, choices = ['True','False'],
                        help='Optional: Replace data? True or False, defaults to False')   
                       
    args = parser.parse_args()
    masterlist_name = args.masterlist_name
    data_label = args.data_label
    ioType = args.iotype
    
    #TODO: Integrate this intothe parser for cleaner code!
    
    try:
        adata_label = args.adata_label
    except AttributeError:
        adata_label = 'roi'
    
    try:
        path = os.path.expanduser(os.path.join(args.path))
    except AttributeError:
        path = os.path.expanduser(os.path.join('~','sdata',masterlist_name))
    
    try:
        force = args.force
        force = True
    except AttributeError:
        force = False

    generate_rois(masterlist_name, data_label, 
                  adata_label, ioType, path, forceVal = force)