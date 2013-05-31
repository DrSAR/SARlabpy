# -*- coding: utf-8 -*-
"""
Created on Thu May 30 17:49:35 2013

@author: fmoosvi
"""

import sarpy
import nibabel
import json
import os

def generate_ROI(masterlist_name, data_label, adata_label = None, 
                 ioType = None, path = None, forceVal = False):

    mdata = os.path.expanduser(os.path.join(
    '~','mdata',masterlist_name+'.json'))
  
    with open(mdata,'r') as master_file:
        master_list = json.load(master_file)

    if path is None:
        path = os.path.expanduser(os.path.join('~','adata'))
    
    for k,v in master_list.iteritems():
        
        try:
            scan = sarpy.Scan(v[data_label][0])
            sdir = scan.shortdirname
            sdir = sdir.replace('/','_')
            
            if ioType == 'export':
                fname = os.path.join(path, sdir + '.nii')
                scan.pdata[0].export2nii(fname)
            
            elif ioType == 'import':
                
                if adata_label is None and forceVal is True:
                    raise IOError(
                    'Please specify an adata_label or remove the force flag')
                    
                elif adata_label is None:
                    adata_label = 'roi'
                    print("saving roi in adata with generic 'roi' label")                    
                
                    fname = os.path.join(path, sdir + 'a.nii')
                    roi = nibabel.load(fname).get_data()[:,:,:,0]
                    scan.store_adata(key=adata_label, data = roi, force = forceVal)

                else:
                    print("Please specify either 'import' or 'export' \
                    for the ioType!")
            
        except IOError:
            
            print('Not found: {0} and {1}'.format(k,data_label) )
            
    print('Nifti images were processed in {0}'.format(path))


ms = 'NecS3'

# First export the ROIs as nifti files
generate_ROI(ms, '0h-IR_A', ioType = 'export')
generate_ROI(ms, '24h-IR_A', ioType = 'export')

# Next, go in to ImageJ and draw you ROIs, when you save the file, append 'a'
# to the end of the file name so your ROIs don't accidentally get overwritten
# This routine requires that there is an 'a' at the end of the filename. Leave
# the rest of the file name as exactly the same
generate_ROI(ms, '0h-IR_A', ioType = 'import', forceVal = True)
generate_ROI(ms, '24h-IR_A', ioType = 'import', forceVal = True)