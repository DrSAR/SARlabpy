# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 21:05:50 2013

@author: firas
"""

## Testing function to take in a mask and compute the distribution of T1s in that mask.

import sarpy
import nibabel
import numpy
import pylab
import scipy

NecS3Hs04a = sarpy.Experiment('NecS3').studies[0].find_scan_by_protocol('04')[0]
NecS3Hs04b = sarpy.Experiment('NecS3').studies[0].find_scan_by_protocol('04')[0]


scan1 = NecS3Hs04a.pdata[0] 
scan2 = NecS3Hs04b.pdata[0]

#scan1.write2nii('scan1.nii')
#scan2.write2nii('scan2.nii')

scan1mask = nibabel.load('scan1ed.nii').get_data()[:,:,:,0]
scan2mask = nibabel.load('scan2ed.nii').get_data()[:,:,:,0]


## magic getting T1 of masked region

def image_to_mask(data):

    mask_val = scipy.percentile(data.flatten(),95)
    data[data == mask_val] = numpy.nan
    data[data != numpy.nan] = 1

    return data
    
scan1masked= image_to_mask(scan1mask)

a =scan1masked * scan1.data