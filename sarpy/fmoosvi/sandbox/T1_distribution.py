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
import sarpy.fmoosvi.wrappers

NecS3Hs04a = sarpy.Experiment('NecS3').studies[0].find_scan_by_protocol('04')[0]
NecS3Hs04b = sarpy.Experiment('NecS3').studies[0].find_scan_by_protocol('04')[0]


scan1 = NecS3Hs04a.adata['T1map_LL'].data.get_data()[0,:,:,:]
scan2 = NecS3Hs04b.adata['T1map_LL'].data.get_data()[0,:,:,:]

#scan1.export2nii('scan1.nii')
#scan2.export2nii('scan2.nii')

scan1mask = nibabel.load('scan1ed.nii').get_data()[:,:,:,0]
scan2mask = nibabel.load('scan2ed.nii').get_data()[:,:,:,0]





 # magic getting T1 of masked region
 
def image_to_mask(data):

    masked_data = data[:] 
    mask_val = scipy.percentile(data.flatten(),95)
    masked_data[masked_data == mask_val] = numpy.nan
    masked_data[numpy.isfinite(masked_data)] = 1
    return masked_data
     
scan1masked= image_to_mask(scan1mask)
 
data_masked = scan1masked * scan1[:,:,:]

####
fig1 = pylab.figure()
fig1.add_subplot(311)
im1 = pylab.imshow(scan1[:,:,0])
im1.set_clim([600, 2500])
pylab.title('Full T1 Map')
pylab.colorbar()
####
fig1.add_subplot(312)
im2 = pylab.imshow(scan1masked[:,:,0])
pylab.colorbar()
pylab.title('Mask')
####
fig1.add_subplot(313)
im3 = pylab.imshow(data_masked[:,:,0])
im3.set_clim([600, 2500])
pylab.colorbar()
pylab.title('Mask applied to T1map')

## Time to put up a histogram of the T1 distribution

sarpy.fmoosvi.wrappers.roi_distribution(scan1, scan1mask, bins = numpy.linspace(600,2500,30), display_histogram = True)

scan2masked= image_to_mask(scan1mask)
data2_masked = scan2masked * scan2[:,:,:]


####
fig1 = pylab.figure()
fig1.add_subplot(311)
im1 = pylab.imshow(scan2[:,:,0])
im1.set_clim([600, 2500])
pylab.title('Full T1 Map')
pylab.colorbar()
####
fig1.add_subplot(312)
im2 = pylab.imshow(scan2masked[:,:,0])
pylab.colorbar()
pylab.title('Mask')
####
fig1.add_subplot(313)
im3 = pylab.imshow(data2_masked[:,:,0])
im3.set_clim([600, 2500])
pylab.colorbar()
pylab.title('Mask applied to T1map')

## Time to put up a histogram of the T1 distribution

sarpy.fmoosvi.wrappers.roi_distribution(scan2, scan2mask, bins = numpy.linspace(600,2500,30), display_histogram = True)