# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 00:57:52 2013

@author: firas
"""
import numpy
import sarpy
import pylab
import scipy
#displaying data correctly





def histeq(image_data,bins=256, cnum= 255):

   #get image histogram
   imhist,bins = numpy.histogram(image_data.flatten(),bins,normed=True)
   cdf = imhist.cumsum() #cumulative distribution function
   cdf = cnum * cdf / cdf[-1] #normalize

   #use linear interpolation of cdf to find new pixel values
   im2 = numpy.interp(image_data.flatten(),bins[:-1],cdf)

   return im2.reshape(image_data.shape), cdf


scn = sarpy.Experiment('NecS3').find_scan_by_protocol('05')[0].pdata[0].data
fig = pylab.figure()
fig_a = fig.add_subplot(131)
scn_i = pylab.imshow(scn[:,:,3], cmap = 'gray'), pylab.colorbar()

##### Method 1: Histogram Normalization

histo_eq,b = histeq(scn, cnum=255)
fig_b = fig.add_subplot(132)
a_i = pylab.imshow(histo_eq[:,:,3], cmap = 'gray')
pylab.colorbar()

##### Method 2: 99.99th Percentile

fig_b = fig.add_subplot(133)
a_ii = pylab.imshow(scn[:,:,3], cmap = 'gray')

lim_low = scipy.percentile(scn,10)
lim_high = scipy.percentile(scn,99.99)
a_ii.set_clim(lim_low, lim_high)
pylab.colorbar()








############
