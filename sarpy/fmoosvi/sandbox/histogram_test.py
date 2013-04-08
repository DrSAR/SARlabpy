# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 00:57:52 2013

@author: firas
"""
import numpy
import sarpy
import pylab
import scipy

pylab.close('all')

# Trying to display data correctly

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
fig_a = fig.add_subplot(221)
scn_i = pylab.imshow(scn[:,:,2], cmap = 'gray'), pylab.colorbar()
fig_a.set_title('Default')
##### Method 1: Histogram Normalization

histo_eq,b = histeq(scn, cnum=255)
fig_b = fig.add_subplot(222)
a_i = pylab.imshow(histo_eq[:,:,2], cmap = 'gray')
fig_b.set_title('Histogram Equalization')

pylab.colorbar()

##### Method 2: 99.99th Percentile

fig_c = fig.add_subplot(223)
a_ii = pylab.imshow(scn[:,:,2], cmap = 'gray')
fig_c.set_title('99.9th percentile')
lim_low = scipy.percentile(scn,10)
lim_high = scipy.percentile(scn,99.9)
a_ii.set_clim(lim_low, lim_high)
pylab.colorbar()

##### Method 3: Ideal

fig_b = fig.add_subplot(224)
a_iii = pylab.imshow(scn[:,:,2], cmap = 'gray')
fig_b.set_title('Ideal Image')
a_iii.set_clim(1E4, 5E5)
pylab.colorbar()
fig.tight_layout

#### Checking Param maps ########

scn1 = sarpy.Experiment('NecS3').find_scan_by_protocol('04')[5].adata['T1map_LL'].data.get_data()
fig1 = pylab.figure()
fig1_a = fig1.add_subplot(221)
scn1_i = pylab.imshow(scn1[0,:,:,2])
pylab.colorbar()
fig1_a.set_title('Default')

##### Method 1: Histogram Normalization

histo1_eq,b = histeq(scn1, cnum=255)
fig1_b = fig1.add_subplot(222)
a1_i = pylab.imshow(histo1_eq[0,:,:,2])
fig1_b.set_title('Histogram Equalization')

pylab.colorbar()

##### Method 2: 99.99th Percentile

fig1_c = fig1.add_subplot(223)
a1_ii = pylab.imshow(scn1[0,:,:,2])
fig1_c.set_title('97th percentile')
lim_low = scipy.percentile(scn1,65)
lim_high = scipy.percentile(scn1,98)
a1_ii.set_clim(lim_low, lim_high)
pylab.colorbar()

##### Method 3: Ideal

fig1_b = fig1.add_subplot(224)
a1_iii = pylab.imshow(scn1[0,:,:,2])
fig1_b.set_title('Ideal Image')
a1_iii.set_clim(800, 3000)
pylab.colorbar()
fig1.tight_layout