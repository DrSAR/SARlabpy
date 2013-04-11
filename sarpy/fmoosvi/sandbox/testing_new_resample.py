# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 01:59:31 2013

@author: firas
"""

## Test new resample onto

import sarpy
import nibabel
import sarpy.ImageProcessing.resample_onto
import os
import pylab

pylab.close('all')

## Scenario 0: Example Scenario

source_pdata = sarpy.Scan("readfidTest.ix1/7").pdata[0]
print source_pdata.data.shape
#        (105, 133, 5, 25)
target_pdata = sarpy.Scan("readfidTest.ix1/3").pdata[0]
print target_pdata.data.shape
#        (105, 133, 25)
res0 = sarpy.ImageProcessing.resample_onto.resample_onto_pdata(source_pdata, target_pdata)
print res0.shape
#        (105, 133, 25, 25)
fig0 = pylab.figure()
pylab.imshow(res0[:,:,0,10])


## Scenario 1: dce and rare both pdata
necs3 = sarpy.Experiment('NecS3').studies[7]
dce_scan = necs3.find_scan_by_protocol('06')[0].pdata[0]
roi_scan = necs3.find_scan_by_protocol('05')[0].pdata[0]

adata_key = 'IR_tumour_rois'

res = sarpy.ImageProcessing.resample_onto.resample_onto_pdata(roi_scan,dce_scan)

fig1= pylab.figure()
pylab.imshow(res[:,:,0]) #results in all nans!

## Scenario 2: LL and RARE
necs3 = sarpy.Experiment('NecS3').studies[7]
ll_scan = necs3.find_scan_by_protocol('04')[0].pdata[0]
roi_scan = necs3.find_scan_by_protocol('05')[0].pdata[0]

fig2 = pylab.figure()

res2 = sarpy.ImageProcessing.resample_onto.resample_onto_pdata(roi_scan,dce_scan)

pylab.imshow(res2[:,:,0]) #results in all nans!

