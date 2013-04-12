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

## Scenario 2: LL and RARE
necs3 = sarpy.Experiment('NecS3').studies[7]
ll_scan = necs3.find_scan_by_protocol('04')[0].pdata[0]
roi_scan = necs3.find_scan_by_protocol('05')[0].pdata[0]

fig2 = pylab.figure()

res2 = sarpy.ImageProcessing.resample_onto.resample_onto_pdata(roi_scan,dce_scan)

pylab.imshow(res2[:,:,0]) #results in all nans!

## Scenario 3: LL and RARE adata (ROI)
necs3 = sarpy.Experiment('NecS3').studies[7]
ll_scan = necs3.find_scan_by_protocol('04')[0].pdata[0]
roi_scan = necs3.find_scan_by_protocol('05')[0].adata['IR_tumour_rois']

fig2 = pylab.figure()

res3= sarpy.ImageProcessing.resample_onto.resample_onto_pdata(roi_scan,dce_scan)

pylab.imshow(res3[:,:,3]) #results in all nans!



