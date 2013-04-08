# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 20:58:56 2013

@author: fmoosvi
"""

import sarpy
import nibabel
import sarpy.ImageProcessing.resample_onto
import os
import pylab

pylab.close('all')

#necs3 = sarpy.Experiment('NecS3').studies[6]
#IR_scans = necs3.find_scan_by_protocol('05')
#scan = IR_scans[0]

################################
# Testing daughter_roi and parent_image
################################
#scan_image1 = scan.pdata[0].export2nii('image_test.nii')
scan_image = nibabel.load('image_test.nii')

#scan_roi1 = scan.adata['IR_tumour_rois'].export2nii('roi_test.nii')
scan_roi = nibabel.load('roi_test.nii')

img_path = os.path.abspath('image_test.nii')
roi_path = os.path.abspath('roi_test.nii')


result = sarpy.ImageProcessing.resample_onto.resample_onto(img_path,roi_path)

fig = pylab.figure()

fig.add_subplot(311)
pylab.imshow(scan_image.get_data()[:,:,3])

fig.add_subplot(312)
pylab.imshow(scan_roi.get_data()[:,:,3,0])

fig.add_subplot(313)
pylab.imshow(result[:,:,3])