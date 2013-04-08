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

necs3 = sarpy.Experiment('NecS3').studies[6]
IR_scans = necs3.find_scan_by_protocol('05')
scan = IR_scans[0]

################################
# Testing daughter_roi and parent_image
################################
scan_image1 = scan.pdata[0].export2nii('image_test.nii')
scan_image = nibabel.load('image_test.nii')

scan_roi1 = scan.adata['IR_tumour_rois'].export2nii('roi_test.nii')
scan_roi = nibabel.load('roi_test.nii')

img_path = os.path.abspath('image_test.nii')
roi_path = os.path.abspath('roi_test.nii')


result = sarpy.ImageProcessing.resample_onto.resample_onto(img_path,roi_path)

fig = pylab.figure()

fig.add_subplot(311)
pylab.imshow(scan_image.get_data()[:,:,3])

fig.add_subplot(312)
pylab.imshow(scan_roi.get_data()[:,:,3])

fig.add_subplot(313)
pylab.imshow(result[:,:,3])

################################
# Testing daughter_roi and parent_image
################################

phantomOrien3D = sarpy.Experiment('Phantom').studies[0].scans[-1].pdata[0]
phantomOrienAx = sarpy.Experiment('Phantom').studies[0].scans[1].pdata[0]

phantomOrien3D.export2nii('p3D.nii')
phantomOrienAx.export2nii('pax.nii')




result1= sarpy.ImageProcessing.resample_onto.resample_onto('p3D.nii','pax.nii')

result2 = sarpy.ImageProcessing.resample_onto.resample_onto('pax.nii','p3D.nii')

fig = pylab.figure()

#fig.add_subplot(311)
#pylab.imshow(phantomOrien3D.data[:,:,3])

fig.add_subplot(221)
pylab.imshow(phantomOrienAx.data[:,:,3])

fig.add_subplot(222)
pylab.imshow(result1[:,:,3])

fig.add_subplot(223)
pylab.imshow(phantomOrien3D.data[:,75,:])

fig.add_subplot(224)
pylab.imshow(result2[:,75,:])



###########
necs3 = sarpy.Experiment('NecS3Hs04')
imp = necs3.studies[-1]

scan1 = imp.scans[7]
scan2 = imp.scans[8]

scan1.pdata[0].export2nii('scan1.nii')
scan2.pdata[0].export2nii('scan2.nii')

result3 = sarpy.ImageProcessing.resample_onto.resample_onto('scan1.nii','scan2.nii')

fig.add_subplot(221)
pylab.imshow(scan1.data[:,:,3,0])

fig.add_subplot(222)
pylab.imshow(scan2.data[:,:,3,0])

fig.add_subplot(223)
pylab.imshow(result3[:,:,3])