# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 20:58:56 2013

@author: fmoosvi
"""

import sarpy

import sarpy.ImageProcessing

necs3 = sarpy.Experiment('NecS3')
IR_scans = necs3.find_scan_by_protocol('05')
scan = IR_scans[6]

################################
# Testing daughter_roi and parent_image
################################
scan_img = IR_scans[6].pdata[0].data
scan_roi = scan.adata['IR_tumour_rois'].data

result = sarpy.ImageProcessing.rs(scan_img,scan_roi)