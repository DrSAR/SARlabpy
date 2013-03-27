# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 15:41:30 2013

@author: fmoosvi
"""
import sarpy
import numpy
import sarpy.fmoosvi.wrappers
import time
import pylab

NecS3_exp = sarpy.Experiment('NecS3')

## AUC Maps

DCE_scans = [study.find_scan_by_protocol('06')[0] for study in NecS3_exp.studies if len(study.find_scan_by_protocol('06'))>0]

auc = sarpy.fmoosvi.wrappers.calculate_AUC(DCE_scans)


aim = pylab.imshow(auc[0][:,:,3])
aim.set_clim([0,80])
pylab.colorbar()
pylab.show()