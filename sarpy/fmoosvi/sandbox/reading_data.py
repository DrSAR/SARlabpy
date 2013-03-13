# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 22:51:48 2013

@author: firas
"""

import sarpy
import numpy
import pylab
import scipy

readfidExp = sarpy.Experiment('readfid')
for study in readfidExp.studies:
    print('-'*40+'\n'+study.subject.SUBJECT_id)
    for scan in study.scans:
        print("  "+scan.acqp.ACQ_protocol_name)

#NecS3Exp= sarpy.Experiment('NecS3')
#for study in NecS3Exp.studies:
#    print('-'*40+'\n'+study.subject.SUBJECT_id)
#    for scan in study.scans:
#        print("  "+scan.acqp.ACQ_protocol_name)


NecS3Exp= sarpy.Experiment('NecS3')
LL_scans = sarpy.Experiment.find_scan_in_experiment(NecS3Exp, '04_ubcLL2')

LL = LL_scans[0]
LLdata = LL_scans[0].pdata[0].data

x = test.pdata[0].data.shape[0]
y = test.pdata[0].data.shape[1]
num_slices = test.pdata[0].data.shape[2]

pts = []

for i in range(25):
    pts.append(numpy.mean(LLdata[20:40,40:50,5,i]))
    
pylab.plot(scipy.array(pts))