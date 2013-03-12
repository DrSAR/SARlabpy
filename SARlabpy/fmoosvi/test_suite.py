# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 14:48:36 2013

Test Suite

@author: fmoosvi
"""
import numpy
import pylab
import SARlabpy as sar
import pdb
import numpy.testing
import os


NecS1dce = cls.Experiment('NecS1').studies[2].scans[6]

# Test functions
test_calculate_AUC(NecS1dce)


def test_calculate_AUC(scan_object, pdata_num = 0, debug = False):
    
    # Method params
    scan_object.method.PVM_RepetitionTime = 1000/50
    scan_object.method.PVM_NRepetitions = 70
    scan_object.method.PVM_EncMatrix = 50
    
    # Data altering
    scan_object.pdata[pdata_num].data = numpy.ones([100,50,6,70])
    scan_object.pdata[pdata_num].data[50:78,22:42,0:6,10] = 7
    scan_object.pdata[pdata_num].data[50:78,22:42,0:6,10:110] = numpy.arange(7,1,-0.1)

    norm_data = sar.analysis.normalize_dce(scan_object)
    inj_point = sar.analysis.inj_point(scan_object)
    auc_template = sar.analysis.calculate_AUC(scan_object,60)
    
    pylab.plot(norm_data[60,35,0,0:69])   
    # Compute the AUC of each slice using triangle area formula
    base = data_dict['header']['method']['PVM_NRepetitions'] - (inj_point)
    height = numpy.round(numpy.max(norm_data[:,:,:,:]),2)
    
    area_triangle = base*height/2
    
    try:
        numpy.testing.assert_almost_equal(area_triangle, numpy.max(auc_template), decimal=3, verbose=True)

    except AssertionError:
        fig = pylab.figure
        pylab.plot(norm_data[60,30,0,:])    
        
