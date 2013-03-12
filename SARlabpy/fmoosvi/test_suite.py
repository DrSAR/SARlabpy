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

dir_name = os.path.expanduser('~/data/NecS1Hs02.hi1/8/');
data_dict = sar.read2dseq(dir_name)

# Test the function
test_calculate_AUC(data_dict)




def test_calculate_AUC(data_dict, debug = False):
     
    data_dict['data'] = numpy.ones([100,50,6,70])
    data_dict['data'][50:78,22:42,0:6,10] = 7
    data_dict['data'][50:78,22:42,0:6,10:110] = numpy.arange(7,1,-0.1)
    data_dict['header']['method']['PVM_RepetitionTime'] = 1000/50
    data_dict['header']['method']['PVM_NRepetitions'] = 70
    data_dict['header']['method']['PVM_EncMatrix'][1] = 50
    
    norm_data = sar.analysis.normalize_dce(data_dict)
    inj_point = sar.analysis.inj_point(data_dict)
    auc_template = sar.analysis.calculate_AUC(data_dict,60)
    
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
        
