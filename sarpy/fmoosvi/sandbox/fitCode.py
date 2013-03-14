# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 14:19:16 2013

@author: firas
"""

# Fitting code
from __future__ import division
import pylab
import scipy
import scipy.optimize
import numpy
import sarpy.fmoosvi.analysis
def fit(x_data, y_data, fit_function, initial_params):
    
    
#   Run this to see how fit function works
#        x = numpy.linspace(0, 2000, 2000)
#        y = numpy.abs(1E3 * (1 - 2*numpy.exp(-x/800)))
#        
#        param = [0,0,0]
#        param[0] = 800
#        param[1] = 2
#        param[2] = 1000
#        
#        fitfunc = lambda param, x: numpy.abs(param[0]*(1 - param[1]*numpy.exp(-x/param[2])))
#        initial_params = [10,2,200]
#        
#        param1 = fit(x,y,fitfunc,initial_params)
#        
#        pylab.plot(x, y, "r-", fitfunc(param1,x),"b--")        

    # Distance to the target function    
#    errfunc = lambda param, time, data: fit_function(param,x_data) - y_data

    # Fit function to data    
#    final_params, success = scipy.optimize.leastsq(errfunc, initial_params[:], args = (x_data, y_data))
    
    return final_params

pdata_num = 0
    
scan_object = sarpy.Experiment('NecS3').studies[0].find_scan_by_protocol('04_ubcLL2')[0]
        
data = scan_object.pdata[pdata_num].data
repetition_time = scan_object.method.PVM_RepetitionTime

fed_data = numpy.mean(numpy.mean(data[55:56,35:36,1,:],axis=0),axis=0)
      
x_data = numpy.linspace(1,\
    scan_object.pdata[pdata_num].data.shape[3]*repetition_time,\
    scan_object.pdata[pdata_num].data.shape[3])
                                                           
                                                              
fitfunc = lambda param, x: numpy.abs(param[0]*(1 - param[1]*numpy.exp(-x/param[2])))
                              
new_params = sarpy.fmoosvi.analysis.h_fit(x_data,fed_data, fitfunc, [3.5E5, 2, 600] )

print new_params
        
pylab.plot(x_data, fed_data, "rx", x_data, fitfunc(new_params,x_data),"b--")
    
