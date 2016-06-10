# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 00:27:53 2013

@author: firas
"""

import sarpy
import sarpy.fmoosvi.analysis
import sarpy.fmoosvi.getters as getters
import numpy
import pylab
import math
import scipy

pylab.close('all')

NecS3 = sarpy.Experiment('NecS3').find_scan_by_protocol('04')

scan_object = NecS3[0]
flip_angle_map = math.radians(scan_object.acqp.ACQ_flip_angle)


params = [20000, 2, 650]
pdata_num = 0

num_slices = getters.get_num_slices(scan_object, pdata_num)
repetition_time = scan_object.method.PVM_RepetitionTime
inversion_time = scan_object.method.PVM_InversionTime
data = scan_object.pdata[pdata_num].data[:]

t_data = numpy.linspace(inversion_time,\
    scan_object.pdata[pdata_num].data.shape[3]*repetition_time,\
    scan_object.pdata[pdata_num].data.shape[3])
    
data_after_fitting = numpy.zeros([1,1,1])
fit_results = numpy.array(data_after_fitting[:], dtype=dict)                       

point = data[63,100,2,:]
fit_dict = {}


fit_params,cov,infodict,mesg,ier = scipy.optimize.leastsq(sarpy.fmoosvi.analysis.h_residual_T1,
                                                          params,args=(point,t_data), full_output=True, \
                                                          maxfev = 200)


## Data 
print fit_params
fit_data = sarpy.fmoosvi.analysis.h_func_T1(fit_params,t_data)
#fir_data = sarpy.fmoosvi.analysis.h_func_T1([1,1,1],t_data)

[M,B,T1] = fit_params

fit_dict = {
            'fit_params': fit_params,
            'cov' : cov,
            'infodict' : infodict,
            'mesg' : mesg,
            'ier' : ier
            }


ss_err=(infodict['fvec']**2).sum()
ss_tot=((point-point.mean())**2).sum()
rsquared=1-(ss_err/ss_tot)
print 'R2', rsquared



pylab.plot(t_data,point,'x')
pylab.plot(t_data,fit_data,'-')

#def h_fit_T1_LL(scan_object, flip_angle_map = 0, pdata_num = 0, 
#                params = []):
#    
#    if len(params) == 0:      
#        params = [3E5, 2, 350]
#
#    num_slices = getters.get_num_slices(scan_object, pdata_num)
#    repetition_time = scan_object.method.PVM_RepetitionTime
#    inversion_time = scan_object.method.PVM_InversionTime
#    
#    if type(flip_angle_map) != numpy.ndarray:
#        flip_angle_map = math.radians(scan_object.acqp.ACQ_flip_angle)
#   
#    # Visu_pars params 
#        
#    data = scan_object.pdata[pdata_num].data[:]
#
#    data_after_fitting = numpy.zeros( [data.shape[0],\
#                                       data.shape[1],\
#                                       data.shape[2]] )
#                                       
#    fit_results = numpy.array(data_after_fitting[:], dtype=dict)                       
#            
#    t_data = numpy.linspace(inversion_time,\
#        scan_object.pdata[pdata_num].data.shape[3]*repetition_time,\
#        scan_object.pdata[pdata_num].data.shape[3])
#                                                       
##    fitfunc = lambda param, x_data: numpy.abs(param[0]*\
##                                (1 - param[1]*numpy.exp(-t_data/param[2])))
#  
#    for x in xrange(data.shape[0]):
#        for y in range(data.shape[1]):
#            for slice in range(num_slices):
#                
#                y_data = data[x,y,slice,:]
#                fit_dict = {}
#                
#                fit_params,cov,infodict,mesg,ier = scipy.optimize.leastsq( h_residual_T1,params,args=(y_data,t_data), full_output=True)
#
#                [M,B,T1] = fit_params
#                fit_dict = {
#                            'fit_params': fit_params,
#                            'cov' : cov,
#                            'infodict' : infodict,
#                            'mesg' : mesg,
#                            'ier' : ier
#                            }
#                            
#                data_after_fitting[x,y,slice] = T1
#                fit_results[x,y,slice] = fit_dict
#                