# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 17:11:44 2013

@author: fmoosvi
"""

#conc from signal


import numpy
import sarpy
import sarpy.fmoosvi.analysis
import sarpy.fmoosvi.getters as getters
import sarpy.ImageProcessing.resample_onto

def h_conc_from_signal(scan_object, scan_object_LL, 
                       adata_label = 'T1map_LL', relaxivity=4.3e-3, 
                       bbox = None, pdata_num = 0):

    ########### Getting and defining parameters
    
    # Data
    data = scan_object.pdata[pdata_num].data
    
    # resample the t1map onto the dce
    data_t1map_pre = scan_object_LL.adata[adata_label]
    
    #data_t1map = sarpy.ImageProcessing.resample_onto.resample_onto_pdata(data_t1map_pre,scan_object.pdata[pdata_num])

    x_size = data.shape[0]
    y_size = data.shape[1]
    num_slices = getters.get_num_slices(scan_object,pdata_num)
    
    # Method params
    #TODO: change this so it doesn't require method file WIHOUT BREAKING IT!
    reps =  scan_object.method.PVM_NRepetitions
    
    TR = scan_object.method.PVM_RepetitionTime
    FA = math.radians(scan_object.acqp.ACQ_flip_angle)
    
    inj_point = sarpy.fmoosvi.analysis.h_inj_point(scan_object, pdata_num = 0)    
    
    # Deal with bounding boxes

    if bbox is None:        
        bbox = numpy.array([0,x_size-1,0,y_size-1])    
       
    if bbox.shape == (4,):            
    
        bbox_mask = numpy.empty([x_size,y_size])
        bbox_mask[:] = numpy.nan        
        bbox_mask[bbox[0]:bbox[1],bbox[2]:bbox[3]] = 1
    
        # First tile for slice
        bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size,1),num_slices)

    else:      
        raise ValueError('Please supply a bbox for h_conc_from_signal')   

    T1 = numpy.empty([x_size,y_size,num_slices,reps])
    T1[:] = numpy.nan
    
    for x in xrange(bbox[0],bbox[1]):
        for y in xrange(bbox[2],bbox[3]):
            for slice in xrange(num_slices):
                               
                baseline_s = numpy.mean(data[x,y,slice,0:inj_point])
                
                E1 = numpy.exp(-TR/data_t1map_pre.data[x,y,slice])
                c = numpy.cos(FA)
                
                T1[x,y,slice,0:inj_point = data_t1map_pre.data[x,y,slice]
                               
                for rep in xrange(inj_point,reps):
                    
                    s = data[x,y,slice,rep] / baseline_s
                    E2 = (-E1*c + E1*s - s + 1) / (E1*s*c - E1*c - s*c +1)
                    
                    T1[x,y,slice,rep] = -TR / numpy.log(E2)
    
    T1baseline = data_t1map_pre.data*bbox_mask
    T1baseline =  numpy.tile(T1baseline.reshape(x_size,y_size,num_slices,1),reps)
    conc = (1/relaxivity) * ( (1/T1) - (1/T1baseline) )
                
    return conc
                


n = sarpy.Experiment('NecS3').studies[13].find_scan_by_protocol('06_')[0]
     
t1map = sarpy.Experiment('NecS3').studies[13].find_scan_by_protocol('04_')[0]
     
rt1 = h_conc_from_signal(n,t1map,relaxivity=1)

f = pylab.figure()

a = pylab.imshow(rt1[:,:,3,10])
a.set_clim(0, 2E-3)
f.canvas.draw()

g = pylab.figure()

b = pylab.imshow(rt1[:,:,3,50])
b.set_clim(0, 2E-3)
g.canvas.draw()

h = pylab.figure()

c = pylab.imshow(rt1[:,:,3,90])
c.set_clim(0, 2E-3)
h.canvas.draw()