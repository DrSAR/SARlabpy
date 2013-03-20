# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 22:41:46 2013

@author: firas
"""

## Porting over LookLocker code from Annie

import sarpy
import sarpy.fmoosvi.analysis
import sarpy.fmoosvi.getters
import sarpy.fmoosvi.wrappers
import numpy
import pylab
import time


######################## Start B1 Map #########################
# This is for one phantom study using LL and BSB1 MSME
APExp = sarpy.Experiment('DiLL.iI1')
BS_scans = APExp.find_scan_by_protocol('09\-2DBSB1map\-MSME_bas\(modified\)')[-4:]

LLdata = sarpy.Experiment('DiLL.iI1').find_scan_by_protocol('MOBILE*')

zero_BS_minus = BS_scans[0]
zero_BS_plus = BS_scans[1]
power_BS_minus = BS_scans[2]
power_BS_plus = BS_scans[3]

b1map = sarpy.fmoosvi.analysis.h_BS_B1map(zero_BS_minus,\
                                          zero_BS_plus,\
                                          power_BS_minus,\
                                          power_BS_plus,\
                                          LLdata[0])
                                          
pylab.imshow(b1map[:,:,0,0], aspect = 0.5)
pylab.colorbar()
######################## End B1 Map #########################

######################## Start Look Locker T1 map #########################


proc_data = sarpy.fmoosvi.wrappers.calculate_T1map(LLdata, flip_angle_map = numpy.fliplr(b1map[:,:,:,0]))

start_time = time.time()
#proc_data = sarpy.fmoosvi.analysis.h_fit_T1(LLdata, FA = numpy.fliplr(b1map[:,:,:,0]))
proc_data_noAngle = sarpy.fmoosvi.analysis.h_fit_T1_LL(LLdata, flip_angle_map = 0)

print 'This run took {0} seconds.'.format(round(time.time() - start_time))

#img = pylab.imshow(proc_data[:,:,0], aspect = 0.5)
#img.set_clim(0.0,1200)
#pylab.colorbar()
#pylab.show()
#
#pylab.figure()
#img = pylab.imshow(proc_data_noAngle[:,:,0], aspect = 0.5)
#img.set_clim(0.0,1200)
#pylab.colorbar()
#pylab.show()
#
#print ('-----Without FA correction')
#print numpy.mean(proc_data_noAngle[45:80,10:20,0])
#print numpy.mean(proc_data_noAngle[45:52,35:55,0])
#print ('-----With FA correction')
#print numpy.mean(proc_data[45:80,10:20,0])
#print numpy.mean(proc_data[45:52,35:55,0])
