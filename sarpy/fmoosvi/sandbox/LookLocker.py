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
import pylab
import time

phantom = sarpy.Experiment('DiLL')

LLdata = phantom.find_scan_by_protocol('MOBILE*')

start_time = time.time()
proc_data = sarpy.fmoosvi.analysis.h_fit_T1(LLdata)

print time.time() - start_time

img = pylab.imshow(proc_data[:,:,0], aspect = 0.5)
img.set_clim(0.0,1200)
pylab.colorbar()
pylab.show()


print numpy.mean(proc_data[45:80,10:20,0])
print numpy.mean(proc_data[45:52,35:55,0])
