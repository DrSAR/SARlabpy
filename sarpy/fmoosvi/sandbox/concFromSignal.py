# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 17:11:44 2013

@author: fmoosvi
"""

#conc from signal


import numpy
import sarpy
import pylab
import math
import scipy
import sarpy.fmoosvi.analysis
import sarpy.fmoosvi.getters as getters
import sarpy.ImageProcessing.resample_onto





n = sarpy.Experiment('NecS3').studies[13].find_scan_by_protocol('06_')[0]
     
t1map = sarpy.Experiment('NecS3').studies[13].find_scan_by_protocol('04_')[0]
     
rt1 = h_conc_from_signal(n,t1map,relaxivity=1)

augc = h_calculate_AUGC(n,rt1)

#f = pylab.figure()
#
#a = pylab.imshow(rt1[:,:,3,10])
#a.set_clim(0, 2E-3)
#f.canvas.draw()
#
#g = pylab.figure()
#
#b = pylab.imshow(rt1[:,:,3,50])
#b.set_clim(0, 2E-3)
#g.canvas.draw()
#
#h = pylab.figure()
#
#c = pylab.imshow(rt1[:,:,3,90])
#c.set_clim(0, 2E-3)
#h.canvas.draw()