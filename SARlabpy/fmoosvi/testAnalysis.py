# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 00:24:57 2013

@author: firas
"""

# testing code

import numpy
import os
import SARlabpy as sar
import pylab
import operator

reload(sar)

dir_name = os.path.expanduser('~/data/NecS1Hs02.hi1/8/');
scandirname = os.path.expanduser('~/data/NecS1Hs02.hi1/');

data_dict = sar.read2dseq(dir_name)

fig1 = pylab.figure();
sar.calculateAUC(scandirname,'8')

cur = sar.enhancementCurve(data_dict)


fig2 = pylab.figure();
pylab.imshow(data_dict['data'][4,:,:,50])

##

a,b = sar.determineInjectionPoint(data_dict)

#fig3 = pylab.figure()



