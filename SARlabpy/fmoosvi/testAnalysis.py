# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 00:24:57 2013

@author: firas
"""
# Add sarlabpy to the di

import sys
sys.path.append("/Users/firas/git/SARlabpy/")


# testing code

import numpy
import os
import SARlabpy as sar
import pylab

reload(sar)

dir_name = os.path.expanduser('~/data/NecS1Hs02.hi1/8/');
scandirname = os.path.expanduser('~/data/NecS1Hs02.hi1/');

data_dict = sar.read2dseq(dir_name)

fig1 = pylab.figure();
auc = sar.calculateAUC(data_dict)
#xlim( (0,20) )
cur = sar.enhancementCurve(data_dict)

fig2 = pylab.figure();
pylab.imshow(auc[:,:,3])
pylab.axis('off')
pylab.colorbar()
pylab.clim(5,30)
pylab.title('Sample tumour AUC')


## Creating an AUC template data set

fig3 = pylab.figure()

data_dict_template = data_dict

data_dict_template['data'] = numpy.ones([128,64,6,110])

filler = numpy.arange

data_dict_template['data'][50:78,22:42,0:6,10] = 11
data_dict_template['data'][50:78,22:42,0:6,10:110] = numpy.arange(11,1,-0.1)


data_dict_template['header']['method']['PVM_RepetitionTime'] = 1000/64
data_dict_template['header']['method']['PVM_NRepetitions'] = 110

inj = sar.enhancementCurve(data_dict_template)
inj_point = sar.determineInjectionPoint(data_dict_template)
auc_template = sar.calculateAUC(data_dict_template)

fig4 = pylab.figure()
pylab.imshow(auc_template[:,:,0], cmap = 'copper')
pylab.colorbar()

# Compute the AUC of each slice as follows:

# Area_triangle = Base*height / 2

base = data_dict_template['header']['method']['PVM_NRepetitions'] - (inj_point+1)
height = numpy.round(inj[0][inj_point+1],2)

area_triangle = base*height/2


print 'The area of the triangle should be 84.0 and it is actually ' + str(area_triangle)

