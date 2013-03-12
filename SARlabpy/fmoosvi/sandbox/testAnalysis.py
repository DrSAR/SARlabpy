# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 00:24:57 2013

@author: firas
"""
# Add sarlabpy to the di

import sys
#sys.path.append("/Users/firas/git/SARlabpy/")
#sys.path.append('/Volumes/Data/Dropboxes/PhD./Dropbox/code/python/')

# testing code

import numpy
import os
import SARlabpy as sar
import pylab

reload(sar)

dir_name = os.path.expanduser('~/data/NecS1Hs02.hi1/8/pdata/1/');
data_dict = sar.read2dseq(dir_name)

fig1 = pylab.figure();
auc = sar.analysis.calculate_AUC(data_dict)
#xlim( (0,20) )

fig2 = pylab.figure();
pylab.imshow(auc[:,:,3])
pylab.axis('off')
pylab.colorbar()
pylab.clim(5,60)
pylab.title('Sample tumour AUC')


## Creating an AUC template data set

fig3 = pylab.figure()










