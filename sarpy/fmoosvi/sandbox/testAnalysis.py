# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 00:24:57 2013

@author: firas
"""
# Add sarpy to the di

import numpy
import os
import sarpy
import pylab

reload(sarpy)


NecS1Exp = sarpy.Experiment('NecS1')

NecS1dce = NecS1Exp.studies[2].scans[6]

auc = sarpy.analysis.calculate_AUC(NecS1dce)
#xlim( (0,20) )

fig1 = pylab.figure();
pylab.imshow(auc[:,:,3])
pylab.axis('off')
pylab.colorbar()
pylab.clim(5,60)
pylab.title('Sample tumour AUC')


## Creating an AUC template data set

#fig3 = pylab.figure()










