# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 00:24:57 2013

@author: firas
"""

# testing code

import numpy
import os
import SARlabpy as sar
import time
import pylab as py

reload(sar)

dir_name = os.path.expanduser('~/data/NecS1Hs02.hi1/8/pdata/1/');

data_dict = sar.read2dseq(dir_name)

sar.enhancementCurve(data_dict)

   
    py.imshow(data_dict['data'][4,:,:,50])


