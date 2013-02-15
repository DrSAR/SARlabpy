# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 13:27:34 2013

# This is a comment
# I am trying to figure out how to use python

@author: fmoosvi
"""

import SARlabpy as sar
# import os as os
import pylab as py

# To get current working directory, make sure os is imported, then do os.chdir(path) or os.getcwd()

path = '~/dBlochSiegert1.gP2/38/pdata/1/'
path1 = '~/Moosvi.ii2/10/pdata/1/'

dataDict = sar.read2dseq(path)
dataDict2 = sar.read2dseq(path1)

data = dataDict['data']
data2 = dataDict['data']


py.imshow(data1.real[:,:,30])

print data1[10,10,10]

#data1=b.fftbruker(data,encoding=[1, 1, 1, 0])
