# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 13:27:34 2013

# This is a comment
# I am trying to figure out how to use python

@author: fmoosvi
"""

import SARlabpy as b
# import os as os
import pylab as p

# To get current working directory, make sure os is imported, then do os.chdir(path) or os.getcwd()

path = '/Volumes/Data/brukerdata-rsync/stefan/nmr/dBlochSiegert1.gP2/38/pdata/1/'

dataDict = b.read2dseq(path)

data = dataDict['data']


data1=b.fftbruker(data,encoding=[1, 1, 1, 0])

p.imshow(data1.real[:,:,30])

data1[10,10,10]