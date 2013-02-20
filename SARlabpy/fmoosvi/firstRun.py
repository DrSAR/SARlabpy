# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 13:27:34 2013

# This is a comment
# I am trying to figure out how to use python

@author: fmoosvi
"""

import SARlabpy as sar
#import os as os
import pylab as py

# To get current working directory, make sure os is imported, then do os.chdir(path) or os.getcwd()
reload(sar)

#########
## Controlling for differences in opinion on certain things in Bruker IO
#########
#
#import getpass # Get the current username of the script executer by getpass.user()
#
#user_name = getpass.getuser()
#
#if user_name == 'fmoosvi' or user_name == 'firas':
#    print 'You are operating as Firas, the order of the matrix will be X,Y,Z instead of Z,Y,X'
#    # Other block of code

########
# Running code on other computers
########

# Change this path depending on where the symbolic link to '../brukerdata-rsync/stefan/nmr/' is located

your_path = '/Users/fmoosvi/data/'
study_name = 'dBlochSiegert1.gP2' # could set this to be input by the system
series_num = '/38'
default_reco = '/pdata/1/'

#path = ''.join([your_path,study_name,series_num,default_reco])
#path2 = '/Users/fmoosvi/data/Moosvi.ii2/22/pdata/1/'

path2 = '/Users/firas/data/Moosvi.ii2/22/pdata/1/'


#dataDict = sar.read2dseq(path,1)
dataDict2 = sar.read2dseq(path2,1)

#data = dataDict['data']
data2 = dataDict2['data']

#py.imshow(data2[:,:,2])

#imgplot = py.imshow(data2[:,:,2],aspect=1)
#imgplot.set_clim(2E3,5E4)


## Editing jcamp to split up parameters stored in the ExcPulse Array

#filename_method = ''.join(['/Users/fmoosvi/data/Moosvi.ii2/22/','method'])
filename_method = ''.join(['/Users/firas/data/Moosvi.ii2/22/pdata/1/','reco'])

reco = sar.readJCAMP(filename_method)

filename_method = ''.join(['/Users/fmoosvi/data/Moosvi.ii2/22/','acqp'])
#filename_method = ''.join(['/Volumes/Data/work/dBS1/dBlochSiegert1.gP2/36/','acqp'])

method = sar.readJCAMP(filename_method)
methodO = sar.readJCAMP(filename_method)
#newDict = sar.io.BRUKERIO.splitParamArrays(method,'ExcP1ulse')

method = sar.io.BRUKERIO.typecastThings(method) # Case 0
