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
path2 = '/Users/fmoosvi/data/Moosvi.ii2/22/pdata/1/'

#dataDict = sar.read2dseq(path,1)
dataDict2 = sar.read2dseq(path2,1)

#data = dataDict['data']
data2 = dataDict2['data']

#py.imshow(data2[:,:,2])

#imgplot = py.imshow(data2[:,:,2],aspect=1)
#imgplot.set_clim(2E3,5E4)


## Editing jcamp to split up parameters stored in the ExcPulse Array

filename_method = ''.join(['/Users/fmoosvi/data/Moosvi.ii2/22/','acqp'])
#filename_method = ''.join(['/Volumes/Data/work/dBS1/dBlochSiegert1.gP2/36/','acqp'])

method = sar.readJCAMP(filename_method)
methodO = sar.readJCAMP(filename_method)
#newDict = sar.io.BRUKERIO.splitParamArrays(method,'ExcP1ulse')

## Write a function that converts strings to intgers/floats or strings as appropriate

## Begin section of code dealing with various cases in the acqp or method file
#
# Purely one integer, or one float have already been dealt with
#
# Case 1: array of ints. 'ACQ_phase_enc_start': ' -1 -1',
# Case 2: array of floats and ints. 'ACQ_spatial_phase_1': ' -1 -0.9583 -0.9166 -0.875 -0.8333 ...
# Case 3: weird combo. 'TPQQ': ' (<hermite.exc>, 16.4645986123031, 0) (<fermi.exc>, 115.8030276379, 0)
#
##
import re

method = sar.io.BRUKERIO.typecastThings(method) # Case 0


# Cases 1-3

for dict_item in method:

    try:

        if isinstance(method[dict_item],str):

            # Case 1: Array of purely integers
            split_string = [s for s in re.split('[(), < >]',method[dict_item]) if s]
            split_string = sar.io.BRUKERIO.typecastThings(split_string)

            if all(isinstance(list_item, int) for list_item in split_string):
                print('you win - integer')
                method[dict_item] = split_string

            # Case 2: Array of floats and ints
            split_string = [s for s in re.split('[(), < >]',method[dict_item]) if s]
            split_string = sar.io.BRUKERIO.typecastThings(split_string)

            if all(isinstance(list_item,(int,float)) for list_item in split_string):
                print('you win - float')
                method[dict_item] = split_string

            # Case 3: Array of floats, ints, and strings with funny brackets and such
            # Example: 'TPQQ': ' (<hermite.exc>, 16.4645986123031, 0) (<fermi.exc>, 115.8030276379, 0)

    except:
        print 'string'







