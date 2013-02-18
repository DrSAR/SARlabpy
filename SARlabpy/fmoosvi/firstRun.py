# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 13:27:34 2013

# This is a comment
# I am trying to figure out how to use python

@author: fmoosvi
"""

import SARlabpy as sar
# import os as os
#import pylab as py

# To get current working directory, make sure os is imported, then do os.chdir(path) or os.getcwd()


########
# Controlling for differences in opinion on certain things in Bruker IO
########

import getpass # Get the current username of the script executer by getpass.user()

user_name = getpass.getuser()

if user_name == 'fmoosvi' or user_name == 'firas':
	print 'You are operating as Firas, the order of the matrix will be X,Y,Z instead of Z,Y,X'
	# Other block of code



########
# Running code on other computers
########

# Change this path depending on where the symbolic link to '../brukerdata-rsync/stefan/nmr/' is located

#your_path = '/Users/fmoosvi/data/'
#study_name = 'dBlochSiegert1.gP2' # could set this to be input by the system
#series_num = '/38'
#default_reco = '/pdata/1/'

#path = ''.join([your_path,study_name,series_num,default_reco])
#path2 = '/Users/fmoosvi/data/Moosvi.ii2/22/pdata/1/'

#dataDict = sar.read2dseq(path)
##dataDict2 = sar.read2dseq(path2)

#data = dataDict['data']
##data2 = dataDict2['data']

#py.imshow(data[:,:,30])