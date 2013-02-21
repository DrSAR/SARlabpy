# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 11:36:19 2013

@author: fmoosvi
"""

import SARlabpy as sar
#import os
import pylab as py
import math
import numpy

reload(sar)

your_path = '/Users/firas/data/'
study_name = 'dBlochSiegert1.gP2' # could set this to be input by the system
default_reco = '/pdata/2/'

series_num = ['/36','/37','/46','/47'] # 36/37 are 0 power. 46/47 are high power
data_dict = []

procdirname_0 = ''.join([your_path,study_name,series_num[0],default_reco])
procdirname_1 = ''.join([your_path,study_name,series_num[1],default_reco])
procdirname_2 = ''.join([your_path,study_name,series_num[2],default_reco])    
procdirname_3 = ''.join([your_path,study_name,series_num[3],default_reco]) 
  
current_data_dict_0 = sar.read2dseq(procdirname_0,typecast=True)
current_data_dict_1 = sar.read2dseq(procdirname_1,typecast=True)
current_data_dict_2 = sar.read2dseq(procdirname_2,typecast=True)
current_data_dict_3 = sar.read2dseq(procdirname_3,typecast=True)

## Loadng in parameters as needed
    
n_read = current_data_dict_0['header']['d3proc']['IM_SIX']
n_phase = current_data_dict_0['header']['d3proc']['IM_SIY']
n_slice = current_data_dict_0['header']['d3proc']['IM_SIZ']

currslice = math.floor(n_slice/2);
 
TPQQ_BS = 20.24;
pulse_width = 2E-3
integral_ratio = 0.071941
TPQQ_POI = 5.00591
KBS = 71.16


# Calculate Offset and deal with phase wrap
offset = current_data_dict_1['data'] - current_data_dict_0['data']
offset_phase_unwrapped = numpy.unwrap(offset)

# Calculate BS phase difference, accounting for offset, deal with phase wrap
phase_diff = current_data_dict_3['data'] - current_data_dict_2['data'] + offset
phase_diff_unwrapped = numpy.unwrap(phase_diff)

# Calculate B1 peak
B1peak = numpy.sqrt(numpy.absolute(phase_diff_unwrapped)//(2*KBS))

# Calculate the flip angle for the pulse of interest 
gamma = 267.513e6

#alpha_BS = gamma*B1peak/10000*10^((TPQQ_BS-TPQQ_POI)/20)*IntegRatio*Pulsewidth*180/pi;
alpha_BS = (gamma*B1peak/10000) * (math.pow(10,TPQQ_BS-TPQQ_POI)/20) * integral_ratio*pulse_width*180/numpy.pi

#py.imshow(data[40,:,:],cmap='Greys')













