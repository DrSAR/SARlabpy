BSanalysis.py # -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 11:36:19 2013

@author: fmoosvi
"""

from __future__ import division
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
data_dict = {}

for i in range(0,len(series_num)):

    procdirname = ''.join([your_path,study_name,series_num[i],default_reco])
    scan_name = ''.join(['scan',str(i)])
    current_data_dict = sar.read2dseq(procdirname,typecast=True)
    
    data_dict[scan_name] = current_data_dict

## Loadng in parameters as needed
    
n_read = current_data_dict['header']['d3proc']['IM_SIX']
n_phase = current_data_dict['header']['d3proc']['IM_SIY']
n_slice = current_data_dict['header']['d3proc']['IM_SIZ']

currslice = math.floor(n_slice/2);
 
TPQQ_BS = float(data_dict['scan2']['header']['method']['BSPulse'][3]); # 20.24
pulse_width = float(data_dict['scan2']['header']['method']['ExcPulse'][0])*1E-3 # actually it's 1e-3, AY had it 2E-3,
integral_ratio = 0.071941
TPQQ_POI = 5.00591
KBS = 71.16


# Calculate Offset and deal with phase wrap
offset = data_dict['scan0']['data'] - data_dict['scan1']['data']
#offset_phase_unwrapped = numpy.unwrap(offset,discont=2*numpy.pi)

# Calculate BS phase difference, accounting for offset, deal with phase wrap
phase_diff = data_dict['scan2']['data'] - data_dict['scan3']['data'] + offset
#phase_diff_unwrapped = numpy.unwrap(phase_diff,discont=2*numpy.pi)

# Calculate B1 peak
B1peak = numpy.sqrt(numpy.absolute(phase_diff)/(2*KBS))

# Calculate the flip angle for the pulse of interest
gamma = 267.513e6

#alpha_BS = gamma*B1peak/10000*10^((TPQQ_BS-TPQQ_POI)/20)*IntegRatio*Pulsewidth*180/pi;
alpha_BS = (gamma*B1peak/10000) * (math.pow(10,(TPQQ_BS-TPQQ_POI)/20)) * integral_ratio*pulse_width*180/numpy.pi

#py.imshow(data[40,:,:],cmap='Greys')

py.imshow(alpha_BS[20,:,:])



