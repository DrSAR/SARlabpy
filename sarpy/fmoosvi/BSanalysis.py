#BSanalysis.py # -*- coding: utf-8 -*d-
"""
Created on Tue Feb 19 11:36:19 2013

@author: fmoosvi
"""

from __future__ import division
import numpy
import pylab
import math

import sarpy
import sarpy.fmoosvi.analysis
reload(sarpy)


dBloch = sarpy.Experiment('dBloch').studies[0].scans[45]
BSparams = sarpy.fmoosvi.analysis.temp_parse_BS(dBloch.method.BSPulse)
Excparams = sarpy.fmoosvi.analysis.temp_parse_BS(dBloch.method.ExcPulse)

# Constants and Parameters

gamma = 267.513e6
TPQQ_POI = 5.00591
TPQQ_BS = BSparams[3] #20.24
pulse_width = float(Excparams[0].replace('(','')) * 1E-3 #
integral_ratio = 0.071941
KBS = 71.16

high_f_power = sarpy.Experiment('dBloch').studies[0].scans[46] #  4000
low_f_power = sarpy.Experiment('dBloch').studies[0].scans[45]  # -4000

high_f_zero = sarpy.Experiment('dBloch').studies[0].scans[36]   # 4000
low_f_zero =sarpy.Experiment('dBloch').studies[0].scans[35]     #-4000

# Calculate Offset
offset = high_f_zero.pdata[1].data - low_f_zero.pdata[1].data

# Calculate BS phasedifference, accounting for offset
phase_diff = high_f_power.pdata[1].data - low_f_power.pdata[1].data

# Calculate B1 peak
B1peak = numpy.sqrt(numpy.absolute(phase_diff)/(2*KBS))

# Calculate Flip Angle for the pulse of interest
##alpha_BS = gamma*B1peak/10000*10^((TPQQ_BS-TPQQ_POI)/20)*IntegRatio*Pulsewidth*180/pi;
alpha_BS = (gamma*B1peak/10000) * (math.pow(10,(TPQQ_BS-TPQQ_POI)/20)) * integral_ratio*pulse_width*180/numpy.pi


##py.imshow(data[40,:,:],cmap='Greys')
#

pylab.imshow(alpha_BS[20,:,:])
pylab.colorbar()
pylab.show()




#your_path = '/Users/firas/data/'
#study_name = 'dBlochSiegert1.gP2' # could set this to be input by the system
#default_reco = '/pdata/2/'
#
#series_num = ['/36','/37','/46','/47'] # 36/37 are 0 power. 46/47 are high power
#data_dict = {}
#
#for i in range(0,len(series_num)):
#
#    procdirname = ''.join([your_path,study_name,series_num[i],default_reco])
#    scan_name = ''.join(['scan',str(i)])
#    current_data_dict = sarpy.read2dseq(procdirname,typecast=True)
#    
#    data_dict[scan_name] = current_data_dict
#
### Loadng in parameters as needed
#    
#n_read = current_data_dict['header']['d3proc']['IM_SIX']
#n_phase = current_data_dict['header']['d3proc']['IM_SIY']
#n_slice = current_data_dict['header']['d3proc']['IM_SIZ']
#
#currslice = math.floor(n_slice/2);
# 
#TPQQ_BS = float(data_dict['scan2']['header']['method']['BSPulse'][3]); # 20.24
#pulse_width = float(data_dict['scan2']['header']['method']['ExcPulse'][0])*1E-3 # actually it's 1e-3, AY had it 2E-3,
#integral_ratio = 0.071941
#TPQQ_POI = 5.00591
#KBS = 71.16
#
#
## Calculate Offset and deal with phase wrap
#offset = data_dict['scan0']['data'] - data_dict['scan1']['data']
##offset_phase_unwrapped = numpy.unwrap(offset,discont=2*numpy.pi)
#
## Calculate BS phase difference, accounting for offset, deal with phase wrap
#phase_diff = data_dict['scan2']['data'] - data_dict['scan3']['data'] + offset
##phase_diff_unwrapped = numpy.unwrap(phase_diff,discont=2*numpy.pi)
#
## Calculate B1 peak
#B1peak = numpy.sqrt(numpy.absolute(phase_diff)/(2*KBS))
#
## Calculate the flip angle for the pulse of interest
#gamma = 267.513e6
#
##alpha_BS = gamma*B1peak/10000*10^((TPQQ_BS-TPQQ_POI)/20)*IntegRatio*Pulsewidth*180/pi;
#alpha_BS = (gamma*B1peak/10000) * (math.pow(10,(TPQQ_BS-TPQQ_POI)/20)) * integral_ratio*pulse_width*180/numpy.pi
#
##py.imshow(data[40,:,:],cmap='Greys')
#
#py.imshow(alpha_BS[20,:,:])



