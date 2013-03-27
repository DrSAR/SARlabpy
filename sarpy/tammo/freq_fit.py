# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 13:03:03 2012

Fit third order transfer function with some injection protocol
to Jen's AIF.

@author: tammo
"""

from __future__ import print_function, division
import tammo_lib as tl
from scipy import optimize
import scipy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import norm
import matplotlib.mlab as mlab

#JEN approx: t=6s, rate=(1/60)*(ml/s) = , vol=100 mikro l


AIF_jen_time = np.load('aif_time.npy')

### Load Jen's AIF
AIF_jen = np.load('aif.npy')


### Load Parker AIF
#p_aif = np.load('p_aif.npy')
#AIF_jen = tl.aif_fit_function(AIF_jen_time, p_aif)


### Load parameters for TDM transfer function
p1 = np.load('p1_TDM.npy')
#Aerts_pars =[ 6.9e-5, 4.7e-2, 0.25, 3.4e-6, 1e-4, 2.3e-5, 1e-6 ]
#initial_pars = [.1,.1,.1,.1,.1,.1,.1]


## Set parameters for injection curve
volume = 20 #(microliter)
slope = .5


### set x_max to the AIF maximum
xmax = AIF_jen_time[AIF_jen.argmax()]#-(ymax/slope) # choose aif peak as inj end
print(max(AIF_jen_time))

### Create time array for TDM. Sampling important for FFT results. FFT symmetric, double time
time_TDM = np.arange(0, max(AIF_jen_time), .1)


### Interpolate AIF to time TDM
AIF_jen = np.interp(time_TDM, AIF_jen_time, AIF_jen)


### Compute injection curve
injection = tl.Trapezoid(time_TDM, [ volume, xmax, slope] )


### Plot injection and AIF
#plt.plot(time_TDM, injection)
#plt.plot(time_TDM, AIF_jen,'x')
#plt.show()


# play around (esp in freq space)
#TDM_AIF, freq = tl.AIF_from_TDM(injection, time_TDM, p1)
#plt.plot(time_TDM, np.fft.ifft(tl.H_TDM(freq, p1)))
#plt.show()
#raise SystemExit



### Perform fit of TDM to given AIF #######################################

def err_func (parms, args):
    time, injection, AIF_jen = args[:]
    out = tl.AIF_from_TDM(injection, time, parms) - AIF_jen
    print(sum(out)/len(time)) # error per data point
    return out

p1, success = optimize.leastsq(err_func, p1[:],\
    args = [time_TDM, injection, AIF_jen], maxfev=2000)

## Save Fit Parameters
np.save('p1_TDM', p1)
print('Fit Parameters:', p1)
np.save('p_TDM', p1)
###########################################################################

## Load Fit parameters (If fit is not performed):
#p1 = initial_pars

## Create AIF from Fit parameters
AIF = tl.AIF_from_TDM( injection , time_TDM,  p1  )


plt.plot(time_TDM, injection)
plt.plot(time_TDM, AIF_jen,'x')
plt.plot(time_TDM, AIF,'--')

plt.show()

#raise SystemExit