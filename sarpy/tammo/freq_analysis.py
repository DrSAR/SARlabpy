# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 16:44:54 2012

Reproduction of Aerts et al.
Very messy! Convolution and multiplication
do not yield identical results.
Normalization of the FFT is an issue.
Time resolution of injection profiles influences magnitude
if the FFT data.

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



### Load AIF
aif_parker_time = np.load('aif_time.npy')
aif_parker_delta = aif_parker_time[1] - aif_parker_time[0]
p_aif_parker = np.load('p_aif.npy')
aif_parker = tl.aif_fit_function(aif_parker_time, p_aif_parker)
#aif_parker = np.load('parker_aif.npy')
###

Aerts_pars =[ 6.9e-5, 4.7e-2, 0.25, 3.4e-6, 1e-4, 2.3e-5, 1e-6 ]

trapez_parms = [ [12, 2e3, 6e3], [21, 2e3, 6e3], [7, 4e3, 6e3] ] # AERTS ET AL
trapez_parms = [ [12, 2e3, 6e3] ]
#~ trapez_parms = [ [16, 2, 15], [8, 4, 15], [4, 8, 15] ]


time_delta = 1
t_i = 0
t_f = 2000

time = np.arange(t_i, t_f, time_delta)
N = len(time)
freq_delta = (1/ (time_delta * N) )

freq = np.fft.fftfreq(N, d=time_delta)

# Caluclate array with several injection curves
injections = [ tl.Trapezoid(time, parms) for parms in trapez_parms ]

# Calculate Fourier Transform of injection profiles
injections_freq = np.fft.fft(injections)


# Calculate TDM from parameters
TDM_freq = tl.H_TDM(freq, Aerts_pars)

# Calculate FFT of the TDM
TDM = np.fft.ifft(TDM_freq) # DEPENDS ON SAMPLING FREQUENCY


# Convolve (multiply) TDM and injection
AIFs_freq = []
for injection_freq in injections_freq:
    AIFs_freq.append(TDM_freq * injection_freq)

AIFs = np.fft.ifft(AIFs_freq)


AIFs_c = []
for injection in injections:
    AIFs_c.append(tl.convolve_aif( TDM, injection, time ))


AIFs_c_np = []
for injection in injections:
    AIFs_c_np.append(np.convolve( TDM, injection))


plt.figure()
for AIF in AIFs_c:
    plt.plot(time[:N/(t_f/100)], AIF[:N/(t_f/100)])

for AIF in AIFs_c_np:
    plt.plot(time[:N/(t_f/100)], AIF[:N/(t_f/100)])


for AIF in AIFs:
    plt.plot(time[:N/(t_f/100)], AIF[:N/(t_f/100)])


plt.show()


raise SystemExit






################# Plot TDM in time and freq domain #####################
#plt.figure()
#plt.plot( np.log10(freq), TDM_freq, 'x')
#plt.title('TDM in frequency domain (xaxis: log10)')
#
#plt.figure()
#plt.plot( time[:500], TDM[:500], 'x')
#plt.title('TDM in time domain')
#
#plt.show()
########################################################################


########## Plot Injections in time and freq domain #####################
#plt.figure()
#plt.title('Injections curves in time domain')
#for injection in injections:
#	plt.plot(time[:], injection[:])
#
#
#plt.figure()
#plt.title('Injection curves in freq domain (xaxis: log10)')
#for injection_freq in injections_freq:
#	plt.plot( np.log10(freq) , ( abs(injection_freq)**2) )
#
#plt.show()
########################################################################

















def C_t_omega_AATH_0 (omega, K, v_p, v_e):
    return K/(1j*omega)

def C_t_omega (omega, K, v_p, v_e):
    return v_p + (K / ( (1j*omega) + (K/v_e) ) )
