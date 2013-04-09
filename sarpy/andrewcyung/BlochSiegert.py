# -*- coding: utf-8 -*-
# Copyright (C) 2012-2013 Andrew Yung and UBC 7T MRI lab
"""
Analysis routines for Bloch-Siegert B1 and flip angle mapping on Bruker data. 

Ref:  B1 mapping by Bloch-Siegert shift. Sacolick LI, Wiesinger F, Hancu I, 
Vogel MW.  Magn Reson Med. 2010 May;63(5):1315-22. 

Brief description of method:  off-resonance pulse applied after spin excitation, far away from Larmor frequency that it does not perturb magnetization, but strong enough to perturb direction of magnetic field.  Resultant phase shift is related to the peak B1 of the off-resonance pulse.  This B1 measurement can be used to map B1 and flip angle of other pulses (with knowledge of relative transmit power levels).

The technique requires at least two scans, repeated at positive and negative offset centered around the Larmor frequency.  The phase difference is proportional to B1, which is independent of B0 inhomogeneity up to the first order.  As of version 1.0, the method requires an additonal two scans with zero-power off-resonant pulses to determine the phase shift originating from switching the frequency synthesizer for the Bloch-Siegert pulse.

The analysis routines assume that the user is interested in using the peak B1 from the Bloch-Siegert experiment to determine the B1 or flip angle of a RF pulse of interest ("POI" from a different scan).  The user must specify a Bruker pulse shape directory to allow calculation of KBS constant.

Several function wrappers exist for calculation of flip angle maps of a pulse of interest (POI) which differ by their input type:
    
1) BS_flipangle_wrt_Bruker2dseq: file path of Bloch-Siegert expno and POI expno 

2) BS_flipangle_wrt_ScanObject: sarlabpy Scan objects

3) calc_BS_flipangle:  arrays for the individual Bloch-Siegert acquired scans

write_BS_flipangle_2dseq creates a 2dseq and visu_pars file for the calculated flip angle map, thus enabling direct visualization of the results in Paravision image display. 
"""

from __future__ import division

import numpy
import os
import pylab
import sarpy
import pdb
import scipy.integrate
import scipy.optimize
import scipy.fftpack
import sarpy.fmoosvi.getters as getters
import math

def BS_flipangle_wrt_Bruker2dseq(BS_2dseq_path, POI_2dseq_path, shapefile_path):

    BS_scan = BRUKER_classes.Scan(BS_2dseq_path)
    POI_scan = BRUKER_classes.Scan(POI_2dseq_path)
    
    flipanglemap = BS_flipangle_wrt_ScanObject(BS_scan, POI_scan, shapefile_path)

    return flipanglemap
    
def BS_flipangle_wrt_ScanObject(BS_scan, POI_scan, shapefile_path):

    try:
        dBpwr_BS = BS_scan.method.BSPulse[3]
        dBpwr_POI = POI_scan.method.ExcPulse[3]
        integralratio_POI = POI_scan.method.ExcPulse[10]
        width_BS = BS_scan.method.BSPulse[0]*1e-3
        width_POI = POI_scan.method.ExcPulse[0]*1e-3
        freqoffset = BS_scan.method.BSOffset
        BS_shape = BS_scan.method.BSPulse[1]
    except:
        print('something wrong happened with pulse parameter query')
        
    try:
        on_BSplus_phase = BS_scan.pdata[:,:,0]
        on_BSminus_phase = BS_scan.pdata[:,:,1]
        off_BSplus_phase = BS_scan.pdata[:,:,2]
        off_BSminus_phase = BS_scan.pdata[:,:,3]
    except:
        print('error in retrieving phase images')
        
    KBS = calc_KBS(freqoffset, width_BS, BS_shape, shapefile_path)

    B1peak_BS = calc_BS_B1peak(KBS, on_BSminus_phase, on_BSplus_phase, off_BSminus_phase, off_BSplus_phase)
    flipanglemap = calc_BS_flipangle(B1peak_BS, dBpwr_BS, dBpwr_POI, integralratio_POI, width_POI)
        
    return flipanglemap
    
def calc_BS_flipangle(B1peak_BS, dBpwr_BS, dBpwr_POI, integralratio_POI, width_POI):
    
    gamma = 267.513e6
    flipangle_POI = (gamma*B1peak_BS/10000) * (math.pow(10,(dBpwr_BS-dBpwr_POI)/20)) *\
                integralratio_POI*width_POI
                
    return flipangle_POI
    
def calc_BS_B1peak(KBS, on_BSminus_phase, on_BSplus_phase, off_BSminus_phase=0, off_BSplus_phase=0):
    
    if (off_BSminus_phase & off_BSplus_phase):
        offset = off_BSplus_phase - off_BSminus_phase
    else:
        offset = 0
    
    phase_diff = on_BSplus_phase - on_BSminus_phase + offset
    B1peak = numpy.sqrt(numpy.absolute(phase_diff)/(2*KBS))

    return B1peak
    
   

def calc_KBS(freqoffset, pulsewidth, pulseshape, shapefile_path, mode='gradientecho'):
    '''
calculates KBS in the units of radians/T^2.  Divide the result by 1e8 if units of   radians/Gauss^2.  See literature reference in the module docstring for details.

    '''
    filename = os.path.join(shapefile_path, pulseshape)
    pulse = sarpy.io.BRUKERIO.readRFshape(filename)
    dt = pulsewidth/pulse['NPOINTS']
    amp = numpy.asarray(pulse['amp'])
    amp = amp/max(amp)

    gamma = 267.513e6 # radians/(s*T)
    integrand = (gamma*amp)**2/(2*2*3.141592654*freqoffset)
    KBS = scipy.integrate.trapz(integrand, None, dt)

    return KBS
    
            
def h_phase_from_fid(scan_object):
    
    phase_data = numpy.angle(scipy.fftpack.fftshift(scipy.fftpack.fftn(scipy.fftpack.fftshift(scan_object.fid))))
    
    return phase_data
    
def h_mag_from_fid(scan_object):

    mag_data = numpy.abs(scipy.fftpack.fftshift(scipy.fftpack.fftn(scipy.fftpack.fftshift(scan_object.fid))))
    
    return mag_data
    
    
    
if __name__ == '__main__':
#    shapefile_path = 'opt/PV5.1/exp/stan/nmr/lists/wave'
    shapefile_path = os.path.expanduser(os.path.join('~','wave'))
#    KBS = calc_KBS(-4000,8e-3,'fermi.exc',shapefile_path)
    KBS = calc_KBS(4000,8e-3,'sech.exc',shapefile_path)
    print KBS/1e8
#    import doctest
#    doctest.testmod()






